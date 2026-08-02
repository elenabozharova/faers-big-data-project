import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import build_faers_chembl_mapping as mapper  # noqa: E402
from build_multi_year_drug_vocabulary import AnnualDrugInput  # noqa: E402
import process_faers_case_level as consumer  # noqa: E402


class ExactChemblMappingTests(unittest.TestCase):
    def _vocabulary(self):
        names = ["aspirin", "brand", "ambiguous", "unknown", "aspirin|brand",
                 "aspirin|unknown", "aspirin|ambiguous"]
        return pd.DataFrame([{ "normalized_drug_name": name, "first_seen_year": 2023,
            "last_seen_year": 2024, "years_seen": [2023, 2024], "year_count": 2,
            "raw_row_count": 2, "distinct_case_count": 1, "distinct_primaryid_count": 1,
            "distinct_raw_name_count": 1, "sample_raw_names": [name]} for name in names])

    def test_exact_precedence_ambiguity_and_combinations(self):
        preferred = {"aspirin": ["CHEMBL25"], "brand": ["CHEMBL_PREF"],
                     "ambiguous": ["CHEMBL1", "CHEMBL2"]}
        synonyms = {"aspirin": ["CHEMBL_OTHER"], "brand": ["CHEMBL_SYN"]}
        result = mapper.build_exact_mapping(self._vocabulary(), preferred, synonyms,
            chembl_version=36, normalization_version="norm_v1", vocabulary_artifact_version="vocab_v1",
            mapping_artifact_version="exact_v1", name_index_artifact_version="index_v1").set_index("normalized_drug_name")
        self.assertEqual(result.loc["aspirin", "chembl_id"], "CHEMBL25")
        self.assertEqual(result.loc["aspirin", "mapping_method"], "exact_pref")
        self.assertEqual(result.loc["ambiguous", "mapping_status"], "unmapped_ambiguous")
        self.assertEqual(result.loc["ambiguous", "mapping_candidates"], ["CHEMBL1", "CHEMBL2"])
        self.assertEqual(result.loc["unknown", "mapping_status"], "unmapped_no_match")
        self.assertEqual(result.loc["aspirin|brand", "mapping_status"], "mapped_combination_components")
        self.assertEqual(result.loc["aspirin|brand", "component_chembl_ids"], ["CHEMBL25", "CHEMBL_PREF"])
        self.assertEqual(result.loc["aspirin|unknown", "mapping_status"], "unmapped_combination_partial")
        self.assertEqual(result.loc["aspirin|ambiguous", "mapping_status"], "unmapped_combination_ambiguous")
        details = json.loads(result.loc["aspirin|ambiguous", "component_mapping_details"])
        self.assertEqual([item["component"] for item in details], ["aspirin", "ambiguous"])

    def test_name_index_preserves_duplicate_normalized_relationships(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "index.parquet"
            pd.DataFrame([
                {"normalized_chembl_name": "same", "chembl_id": "CHEMBL1", "name_type": "preferred", "chembl_version": 36, "normalization_version": "norm_v1"},
                {"normalized_chembl_name": "same", "chembl_id": "CHEMBL2", "name_type": "preferred", "chembl_version": 36, "normalization_version": "norm_v1"},
            ]).to_parquet(path, index=False)
            preferred, _ = mapper.load_name_index(path, chembl_version=36, normalization_version="norm_v1")
        self.assertEqual(preferred["same"], ["CHEMBL1", "CHEMBL2"])

    def test_consumer_uses_explicit_status_and_does_not_expose_combo_id(self):
        mapping = pd.DataFrame([
            {"normalized_drug_name": "aspirin", "chembl_id": "CHEMBL25", "mapping_status": "mapped", "mapping_method": "exact_pref"},
            {"normalized_drug_name": "ambiguous", "chembl_id": None, "mapping_status": "unmapped_ambiguous", "mapping_method": "exact_synonym"},
            {"normalized_drug_name": "aspirin|brand", "chembl_id": None, "mapping_status": "mapped_combination_components", "mapping_method": "combo_all_exact"},
        ])
        drugs = pd.DataFrame([
            {"caseid": "1", "primaryid": "11", "prod_ai": "aspirin", "drugname": "", "role_cod": "PS"},
            {"caseid": "2", "primaryid": "22", "prod_ai": "ambiguous", "drugname": "", "role_cod": "PS"},
            {"caseid": "3", "primaryid": "33", "prod_ai": "aspirin / brand", "drugname": "", "role_cod": "PS"},
        ])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mapping.parquet"; mapping.to_parquet(path, index=False)
            long, case_drugs, _, _ = consumer.build_case_drugs(drugs, str(path))
        self.assertEqual(long["mapping_status"].tolist(), ["mapped", "unmapped_ambiguous", "mapped_combination_components"])
        self.assertEqual(case_drugs["mapped_drug_id"].tolist(), ["CHEMBL25"])

    def test_legacy_ambiguous_null_id_is_not_no_match(self):
        mapping = pd.DataFrame([{"normalized": "ambiguous", "chembl_id": None,
                                 "is_ambiguous": True, "match_method": "exact_synonym"}])
        drugs = pd.DataFrame([{"caseid": "1", "primaryid": "1", "prod_ai": "ambiguous",
                               "drugname": "", "role_cod": "PS"}])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "old.parquet"; mapping.to_parquet(path, index=False)
            long, _, _, _ = consumer.build_case_drugs(drugs, str(path))
        self.assertEqual(long.loc[0, "mapping_status"], "unmapped_ambiguous")

    def test_coverage_has_union_all_row_and_separate_combinations(self):
        preferred = {"aspirin": ["CHEMBL25"], "brand": ["CHEMBL2"],
                     "ambiguous": ["CHEMBL1", "CHEMBL3"]}
        mapping = mapper.build_exact_mapping(self._vocabulary(), preferred, {}, chembl_version=36,
            normalization_version="norm_v1", vocabulary_artifact_version="vocab_v1",
            mapping_artifact_version="exact_v1", name_index_artifact_version="index_v1")
        rows = [{"caseid": str(index), "primaryid": str(index), "prod_ai": name,
                 "drugname": ""} for index, name in enumerate(self._vocabulary()["normalized_drug_name"], 1)]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "resolved.parquet"; pd.DataFrame(rows).to_parquet(path, index=False)
            coverage, status = mapper.build_coverage(mapping, [AnnualDrugInput(2023, path)], batch_size=2)
        all_row = coverage[coverage["year"].astype(str) == "ALL"].iloc[0]
        self.assertEqual(all_row["eligible_unique_normalized_names"], len(mapping))
        self.assertEqual(all_row["single_id_mapped_unique_names"], 2)
        self.assertEqual(all_row["mapped_combination_unique_names"], 1)
        self.assertEqual(all_row["mapped_combination_rows"], 1)
        self.assertEqual(all_row["ambiguous_unique_names"], 1)
        self.assertEqual(all_row["ambiguous_combination_unique_names"], 1)
        self.assertIn("ALL", set(status["year"].astype(str)))

    def test_versioned_outputs_are_portable_and_protected(self):
        preferred = {"aspirin": ["CHEMBL25"], "brand": ["CHEMBL2"],
                     "ambiguous": ["CHEMBL1", "CHEMBL3"]}
        mapping = mapper.build_exact_mapping(self._vocabulary(), preferred, {}, chembl_version=36,
            normalization_version="norm_v1", vocabulary_artifact_version="vocab_v1",
            mapping_artifact_version="exact_v1", name_index_artifact_version="index_v1")
        coverage = pd.DataFrame([{"year": "ALL", "eligible_drug_rows": 14,
            "invalid_or_ineligible_rows": 0, "eligible_unique_normalized_names": 7,
            "single_id_mapped_rows": 4, "single_id_unmapped_rows": 10,
            "single_id_mapped_unique_names": 2, "single_id_unmapped_unique_names": 5,
            "single_id_row_mapping_percentage": 28.57, "single_id_name_mapping_percentage": 28.57,
            "exact_pref_unique_names": 2, "exact_synonym_unique_names": 0,
            "mapped_combination_unique_names": 1, "mapped_combination_rows": 2,
            "ambiguous_unique_names": 2, "ambiguous_rows": 4,
            "unmatched_unique_names": 1, "unmatched_rows": 2,
            "partial_combination_unique_names": 1, "ambiguous_combination_unique_names": 1}])
        status = pd.DataFrame([{"year": "ALL", "mapping_status": "mapped",
            "mapping_method": "exact_pref", "unique_normalized_names": 2,
            "raw_row_count": 4, "distinct_case_count": 2}])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); vocabulary = root / "vocab.parquet"; index = root / "index.parquet"
            vocabulary.write_bytes(b"vocabulary"); index.write_bytes(b"index")
            output = root / "mapping.parquet"; reports = root / "reports"
            summary = mapper.write_outputs(mapping, coverage, status, output=output, report_dir=reports,
                vocabulary_path=vocabulary, index_path=index, requested_years=[2023, 2024],
                normalization_version="norm_v1", mapping_artifact_version="exact_v1", overwrite=False)
            self.assertFalse(summary["absolute_paths_persisted"])
            self.assertEqual(summary["exact_pref_count"] + summary["exact_synonym_count"], 2)
            self.assertEqual(summary["ambiguous_count"], 1)
            self.assertFalse(list(root.rglob("*.part")))
            with self.assertRaises(FileExistsError):
                mapper.write_outputs(mapping, coverage, status, output=output, report_dir=reports,
                    vocabulary_path=vocabulary, index_path=index, requested_years=[2023, 2024],
                    normalization_version="norm_v1", mapping_artifact_version="exact_v1", overwrite=False)


if __name__ == "__main__": unittest.main()
