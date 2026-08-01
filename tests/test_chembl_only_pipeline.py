"""Regression tests for the ChEMBL-only FAERS case-level pipeline."""

import inspect
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import map_faers_to_chembl  # noqa: E402,F401
import normalize_drug_names  # noqa: E402,F401
import process_faers_case_level as case_level  # noqa: E402
import process_faers_case_level_annual as annual  # noqa: E402


class ChemblOnlyPipelineTests(unittest.TestCase):
    def _mapping_path(self, directory: str) -> Path:
        mapping = pd.DataFrame(
            [
                {
                    "normalized": "aspirin",
                    "chembl_id": "CHEMBL25",
                    "is_ambiguous": False,
                    "match_method": "exact_pref",
                },
                {
                    "normalized": "ambiguous drug",
                    "chembl_id": "CHEMBL999",
                    "is_ambiguous": True,
                    "match_method": "exact_synonym",
                },
                {
                    "normalized": "fuzzy drug",
                    "chembl_id": "CHEMBL123",
                    "is_ambiguous": False,
                    "match_method": "fuzzy",
                },
            ]
        )
        path = Path(directory) / "mapping.parquet"
        mapping.to_parquet(path, index=False)
        return path

    def test_case_drugs_use_chembl_ids_and_preserve_unmapped_statuses(self):
        resolved = pd.DataFrame(
            [
                {"caseid": "1", "primaryid": "11", "role_cod": "PS", "drugname": "Aspirin", "prod_ai": ""},
                {"caseid": "1", "primaryid": "11", "role_cod": "C", "drugname": "Unknown medicine", "prod_ai": ""},
                {"caseid": "2", "primaryid": "22", "role_cod": "SS", "drugname": "Ambiguous drug", "prod_ai": ""},
                {"caseid": "3", "primaryid": "33", "role_cod": "PS", "drugname": "Fuzzy drug", "prod_ai": ""},
                {"caseid": "4", "primaryid": "44", "role_cod": "PS", "drugname": "", "prod_ai": ""},
            ]
        )

        with tempfile.TemporaryDirectory() as directory:
            mapping_path = self._mapping_path(directory)
            long_df, case_drugs, diagnostics, _ = case_level.build_case_drugs(
                resolved, str(mapping_path)
            )

        self.assertEqual(case_drugs["mapped_drug_id"].tolist(), ["CHEMBL25"])
        mode_field = "drug_" + "id_mode"
        self.assertNotIn(mode_field, long_df.columns)
        self.assertNotIn(mode_field, case_drugs.columns)
        self.assertNotIn("n_removed_outside_drug_universe", diagnostics)
        self.assertEqual(
            long_df["mapping_status"].tolist(),
            [
                "mapped",
                "unmapped_no_match",
                "unmapped_ambiguous",
                "excluded_fuzzy",
                "unmapped_blank_name",
            ],
        )
        case_level.run_case_drugs_validation(resolved, long_df, case_drugs)

    def test_canonical_chembl_pairs_are_deterministic_and_not_self_pairs(self):
        case_drugs = pd.DataFrame(
            [
                {"caseid": "1", "mapped_drug_id": "CHEMBL2"},
                {"caseid": "1", "mapped_drug_id": "CHEMBL10"},
                {"caseid": "1", "mapped_drug_id": "CHEMBL1"},
                {"caseid": "1", "mapped_drug_id": "CHEMBL2"},
                {"caseid": "2", "mapped_drug_id": "CHEMBL7"},
            ]
        )

        pairs, _, _ = case_level.build_case_drug_pairs(case_drugs)

        self.assertEqual(
            list(pairs[["drug_a", "drug_b"]].itertuples(index=False, name=None)),
            [
                ("CHEMBL1", "CHEMBL10"),
                ("CHEMBL1", "CHEMBL2"),
                ("CHEMBL10", "CHEMBL2"),
            ],
        )
        self.assertFalse((pairs["drug_a"] == pairs["drug_b"]).any())
        case_level.run_case_drug_pairs_validation(case_drugs, pairs)

    def test_no_graph_index_mode_or_loader_remains(self):
        source = inspect.getsource(case_level).lower()
        mode_field = "drug_" + "id_mode"
        graph_index_name = "drug_" + "to_idx"
        self.assertNotIn(mode_field, source)
        self.assertNotIn(graph_index_name, source)
        self.assertFalse(hasattr(case_level, "load_drug_universe"))
        self.assertNotIn(mode_field, inspect.signature(annual.process_year).parameters)

    def test_pair_reaction_outputs_do_not_require_graph_indices(self):
        selected_pairs = pd.DataFrame(
            [
                {
                    "caseid": "1",
                    "drug_a": "CHEMBL1",
                    "drug_b": "CHEMBL2",
                    "n_mapped_drugs_in_case": 2,
                    "pair_role_filter": "all_pairs",
                    "is_cc_pair": False,
                    "has_ps": True,
                    "has_ss": False,
                    "has_i": False,
                    "has_suspect_or_interacting": True,
                }
            ]
        )
        reactions = pd.DataFrame([{"caseid": "1", "reaction_pt": "headache"}])

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "pair_reactions.parquet"
            total, counts, _, validation = case_level.build_case_pair_reactions_chunked(
                selected_pairs,
                reactions,
                str(output),
                "all_pairs",
                max_rows_per_chunk=10,
            )

        self.assertEqual(total, 1)
        self.assertEqual(counts.loc[0, "n_cases"], 1)
        self.assertEqual(counts.loc[0, "drug_a"], "CHEMBL1")
        self.assertNotIn("drug_" + "id_mode", counts.columns)
        case_level.run_case_pair_reactions_validation(validation, "all_pairs")

    def test_active_cli_help_has_no_identifier_mode(self):
        for script in ("process_faers_case_level.py", "process_faers_case_level_annual.py"):
            result = subprocess.run(
                [sys.executable, str(SRC / script), "--help"],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertNotIn("--drug-" + "id-mode", result.stdout)

    def test_tiny_faers_zip_runs_through_distinct_case_support_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            zip_path = directory_path / "faers.zip"
            with zipfile.ZipFile(zip_path, "w") as archive:
                archive.writestr(
                    "ASCII/DEMO.txt",
                    "PRIMARYID$CASEID$CASEVERSION$FDA_DT\n"
                    "100$10$1$20230101\n"
                    "101$10$2$20230401\n"
                    "200$20$1$20230201\n"
                    "300$30$1$20230301\n",
                )
                archive.writestr(
                    "ASCII/DRUG.txt",
                    "PRIMARYID$CASEID$ROLE_COD$DRUGNAME$PROD_AI\n"
                    "100$10$PS$Superseded drug$\n"
                    "101$10$PS$Aspirin$\n"
                    "101$10$C$ASPIRIN$\n"
                    "101$10$SS$Ibuprofen$\n"
                    "200$20$PS$Aspirin$\n"
                    "300$30$PS$Aspirin$\n"
                    "300$30$C$Ibuprofen$\n",
                )
                archive.writestr(
                    "ASCII/REAC.txt",
                    "PRIMARYID$CASEID$PT\n"
                    "100$10$Superseded reaction\n"
                    "101$10$Headache\n"
                    "101$10$  headache  \n"
                    "101$10$Nausea\n"
                    "300$30$Headache\n"
                    "300$30$HEADACHE\n",
                )

            mapping_path = directory_path / "mapping.parquet"
            pd.DataFrame(
                [
                    {
                        "normalized": "aspirin",
                        "chembl_id": "CHEMBL25",
                        "is_ambiguous": False,
                        "match_method": "exact_pref",
                    },
                    {
                        "normalized": "ibuprofen",
                        "chembl_id": "CHEMBL521",
                        "is_ambiguous": False,
                        "match_method": "exact_pref",
                    },
                ]
            ).to_parquet(mapping_path, index=False)

            demo = case_level.load_demo_table(str(zip_path), "", "ASCII/DEMO.txt")
            drug = case_level.load_drug_table(str(zip_path), "", "ASCII/DRUG.txt")
            reac = case_level.load_reac_table(str(zip_path), "", "ASCII/REAC.txt")
            retained = case_level.resolve_case_versions(demo)
            resolved_drug = case_level.filter_to_retained_primaryids(drug, retained, "DRUG")
            resolved_reac = case_level.filter_to_retained_primaryids(reac, retained, "REAC")

            retained_case_10 = retained.loc[retained["caseid"] == "10", "retained_primaryid"]
            self.assertEqual(retained_case_10.tolist(), ["101"])
            self.assertNotIn("100", resolved_drug["primaryid"].tolist())
            self.assertNotIn("100", resolved_reac["primaryid"].tolist())

            reactions, _ = case_level.build_case_reactions(resolved_reac)
            self.assertEqual(
                set(reactions[["caseid", "reaction_pt"]].itertuples(index=False, name=None)),
                {("10", "headache"), ("10", "nausea"), ("30", "headache")},
            )

            _, case_drugs, _, _ = case_level.build_case_drugs(
                resolved_drug, str(mapping_path)
            )
            case_10_drugs = case_drugs[case_drugs["caseid"] == "10"]
            self.assertEqual(set(case_10_drugs["mapped_drug_id"]), {"CHEMBL25", "CHEMBL521"})
            aspirin_rows = case_10_drugs.loc[
                case_10_drugs["mapped_drug_id"] == "CHEMBL25", "n_raw_rows_for_case_drug"
            ]
            self.assertEqual(aspirin_rows.tolist(), [2])

            pairs, _, _ = case_level.build_case_drug_pairs(case_drugs)
            self.assertEqual(
                set(pairs[["caseid", "drug_a", "drug_b"]].itertuples(index=False, name=None)),
                {
                    ("10", "CHEMBL25", "CHEMBL521"),
                    ("30", "CHEMBL25", "CHEMBL521"),
                },
            )

            role_audit, _ = case_level.build_case_drug_pair_role_audit(pairs, case_drugs)
            selected, _ = case_level.apply_pair_role_filter(pairs, role_audit, "all_pairs")
            output_path = directory_path / "case_pair_reactions.parquet"
            total, counts, _, validation = case_level.build_case_pair_reactions_chunked(
                selected,
                reactions,
                str(output_path),
                "all_pairs",
                max_rows_per_chunk=2,
            )

            self.assertEqual(total, 3)
            supports = counts.set_index("reaction_pt")["n_cases"].to_dict()
            self.assertEqual(supports, {"headache": 2, "nausea": 1})
            case_level.run_case_pair_reactions_validation(validation, "all_pairs")

    def test_valid_chembl_name_cache_is_used_without_a_connection(self):
        with tempfile.TemporaryDirectory() as directory:
            cache_path = Path(directory) / "chembl_name_index.parquet"
            pd.DataFrame(
                [
                    {
                        "chembl_id": "CHEMBL25",
                        "pref_name": "Aspirin",
                        "molecule_type": "Small molecule",
                        "synonym": "Acetylsalicylic acid",
                        "syn_type": "INN",
                    }
                ]
            ).to_parquet(cache_path, index=False)

            pref, syn = map_faers_to_chembl.build_chembl_name_index(None, directory)

        self.assertEqual(pref["aspirin"][0]["chembl_id"], "CHEMBL25")
        self.assertEqual(syn["acetylsalicylic acid"][0]["chembl_id"], "CHEMBL25")

    def test_missing_or_forced_chembl_name_cache_requires_a_connection(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(RuntimeError, "database connection is required"):
                map_faers_to_chembl.build_chembl_name_index(None, directory)

            cache_path = Path(directory) / "chembl_name_index.parquet"
            pd.DataFrame(
                columns=["chembl_id", "pref_name", "molecule_type", "synonym", "syn_type"]
            ).to_parquet(cache_path, index=False)
            with self.assertRaisesRegex(RuntimeError, "database connection is required"):
                map_faers_to_chembl.build_chembl_name_index(
                    None, directory, force_rebuild=True
                )

    def test_invalid_chembl_name_cache_fails_with_a_clear_error(self):
        with tempfile.TemporaryDirectory() as directory:
            cache_path = Path(directory) / "chembl_name_index.parquet"
            pd.DataFrame([{"chembl_id": "CHEMBL25", "pref_name": "Aspirin"}]).to_parquet(
                cache_path, index=False
            )
            with self.assertRaisesRegex(ValueError, "missing required column"):
                map_faers_to_chembl.build_chembl_name_index(None, directory)

            pd.DataFrame(
                [
                    {
                        "chembl_id": None,
                        "pref_name": "Aspirin",
                        "synonym": None,
                        "syn_type": None,
                    }
                ]
            ).to_parquet(cache_path, index=False)
            with self.assertRaisesRegex(ValueError, "missing or blank chembl_id"):
                map_faers_to_chembl.build_chembl_name_index(None, directory)


if __name__ == "__main__":
    unittest.main()
