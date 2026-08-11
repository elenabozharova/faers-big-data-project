"""End-to-end synthetic smoke test for the resumable MS6 production path."""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import rebuild_global_downstream as rebuild  # noqa: E402


def fixture(root: Path, year: int = 2099):
    boundary = root / "boundary"
    year_dir = boundary / str(year)
    year_dir.mkdir(parents=True)
    retained_rows = []
    drug_rows = []
    reac_rows = []
    assignment_rows = []
    for index in range(1, 11):
        caseid = str(index)
        primaryid = str(1000 + index)
        retained_rows.append({"caseid": caseid, "retained_primaryid": primaryid})
        assignment_rows.append({
            "caseid": caseid, "selected_primaryid": primaryid,
            "assigned_report_year": year, "assignment_status": "included",
            "selected_source_year": year, "selected_source_quarter": 1,
        })
        drug_rows.extend([
            {"caseid": caseid, "primaryid": primaryid, "role_cod": "PS", "drugname": "Aspirin", "prod_ai": ""},
            {"caseid": caseid, "primaryid": primaryid, "role_cod": "SS", "drugname": "Ibuprofen", "prod_ai": ""},
        ])
        reac_rows.append({"caseid": caseid, "primaryid": primaryid, "pt": "Headache"})
    # Duplicate input rows must deduplicate at case level.
    drug_rows.append(drug_rows[0].copy())
    reac_rows.append(reac_rows[0].copy())

    extra = [
        ("11", "1011", "Aspirin", "Nausea"),          # one mapped drug
        ("12", "1012", "Aspirin", None),              # DRUG but no REAC
        ("13", "1013", None, "Headache"),              # REAC but no DRUG
        ("14", "1014", "Unknown", "Headache"),       # unmatched
        ("15", "1015", "Ambiguous", "Headache"),     # ambiguous
        ("16", "1016", "Aspirin / Ibuprofen", "Headache"),
        ("17", "1017", "Aspirin / Unknown", "Headache"),
        ("18", "1018", "Ambiguous / Aspirin", "Headache"),
    ]
    for caseid, primaryid, drugname, reaction in extra:
        retained_rows.append({"caseid": caseid, "retained_primaryid": primaryid})
        assignment_rows.append({
            "caseid": caseid, "selected_primaryid": primaryid,
            "assigned_report_year": year, "assignment_status": "included",
            "selected_source_year": year, "selected_source_quarter": 1,
        })
        if drugname:
            drug_rows.append({"caseid": caseid, "primaryid": primaryid,
                              "role_cod": "PS", "drugname": drugname, "prod_ai": ""})
        if reaction:
            reac_rows.append({"caseid": caseid, "primaryid": primaryid, "pt": reaction})

    pd.DataFrame(retained_rows).to_parquet(year_dir / "retained_cases.parquet", index=False)
    pd.DataFrame(drug_rows).to_parquet(year_dir / "resolved_drug_rows.parquet", index=False)
    pd.DataFrame(reac_rows).to_parquet(year_dir / "resolved_reac_rows.parquet", index=False)
    assignments = pd.DataFrame(assignment_rows)
    mapping = pd.DataFrame([
        {"normalized_drug_name": "aspirin", "chembl_id": "CHEMBL25",
         "mapping_status": "mapped", "mapping_method": "exact_pref"},
        {"normalized_drug_name": "ibuprofen", "chembl_id": "CHEMBL521",
         "mapping_status": "mapped", "mapping_method": "exact_synonym"},
        {"normalized_drug_name": "unknown", "chembl_id": None,
         "mapping_status": "unmapped_no_match", "mapping_method": "unmatched"},
        {"normalized_drug_name": "ambiguous", "chembl_id": None,
         "mapping_status": "unmapped_ambiguous", "mapping_method": "exact_synonym"},
        {"normalized_drug_name": "aspirin|ibuprofen", "chembl_id": None,
         "mapping_status": "mapped_combination_components", "mapping_method": "combo_all_exact"},
        {"normalized_drug_name": "aspirin|unknown", "chembl_id": None,
         "mapping_status": "unmapped_combination_partial", "mapping_method": "combo_exact_components"},
        {"normalized_drug_name": "ambiguous|aspirin", "chembl_id": None,
         "mapping_status": "unmapped_combination_ambiguous", "mapping_method": "combo_exact_components"},
    ])
    mapping_path = root / "mapping.parquet"
    mapping.to_parquet(mapping_path, index=False)
    return boundary, assignments, mapping_path


class FullGlobalDownstreamTests(unittest.TestCase):
    def test_complete_boundary_to_threshold_outputs_and_safe_publish(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            boundary, assignments, mapping = fixture(root)
            output = root / "final"
            before = {name: rebuild.sha256_file(boundary / "2099" / name)
                      for name in rebuild.BOUNDARY_FILES}
            summary, performance = rebuild.rebuild_year(
                2099, boundary_root=boundary, output_root=output, mapping_path=mapping,
                assignments=assignments, pair_role_filter="all_pairs",
                max_rows_per_join_chunk=3, overwrite=False,
            )
            final = output / "2099"
            required = [*rebuild.BOUNDARY_FILES, "case_reactions.parquet", "case_drugs_long.parquet",
                        "case_drugs.parquet", "case_drug_pairs.parquet",
                        "case_drug_pair_role_audit.parquet", "case_drug_pairs_selected.parquet",
                        "case_pair_reactions.parquet", "drug_pair_reaction_counts.parquet",
                        "drug_pair_reaction_counts_min3.parquet", "drug_pair_reaction_counts_min5.parquet",
                        "drug_pair_reaction_counts_min10.parquet", "annual_summary.json"]
            self.assertTrue(all((final / name).is_file() for name in required))
            counts = pd.read_parquet(final / "drug_pair_reaction_counts.parquet")
            self.assertEqual(counts.loc[0, "n_cases"], 10)
            for threshold in (3, 5, 10):
                self.assertEqual(len(pd.read_parquet(final / f"drug_pair_reaction_counts_min{threshold}.parquet")), 1)
            long = pd.read_parquet(final / "case_drugs_long.parquet")
            excluded = long[long["mapping_status"] != "mapped"]
            self.assertTrue(excluded["mapped_drug_id"].isna().all())
            self.assertEqual(summary["case_usability"]["cases_with_at_least_two_single_id_mapped_drugs"], 10)
            self.assertGreater(len(performance), 6)
            self.assertEqual(before, {name: rebuild.sha256_file(boundary / "2099" / name)
                                      for name in rebuild.BOUNDARY_FILES})
            self.assertFalse(list(root.glob(".final-2099-*")))
            with self.assertRaises(FileExistsError):
                rebuild.rebuild_year(
                    2099, boundary_root=boundary, output_root=output, mapping_path=mapping,
                    assignments=assignments, pair_role_filter="all_pairs",
                    max_rows_per_join_chunk=3, overwrite=False,
                )


if __name__ == "__main__":
    unittest.main()
