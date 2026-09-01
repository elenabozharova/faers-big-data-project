"""Focused tests for the MS4B global-assignment boundary."""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import integrate_global_case_assignments as ms4b  # noqa: E402
import process_faers_case_level as pipeline  # noqa: E402


def assignments() -> pd.DataFrame:
    return pd.DataFrame([
        {"caseid": "A", "selected_primaryid": "1001", "assigned_report_year": 2023,
         "assignment_status": "included", "selected_source_year": 2023, "selected_source_quarter": 1,
         "selected_caseversion": 1, "assignment_method": "init_fda_dt",
         "has_later_year_followup_in_loaded_data": True, "later_followup_primaryids": ["2001"]},
        {"caseid": "B", "selected_primaryid": "0010", "assigned_report_year": 2024,
         "assignment_status": "included", "selected_source_year": 2023, "selected_source_quarter": 2,
         "selected_caseversion": 2, "assignment_method": "init_fda_dt",
         "has_later_year_followup_in_loaded_data": True, "later_followup_primaryids": []},
        {"caseid": "X", "selected_primaryid": "9000", "assigned_report_year": 2024,
         "assignment_status": "assigned_outside_study_window", "selected_source_year": 2024,
         "selected_source_quarter": 1, "has_later_year_followup_in_loaded_data": False},
    ])


class Ms4bTests(unittest.TestCase):
    def test_identifier_normalization_preserves_leading_zero(self):
        result = ms4b.normalize_identifier(pd.Series([" 0010 ", "", None]))
        self.assertEqual(result.iloc[0], "0010")
        self.assertTrue(result.iloc[1:].isna().all())
        self.assertNotIn(".0", result.iloc[0])

    def test_schema_loading_and_duplicate_selected_primaryid_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "assignments.parquet"
            assignments().to_parquet(path, index=False)
            loaded = ms4b.load_global_assignments(path)
            self.assertEqual(len(loaded), 3)
            bad = assignments()
            bad.loc[1, "selected_primaryid"] = "1001"
            bad.to_parquet(path, index=False)
            with self.assertRaisesRegex(ValueError, "multiple CASEIDs"):
                ms4b.load_global_assignments(path)

    def test_inclusion_ignores_audit_flags_and_adapts_contract(self):
        selected = ms4b.selected_cases_for_year(assignments(), 2023)
        self.assertEqual(selected["caseid"].tolist(), ["A"])
        retained = ms4b.adapt_retained_cases(selected)
        self.assertEqual(retained["retained_primaryid"].tolist(), ["1001"])

    def test_many_to_one_filter_attaches_authoritative_caseid_and_counts_blank(self):
        selected = ms4b.selected_cases_for_year(assignments(), 2023)
        raw = pd.DataFrame([
            {"primaryid": "1001", "caseid": "A", "value": 1},
            {"primaryid": "1001", "caseid": "", "value": 2},
            {"primaryid": "2001", "caseid": "A", "value": 3},
        ])
        filtered, stats = ms4b.filter_selected_rows(raw, selected, "DRUG")
        self.assertEqual(filtered["primaryid"].tolist(), ["1001", "1001"])
        self.assertEqual(filtered["caseid"].tolist(), ["A", "A"])
        self.assertEqual(stats["raw_caseid_blank_count"], 1)

    def test_conflicting_raw_caseid_fails(self):
        selected = ms4b.selected_cases_for_year(assignments(), 2023)
        with self.assertRaisesRegex(ValueError, "conflict"):
            ms4b.filter_selected_rows(
                pd.DataFrame([{"primaryid": "1001", "caseid": "OTHER"}]), selected, "REAC"
            )

    def test_two_year_routing_missing_tables_and_later_followup_exclusion(self):
        tables = {
            (2023, 1, "drug"): pd.DataFrame([{"primaryid": "1001", "caseid": "A", "drugname": "d"},
                                                    {"primaryid": "2001", "caseid": "A", "drugname": "later"}]),
            (2023, 1, "reac"): pd.DataFrame([{"primaryid": "1001", "caseid": "A", "pt": "p"},
                                                    {"primaryid": "2001", "caseid": "A", "pt": "later"}]),
            (2023, 2, "drug"): pd.DataFrame([{"primaryid": "0010", "caseid": "B", "drugname": "d"}]),
            (2023, 2, "reac"): pd.DataFrame(columns=["primaryid", "caseid", "pt"]),
        }
        loader = lambda sy, sq, name: tables[(sy, sq, name)].copy()
        retained23, drug23, reac23, report23 = ms4b.integrate_year(assignments(), 2023, loader)
        retained24, drug24, reac24, report24 = ms4b.integrate_year(assignments(), 2024, loader)
        self.assertEqual(retained23["caseid"].tolist(), ["A"])
        self.assertEqual(retained24["caseid"].tolist(), ["B"])
        self.assertEqual(set(drug23["primaryid"]), {"1001"})
        self.assertEqual(set(reac23["primaryid"]), {"1001"})
        self.assertNotIn("2001", set(drug23["primaryid"]) | set(reac23["primaryid"]))
        self.assertEqual(set(drug24["primaryid"]), {"0010"})
        self.assertTrue(reac24.empty)
        self.assertEqual(report24["selected_cases_missing_reac"], 1)
        self.assertEqual(report24["selected_cases_missing_both"], 0)
        self.assertEqual(report23["later_followup_primaryids_excluded"], 1)

    def test_self_referential_later_followup_primaryid_is_not_a_leak(self):
        """Stage 4 (2013-2022): a rare, real FAERS anomaly (~340/16.4M included
        cases) republishes the identical (CASEID, CASEVERSION) row -- hence the
        identical PRIMARYID -- in a later year's export with no version
        chronology conflict. The winner's own PRIMARYID then legitimately
        appears in its own later_followup_primaryids list; that must not be
        treated as a leaked followup reaching downstream."""
        rows = assignments()
        rows.loc[rows["caseid"] == "A", "later_followup_primaryids"] = pd.Series(
            [["1001"]], index=rows.index[rows["caseid"] == "A"]
        )
        tables = {
            (2023, 1, "drug"): pd.DataFrame([{"primaryid": "1001", "caseid": "A", "drugname": "d"}]),
            (2023, 1, "reac"): pd.DataFrame([{"primaryid": "1001", "caseid": "A", "pt": "p"}]),
            (2023, 2, "drug"): pd.DataFrame([{"primaryid": "0010", "caseid": "B", "drugname": "d"}]),
            (2023, 2, "reac"): pd.DataFrame(columns=["primaryid", "caseid", "pt"]),
        }
        loader = lambda sy, sq, name: tables[(sy, sq, name)].copy()
        retained, drug, reac, report = ms4b.integrate_year(rows, 2023, loader)
        self.assertEqual(retained["caseid"].tolist(), ["A"])
        self.assertEqual(set(drug["primaryid"]), {"1001"})
        self.assertEqual(report["later_followup_primaryids_excluded"], 0)

    def test_selected_case_missing_both_is_retained(self):
        empty = lambda *_: pd.DataFrame(columns=["primaryid", "caseid"])
        retained, drug, reac, report = ms4b.integrate_year(assignments(), 2023, empty)
        self.assertEqual(len(retained), 1)
        self.assertTrue(drug.empty and reac.empty)
        self.assertEqual(report["selected_cases_missing_both"], 1)

    def test_existing_downstream_reaction_result_is_unchanged_for_identical_input(self):
        resolved = pd.DataFrame([{"caseid": "A", "primaryid": "1001", "pt": " Headache "}])
        before, before_diag = pipeline.build_case_reactions(resolved.copy())
        after, after_diag = pipeline.build_case_reactions(resolved.copy())
        pd.testing.assert_frame_equal(before, after)
        self.assertEqual(before_diag, after_diag)


if __name__ == "__main__":
    unittest.main()
