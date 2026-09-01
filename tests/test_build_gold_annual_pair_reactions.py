import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import build_gold_annual_pair_reactions as gold  # noqa: E402


class GoldDatasetTests(unittest.TestCase):
    def test_parse_year_ranges(self):
        self.assertEqual(gold.parse_years(["2023-2025", "2024"]), [2023, 2024, 2025])

    def test_builds_strict_partitioned_contract(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source_dir = root / "input" / "2024"
            source_dir.mkdir(parents=True)
            pd.DataFrame([
                {"year": 2024, "drug_a": "CHEMBL1", "drug_b": "CHEMBL2", "reaction_pt": "a", "n_cases": 1, "pair_role_filter": "all_pairs"},
                {"year": 2024, "drug_a": "CHEMBL1", "drug_b": "CHEMBL3", "reaction_pt": "b", "n_cases": 2, "pair_role_filter": "all_pairs"},
            ]).to_parquet(source_dir / "drug_pair_reaction_counts.parquet", index=False)
            assignments = root / "assignments.parquet"
            pd.DataFrame([
                {"caseid": "1", "assigned_report_year": 2024, "assignment_status": "included"},
            ]).to_parquet(assignments, index=False)
            output = root / "gold"
            report = root / "report.json"
            result = gold.build_gold([2024], root / "input", output, assignments, report, 1, False)
            frame = pd.read_parquet(output / "year=2024" / "part-00000.parquet")
            self.assertEqual(list(frame.columns), gold.GOLD_COLUMNS)
            self.assertEqual(result["total_rows"], 2)
            self.assertTrue((output / "_SUCCESS").exists())
            self.assertTrue(report.exists())

    def test_additive_year_run_preserves_other_years_partitions(self):
        """Stage 4 (2013-2022): building Gold for new years must never delete
        or touch already-published years' partitions. A prior implementation
        staged and swapped the *entire* output_root, which would have wiped
        out every other year's partition the moment this ran for a
        non-overlapping year range (e.g. adding 2013-2022 after 2023-2025)."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            assignments_rows = []
            for year in (2024, 2013):
                source_dir = root / "input" / str(year)
                source_dir.mkdir(parents=True)
                pd.DataFrame([
                    {"year": year, "drug_a": "CHEMBL1", "drug_b": "CHEMBL2", "reaction_pt": "a",
                     "n_cases": 1, "pair_role_filter": "all_pairs"},
                ]).to_parquet(source_dir / "drug_pair_reaction_counts.parquet", index=False)
                assignments_rows.append({"caseid": str(year), "assigned_report_year": year,
                                          "assignment_status": "included"})
            assignments = root / "assignments.parquet"
            pd.DataFrame(assignments_rows).to_parquet(assignments, index=False)
            output = root / "gold"
            report = root / "report.json"

            gold.build_gold([2024], root / "input", output, assignments, report, 1, False)
            self.assertTrue((output / "year=2024" / "part-00000.parquet").exists())

            gold.build_gold([2013], root / "input", output, assignments, report, 1, False)
            self.assertTrue((output / "year=2024" / "part-00000.parquet").exists(),
                             "adding 2013 must not delete the existing 2024 partition")
            self.assertTrue((output / "year=2013" / "part-00000.parquet").exists())
            with (output / "_metadata.json").open() as handle:
                metadata = json.load(handle)
            self.assertEqual(metadata["years"], [2013, 2024])
            self.assertEqual(metadata["total_rows"], 2)

    def test_rejects_invalid_pair_order(self):
        frame = pd.DataFrame([{
            "year": 2024, "drug_a": "CHEMBL2", "drug_b": "CHEMBL1",
            "reaction_pt": "x", "n_cases": 1,
        }])
        with self.assertRaisesRegex(ValueError, "strictly ordered"):
            gold._validate_batch(frame, 2024)

    def test_rejects_caseid_assigned_to_multiple_years(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "assignments.parquet"
            pd.DataFrame([
                {"caseid": "1", "assigned_report_year": 2023, "assignment_status": "included"},
                {"caseid": "1", "assigned_report_year": 2024, "assignment_status": "included"},
            ]).to_parquet(path, index=False)
            with self.assertRaisesRegex(ValueError, "more than one assigned year"):
                gold.validate_global_assignments(path, [2023, 2024])


if __name__ == "__main__":
    unittest.main()
