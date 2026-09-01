import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import build_multi_year_drug_vocabulary as vocab  # noqa: E402


class MultiYearVocabularyTests(unittest.TestCase):
    def _fixture(self, directory: str):
        root = Path(directory) / "annual"
        rows = {
            2023: [
                {"caseid": "01", "primaryid": "001", "prod_ai": " Aspirin 81mg ", "drugname": "Brand", "source_quarter": "2099Q4"},
                {"caseid": "01", "primaryid": "001", "prod_ai": "ASPIRIN", "drugname": "Brand", "source_quarter": "2099Q4"},
                {"caseid": "02", "primaryid": "002", "prod_ai": "", "drugname": "Amoxicillin / Clavulanate", "source_quarter": "2099Q4"},
                {"caseid": "02", "primaryid": "002", "prod_ai": "", "drugname": "", "source_quarter": "2099Q4"},
            ],
            2025: [
                {"caseid": "03", "primaryid": "003", "prod_ai": "aspirin", "drugname": "Different", "source_quarter": "2010Q1"},
                {"caseid": "04", "primaryid": "004", "prod_ai": "12", "drugname": "Ibuprofen 200mg", "source_quarter": "2010Q1"},
                {"caseid": "05", "primaryid": "005", "prod_ai": "", "drugname": "None!", "source_quarter": "2010Q1"},
            ],
        }
        assignments = []
        for year, values in rows.items():
            (root / str(year)).mkdir(parents=True)
            pd.DataFrame(values).to_parquet(root / str(year) / "resolved_drug_rows.parquet", index=False)
            for value in values:
                key = (value["caseid"], value["primaryid"])
                if not any((x["caseid"], x["selected_primaryid"]) == key for x in assignments):
                    assignments.append({"caseid": value["caseid"], "selected_primaryid": value["primaryid"],
                                        "assigned_report_year": year, "assignment_status": "included"})
        assignment_path = Path(directory) / "assignments.parquet"
        pd.DataFrame(assignments).to_parquet(assignment_path, index=False)
        return root, assignment_path

    def test_year_parsing(self):
        self.assertEqual(vocab.parse_years(["2023-2025"]), [2023, 2024, 2025])
        self.assertEqual(vocab.parse_years(["2018", "2020", "2023-2025"]), [2018, 2020, 2023, 2024, 2025])
        for values in ([], ["2025-2023"], ["bad"], ["2023", "2023"], ["2023-2024", "2024"]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                vocab.parse_years(values)

    def test_missing_prod_ai_column_treated_as_blank(self):
        """Stage 4 (2013-2022): prod_ai was added to the FAERS DRUG table
        after 2013, so 2013's resolved_drug_rows.parquet genuinely has no
        prod_ai column (not just null values). resolve_annual_drug_inputs
        must not hard-require it, and the vocabulary build must fall back to
        drugname exactly like _choose_raw_drug_name() already does elsewhere."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "annual"
            (root / "2013").mkdir(parents=True)
            pd.DataFrame([
                {"caseid": "01", "primaryid": "001", "drugname": "Ibuprofen 200mg"},
            ]).to_parquet(root / "2013" / "resolved_drug_rows.parquet", index=False)
            assignments = pd.DataFrame([
                {"caseid": "01", "selected_primaryid": "001", "assigned_report_year": 2013,
                 "assignment_status": "included"},
            ])
            assignment_path = Path(directory) / "assignments.parquet"
            assignments.to_parquet(assignment_path, index=False)

            inputs = vocab.resolve_annual_drug_inputs(root, [2013], "resolved_drug_rows.parquet")
            result, annual, metadata = vocab.build_multi_year_vocabulary(
                inputs, assignments_path=assignment_path, sample_raw_names_limit=2, batch_size=2,
            )
        self.assertIn("ibuprofen", set(result["normalized_drug_name"]))
        self.assertEqual(annual.set_index("year").loc[2013, "resolved_drug_rows"], 1)

    def test_input_resolution_failures(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(FileNotFoundError):
                vocab.resolve_annual_drug_inputs(root, [2024], "resolved.parquet")
            (root / "2024").mkdir()
            pd.DataFrame({"caseid": ["1"]}).to_parquet(root / "2024" / "resolved.parquet")
            with self.assertRaisesRegex(ValueError, "missing required"):
                vocab.resolve_annual_drug_inputs(root, [2024], "resolved.parquet")

    def test_multi_year_counts_blanks_combinations_and_analytical_year(self):
        with tempfile.TemporaryDirectory() as directory:
            root, assignments = self._fixture(directory)
            inputs = vocab.resolve_annual_drug_inputs(root, [2023, 2025], "resolved_drug_rows.parquet")
            result, annual, metadata = vocab.build_multi_year_vocabulary(
                inputs, assignments_path=assignments, sample_raw_names_limit=2, batch_size=2,
            )
        by_name = result.set_index("normalized_drug_name")
        self.assertEqual(by_name.loc["aspirin", "years_seen"], [2023, 2025])
        self.assertEqual(by_name.loc["aspirin", "first_seen_year"], 2023)
        self.assertEqual(by_name.loc["aspirin", "last_seen_year"], 2025)
        self.assertEqual(by_name.loc["aspirin", "year_count"], 2)
        self.assertEqual(by_name.loc["aspirin", "raw_row_count"], 3)
        self.assertEqual(by_name.loc["aspirin", "distinct_case_count"], 2)
        self.assertEqual(by_name.loc["aspirin", "sample_raw_names"], ["ASPIRIN", "Aspirin 81mg"])
        self.assertIn("amoxicillin|clavulanate", by_name.index)
        self.assertEqual(metadata["blank_normalized_name_rows"], 2)
        self.assertNotIn("none", by_name.index)
        self.assertEqual(annual.set_index("year").loc[2023, "blank_chosen_name_rows"], 1)
        self.assertEqual(set(annual["year"]), {2023, 2025})
        self.assertNotIn(2099, set().union(*result["years_seen"]))

    def test_output_is_independent_of_input_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root, assignments = self._fixture(directory)
            inputs = vocab.resolve_annual_drug_inputs(root, [2023, 2025], "resolved_drug_rows.parquet")
            first, _, _ = vocab.build_multi_year_vocabulary(inputs, assignments_path=assignments, batch_size=2)
            second, _, _ = vocab.build_multi_year_vocabulary(list(reversed(inputs)), assignments_path=assignments, batch_size=3)
            pd.testing.assert_frame_equal(first, second)

    def test_assignment_validation_rejects_unexpected_primaryid_and_caseid(self):
        with tempfile.TemporaryDirectory() as directory:
            root, assignments = self._fixture(directory)
            path = root / "2023" / "resolved_drug_rows.parquet"
            frame = pd.read_parquet(path); frame.loc[0, "primaryid"] = "999"; frame.to_parquet(path, index=False)
            inputs = vocab.resolve_annual_drug_inputs(root, [2023], "resolved_drug_rows.parquet")
            with self.assertRaisesRegex(ValueError, "unexpected selected PRIMARYID"):
                vocab.build_multi_year_vocabulary(inputs, assignments_path=assignments)
            frame.loc[0, "primaryid"] = "001"; frame.loc[0, "caseid"] = "wrong"; frame.to_parquet(path, index=False)
            with self.assertRaisesRegex(ValueError, "assignment mismatch"):
                vocab.build_multi_year_vocabulary(inputs, assignments_path=assignments)

    def test_reports_and_existing_output_protection(self):
        with tempfile.TemporaryDirectory() as directory:
            root, assignments = self._fixture(directory)
            inputs = vocab.resolve_annual_drug_inputs(root, [2023], "resolved_drug_rows.parquet")
            result, annual, metadata = vocab.build_multi_year_vocabulary(inputs, assignments_path=assignments, batch_size=2)
            output = Path(directory) / "output.parquet"; reports = Path(directory) / "reports"
            paths = vocab.write_vocabulary_outputs(
                result, annual, metadata, output=output, report_dir=reports,
                normalization_version="norm_v1", artifact_version="vocab_v1",
                sample_raw_names_limit=10, overwrite=False,
            )
            self.assertTrue(all(path.exists() for path in paths))
            collisions = pd.read_csv(reports / "ms5_normalization_collisions.csv")
            combos = pd.read_csv(reports / "ms5_combination_names.csv")
            self.assertTrue((collisions["distinct_raw_name_count"] > 1).all())
            self.assertEqual(combos["normalized_drug_name"].tolist(), ["amoxicillin|clavulanate"])
            self.assertFalse(list(Path(directory).rglob("*.part")))
            self.assertFalse(list(Path(directory).rglob("*.sqlite")))
            with self.assertRaises(FileExistsError):
                vocab.write_vocabulary_outputs(
                    result, annual, metadata, output=output, report_dir=reports,
                    normalization_version="norm_v1", artifact_version="vocab_v1",
                    sample_raw_names_limit=10, overwrite=False,
                )


if __name__ == "__main__":
    unittest.main()
