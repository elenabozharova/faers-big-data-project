import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import audit_missing_selected_reac as audit  # noqa: E402
from download_faers import Source, sha256_file  # noqa: E402


class MissingReacAuditTests(unittest.TestCase):
    def _source(self, year, quarter, filename):
        return Source(year, quarter, "FAERS", "ASCII", "https://example.test/a.zip", filename,
                      filename, "ASCII/DEMO.txt", "ASCII/DRUG.txt", "ASCII/REAC.txt")

    def _zip(self, directory, filename, rows):
        path = Path(directory) / filename
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("ASCII/REAC.txt", "PRIMARYID$CASEID$PT\n" + rows)
        return path

    def test_identifies_missing_and_preserves_leading_zero(self):
        assignments = pd.DataFrame([
            {"caseid": "01", "selected_primaryid": "001", "assigned_report_year": 2024,
             "assignment_status": "included", "selected_source_year": 2024, "selected_source_quarter": 1},
            {"caseid": "02", "selected_primaryid": "002", "assigned_report_year": 2024,
             "assignment_status": "included", "selected_source_year": 2024, "selected_source_quarter": 1},
        ])
        result = audit.identify_missing_selected_reac(
            assignments, pd.DataFrame({"primaryid": [" 001 "]}), 2024
        )
        self.assertEqual(result[["caseid", "selected_primaryid"]].values.tolist(), [["02", "002"]])

    def test_chunked_scan_finds_expected_other_and_same_case_alternative_without_modifying_zip(self):
        with tempfile.TemporaryDirectory() as directory:
            sources = [self._source(2023, 1, "a.zip"), self._source(2024, 1, "b.zip")]
            first = self._zip(directory, "a.zip", "002$02$other partition\n009$01$alternative\n")
            self._zip(directory, "b.zip", "001$$expected blank case\n003$99$irrelevant\n")
            before = sha256_file(first)
            occurrences, case_counts, chunks, _ = audit.scan_reac_sources(
                sources, Path(directory), {"001", "002", "009"}, {"01", "02"}, 1,
            )
            self.assertEqual(occurrences[("001", 2024, 1)]["rows"], 1)
            self.assertEqual(occurrences[("002", 2023, 1)]["rows"], 1)
            self.assertEqual(case_counts[("01", "009")], 1)
            self.assertGreater(chunks, 1)
            self.assertEqual(before, sha256_file(first))
            self.assertEqual(occurrences[("001", 2024, 1)]["caseids"], set())

    def test_deterministic_classifications_and_no_substitution(self):
        self.assertEqual(audit.classify_case(False, False, False), "confirmed_missing_in_raw_reac")
        self.assertEqual(audit.classify_case(True, False, False), "found_expected_partition_integration_error")
        self.assertEqual(audit.classify_case(False, True, False), "found_different_partition_provenance_error")
        self.assertEqual(
            audit.classify_case(False, False, True),
            "selected_version_missing_reac_nonselected_version_available",
        )
        self.assertEqual(audit.classify_case(True, False, False, True), "raw_caseid_conflict")

    def test_duplicate_or_disagreeing_selected_provenance_fails(self):
        missing = pd.DataFrame([{
            "caseid": "1", "selected_primaryid": "01", "selected_source_year": 2024,
            "selected_source_quarter": 1,
        }])
        versions = pd.DataFrame([{
            "caseid": "1", "primaryid": "01", "source_year": 2024,
            "source_quarter": 2, "is_selected": True,
        }])
        with self.assertRaisesRegex(ValueError, "source_quarter"):
            audit.validate_selected_provenance(missing, versions)
        duplicate = pd.concat([versions, versions], ignore_index=True)
        with self.assertRaises(ValueError):
            audit.validate_selected_provenance(missing, duplicate)

    def test_repository_report_path_is_relative_posix(self):
        value = audit.repository_path(ROOT / "reports" / "data_quality" / "audit.json")
        self.assertEqual(value, "reports/data_quality/audit.json")
        self.assertNotIn("\\", value)
        self.assertNotRegex(value, r"^[A-Za-z]:")


if __name__ == "__main__":
    unittest.main()
