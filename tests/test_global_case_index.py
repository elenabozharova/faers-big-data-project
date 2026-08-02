"""Synthetic tests for optimized global CASEID resolution."""

import csv
import json
import sqlite3
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import build_global_case_index as global_index
from download_faers import MANIFEST_COLUMNS


class GlobalCaseIndexTests(unittest.TestCase):
    def observation(
        self,
        row_id=1,
        caseid="1",
        primaryid="10",
        version="1",
        source_year=2024,
        source_quarter=1,
        fda_date="2024-01-01",
        initial_date="2024-01-01",
    ):
        version_number = global_index.parse_integer(version)
        return global_index.Observation(
            id=row_id,
            caseid=caseid,
            primaryid=primaryid,
            caseversion_raw=version,
            caseversion_numeric=version_number,
            caseversion_status=(
                "missing"
                if not version
                else "valid" if version_number is not None else "invalid"
            ),
            fda_dt_raw=fda_date.replace("-", "") if fda_date else "",
            fda_dt_parsed=fda_date,
            fda_dt_status="valid" if fda_date else "missing",
            init_fda_dt_raw=initial_date.replace("-", "") if initial_date else "",
            init_fda_dt_parsed=initial_date,
            init_fda_dt_status="valid" if initial_date else "missing",
            source_year=source_year,
            source_quarter=source_quarter,
            source_zip="data/raw/faers/q.zip",
            source_inner_filename="ASCII/DEMO.txt",
            source_row_number=row_id,
        )

    def test_approved_assignment_eligibility_and_winner_rules(self):
        rows = [
            self.observation(version="1", source_quarter=1, row_id=1),
            self.observation(primaryid="11", version="2", source_quarter=4, row_id=2),
            self.observation(primaryid="12", version="3", source_year=2025, row_id=3),
        ]
        assignment, winner, eligible, later = global_index.resolve_case(
            rows, set(), [2024, 2025]
        )
        self.assertEqual(assignment["assigned_report_year"], 2024)
        self.assertEqual(winner.primaryid, "11")
        self.assertEqual(eligible, {1, 2})
        self.assertEqual(later, {3})
        self.assertTrue(assignment["has_later_followup_in_loaded_data"])

    def test_valid_version_beats_missing_and_invalid(self):
        rows = [
            self.observation(primaryid="100", version="", row_id=1),
            self.observation(primaryid="200", version="bad", row_id=2),
            self.observation(primaryid="3", version="1", row_id=3),
        ]
        assignment, winner, _, _ = global_index.resolve_case(rows, set(), [2024])
        self.assertEqual(winner.primaryid, "3")
        self.assertEqual(assignment["missing_caseversion_count"], 1)
        self.assertEqual(assignment["invalid_caseversion_count"], 1)

    def test_same_year_version_is_additional_but_not_later_year(self):
        rows = [
            self.observation(version="1", source_quarter=1, row_id=1),
            self.observation(primaryid="11", version="2", source_quarter=4, row_id=2),
        ]
        assignment, _, _, _ = global_index.resolve_case(rows, set(), [2024])
        self.assertTrue(assignment["has_additional_version_in_loaded_data"])
        self.assertFalse(assignment["has_later_year_followup_in_loaded_data"])
        self.assertFalse(assignment["has_later_followup_in_loaded_data"])

    def test_later_year_version_is_additional_and_cross_year(self):
        rows = [
            self.observation(version="1", source_year=2024, row_id=1),
            self.observation(primaryid="11", version="2", source_year=2025, row_id=2),
        ]
        assignment, _, _, _ = global_index.resolve_case(rows, set(), [2024, 2025])
        self.assertTrue(assignment["has_additional_version_in_loaded_data"])
        self.assertTrue(assignment["has_later_year_followup_in_loaded_data"])
        self.assertTrue(assignment["has_later_followup_in_loaded_data"])

    def test_conflicts_gaps_chronology_and_outside_window(self):
        rows = [
            self.observation(version="3", initial_date="2024-03-01", row_id=1),
            self.observation(
                primaryid="11",
                version="1",
                source_quarter=2,
                initial_date="2024-01-01",
                row_id=2,
            ),
        ]
        assignment, _, _, _ = global_index.resolve_case(rows, set(), [2024])
        self.assertTrue(assignment["init_fda_dt_conflict"])
        self.assertTrue(assignment["version_gap_detected"])
        self.assertTrue(assignment["version_chronology_conflict"])
        outside = [self.observation(initial_date="2022-01-01")]
        assignment, winner, _, _ = global_index.resolve_case(outside, set(), [2024])
        self.assertEqual(
            assignment["assignment_status"], "assigned_outside_study_window"
        )
        self.assertIsNone(winner)

    def test_numeric_and_lexical_primaryid_ties(self):
        numeric = [
            self.observation(primaryid="9", row_id=1),
            self.observation(primaryid="10", row_id=2),
        ]
        _, winner, _, _ = global_index.resolve_case(numeric, set(), [2024])
        self.assertEqual(winner.primaryid, "10")
        lexical = [
            self.observation(primaryid="alpha", row_id=1),
            self.observation(primaryid="beta", row_id=2),
        ]
        _, winner, _, _ = global_index.resolve_case(lexical, set(), [2024])
        self.assertEqual(winner.primaryid, "beta")

    def test_fixed_schemas_and_direct_parquet_batches(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "assignments.parquet"
            writer = global_index.AtomicParquetWriter(
                path, global_index.ASSIGNMENT_SCHEMA, 1
            )
            base = {field.name: None for field in global_index.ASSIGNMENT_SCHEMA}
            first = base | {
                "caseid": "1",
                "selected_caseversion": 1,
                "quality_warnings": [],
            }
            second = base | {
                "caseid": "2",
                "selected_caseversion": None,
                "quality_warnings": [],
            }
            writer.append(first)
            writer.append(second)
            writer.publish()
            result = pd.read_parquet(path)
            self.assertEqual(len(result), 2)
            self.assertEqual(pq_schema(path), global_index.ASSIGNMENT_SCHEMA)

    def test_sqlite_has_no_json_assignment_or_version_staging(self):
        connection = sqlite3.connect(":memory:")
        self.addCleanup(connection.close)
        global_index.configure_staging_database(connection)
        global_index.create_observation_table(connection)
        tables = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        self.assertEqual(tables, {"observations"})

    def test_secondary_indexes_are_created_after_ingestion_setup(self):
        connection = sqlite3.connect(":memory:")
        self.addCleanup(connection.close)
        global_index.configure_staging_database(connection)
        global_index.create_observation_table(connection)
        before = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='index'"
            )
        }
        self.assertNotIn("ix_observations_case", before)
        global_index.create_secondary_indexes(connection)
        after = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='index'"
            )
        }
        self.assertIn("ix_observations_case", after)

    def test_followup_comparison_uses_actual_source_year_presence(self):
        aggregates = global_index.RunAggregates()
        assignment = assignment_fixture(2023)
        rows = [
            self.observation(source_year=2023),
            self.observation(source_year=2025, row_id=2),
        ]
        global_index.update_aggregates(aggregates, assignment, rows, {2023, 2025})
        self.assertEqual(aggregates.removed_followups[2024], 0)
        self.assertEqual(aggregates.removed_followups[2025], 1)

    def test_removed_followups_are_annual_appearances(self):
        aggregates = global_index.RunAggregates()
        assignment = assignment_fixture(2023)
        rows = [
            self.observation(source_year=2023, row_id=1),
            self.observation(source_year=2024, row_id=2),
            self.observation(source_year=2025, row_id=3),
        ]
        global_index.update_aggregates(aggregates, assignment, rows, {2023, 2024, 2025})
        self.assertEqual(sum(aggregates.removed_followups.values()), 2)
        self.assertEqual(aggregates.counts["cases_seen_in_multiple_years"], 1)

    def test_summary_followups_count_unique_cases_not_rows(self):
        aggregates = global_index.RunAggregates()
        assignment = assignment_fixture(2023)
        assignment["has_additional_version_in_loaded_data"] = True
        assignment["has_later_year_followup_in_loaded_data"] = True
        rows = [
            self.observation(source_year=2023, row_id=1),
            self.observation(source_year=2025, row_id=2),
            self.observation(source_year=2025, row_id=3),
        ]
        global_index.update_aggregates(aggregates, assignment, rows, {2023, 2025})
        self.assertEqual(
            aggregates.counts["cases_with_later_year_followups_in_loaded_data"], 1
        )
        self.assertEqual(aggregates.counts["cases_seen_in_multiple_years"], 1)
        self.assertLessEqual(
            aggregates.counts["cases_with_later_year_followups_in_loaded_data"],
            aggregates.counts["cases_seen_in_multiple_years"],
        )

    def test_streaming_counters_and_invariants_match_outputs(self):
        aggregates = global_index.RunAggregates()
        assignment = assignment_fixture(2024)
        rows = [self.observation()]
        global_index.update_aggregates(aggregates, assignment, rows, {2024})
        aggregates.selected_version_rows = 1
        results = global_index.invariant_results(aggregates, staged_rows=1)
        self.assertTrue(all(results.values()))

    def test_phase_runtime_fields_are_reported_by_sample_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, raw_dir = write_fixture(root)
            out_dir, report_dir = root / "out", root / "reports"
            result = global_index.main(
                [
                    "--years",
                    "2024",
                    "--quarters",
                    "1",
                    "--manifest",
                    str(manifest),
                    "--raw-dir",
                    str(raw_dir),
                    "--out-dir",
                    str(out_dir),
                    "--report-dir",
                    str(report_dir),
                    "--work-dir",
                    str(root / "work"),
                    "--chunk-size",
                    "1",
                ]
            )
            self.assertEqual(result, 0)
            summary = json.loads(
                (report_dir / "cross_year_case_resolution_summary.json").read_text()
            )
            required = {
                "ingestion_seconds",
                "index_build_seconds",
                "resolution_seconds",
                "assignment_export_seconds",
                "version_export_seconds",
                "report_seconds",
                "total_runtime_seconds",
            }
            self.assertTrue(required.issubset(summary["phase_runtimes"]))
            self.assertTrue(all(summary["invariant_results"].values()))
            self.assertIn("followup_metric_scopes", summary)
            self.assertIn(
                "all assignment statuses",
                summary["followup_metric_scopes"][
                    "cases_with_additional_versions_in_loaded_data"
                ],
            )
            self.assertFalse(any((root / "work").glob("global-case-*")))


def assignment_fixture(year):
    return {
        "assigned_report_year": year,
        "assignment_method": "init_fda_dt",
        "assignment_status": "included",
        "in_study_window": True,
        "selected_primaryid": "10",
        "selected_using_valid_caseversion": True,
        "primaryids_observed": ["10"],
        "valid_caseversions_observed": [1],
        "has_later_followup_in_loaded_data": True,
        "has_additional_version_in_loaded_data": True,
        "has_later_year_followup_in_loaded_data": False,
        "init_fda_dt_conflict": False,
        "missing_caseversion_count": 0,
        "invalid_caseversion_count": 0,
        "version_gap_detected": False,
        "version_chronology_conflict": False,
        "critical_identifier_conflict": False,
    }


def pq_schema(path):
    import pyarrow.parquet as pq

    return pq.read_schema(path)


def write_fixture(root):
    raw_dir = root / "raw"
    raw_dir.mkdir()
    archive_path = raw_dir / "q1.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr(
            "ASCII/DEMO24Q1.txt",
            "CASEID$PRIMARYID$CASEVERSION$FDA_DT$INIT_FDA_DT\n"
            "1$10$1$20240102$20240101\n"
            "1$11$2$20240202$20240101\n",
        )
    manifest = root / "manifest.csv"
    row = {
        "year": 2024,
        "quarter": 1,
        "source_system": "FAERS",
        "format": "ASCII",
        "url": "https://example.test/q1.zip",
        "source_filename": "q1.zip",
        "local_filename": "q1.zip",
        "expected_inner_demo": "ASCII/DEMO24Q1.txt",
        "expected_inner_drug": "ASCII/DRUG24Q1.txt",
        "expected_inner_reac": "ASCII/REAC24Q1.txt",
        "notes": "test",
    }
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        writer.writerow(row)
    return manifest, raw_dir


if __name__ == "__main__":
    unittest.main()
