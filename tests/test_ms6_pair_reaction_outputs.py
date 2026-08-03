"""Focused MS6 usability and report-contract tests."""

import sys
import os
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import rebuild_global_downstream as rebuild  # noqa: E402


class Ms6OutputTests(unittest.TestCase):
    def test_memory_observation_is_available_on_windows(self):
        observed = rebuild._observed_memory_mb()
        if os.name == "nt":
            self.assertIsNotNone(observed)
            self.assertGreater(observed, 0)

    def test_case_usability_deduplicates_ids_and_separates_combination_ambiguity(self):
        retained = pd.DataFrame({"caseid": ["1", "2", "3", "4"]})
        long = pd.DataFrame([
            {"caseid": "1", "normalized_drug_name": "a", "mapping_status": "mapped"},
            {"caseid": "1", "normalized_drug_name": "a", "mapping_status": "mapped"},
            {"caseid": "2", "normalized_drug_name": "a|b", "mapping_status": "mapped_combination_components"},
            {"caseid": "3", "normalized_drug_name": "amb", "mapping_status": "unmapped_ambiguous"},
        ])
        drugs = pd.DataFrame([
            {"caseid": "1", "mapped_drug_id": "CHEMBL1"},
        ])
        result = rebuild.build_case_usability(2023, retained, long, drugs)
        self.assertEqual(result["cases_with_one_single_id_mapped_drug"], 1)
        self.assertEqual(result["cases_with_zero_single_id_mapped_drugs"], 3)
        self.assertEqual(result["cases_with_combination_only_mappings"], 1)
        self.assertEqual(result["cases_with_ambiguous_only_mappings"], 1)


if __name__ == "__main__":
    unittest.main()
