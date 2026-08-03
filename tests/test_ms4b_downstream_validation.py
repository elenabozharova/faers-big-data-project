"""Focused MS4B boundary-contract tests for the downstream rebuild."""

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import rebuild_global_downstream as rebuild  # noqa: E402
from tests.test_full_global_downstream_rebuild import fixture  # noqa: E402


class Ms4bDownstreamValidationTests(unittest.TestCase):
    def test_nonselected_primaryid_cannot_enter_reused_boundary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            boundary, assignments, _mapping = fixture(root)
            drug_path = boundary / "2099/resolved_drug_rows.parquet"
            drug = rebuild.pd.read_parquet(drug_path)
            drug.loc[0, "primaryid"] = "later-nonselected"
            drug.to_parquet(drug_path, index=False)
            with self.assertRaisesRegex(AssertionError, "nonselected PRIMARYID"):
                rebuild.load_boundary(boundary, 2099, assignments)

    def test_output_root_must_not_be_boundary_root(self):
        self.assertNotEqual(rebuild.DEFAULT_BOUNDARY_ROOT, rebuild.DEFAULT_OUTPUT_ROOT)


if __name__ == "__main__":
    unittest.main()
