import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import build_temporal_features as temporal  # noqa: E402


class TemporalFeatureTests(unittest.TestCase):
    def _write_year(self, root, year, rows):
        partition = root / f"year={year}"
        partition.mkdir(parents=True)
        pd.DataFrame(rows, columns=["year", "drug_a", "drug_b", "reaction_pt", "n_cases"]).to_parquet(
            partition / "part-00000.parquet", index=False,
        )

    def test_builds_leakage_safe_pilot_features(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            gold = root / "gold"
            self._write_year(gold, 2023, [
                (2023, "CHEMBL1", "CHEMBL2", "r1", 2),
                (2023, "CHEMBL1", "CHEMBL3", "r2", 4),
            ])
            self._write_year(gold, 2024, [
                (2024, "CHEMBL1", "CHEMBL2", "r1", 3),
                (2024, "CHEMBL2", "CHEMBL3", "new", 5),
            ])
            self._write_year(gold, 2025, [
                (2025, "CHEMBL1", "CHEMBL3", "r2", 1),
            ])
            output = root / "features"
            result = temporal.build_features(
                [2023, 2024, 2025], [2024, 2025], gold, output,
                root / "report.json", False, "512MB", 1,
            )
            y2024 = pd.read_parquet(output / "target_year=2024" / "part-00000.parquet")
            self.assertEqual(len(y2024), 2)
            self.assertEqual(y2024["target"].sum(), 1)
            self.assertTrue(y2024["count_t2"].isna().all())
            y2025 = pd.read_parquet(output / "target_year=2025" / "part-00000.parquet")
            repeated = y2025.loc[y2025["reaction_pt"].eq("r1")].iloc[0]
            self.assertEqual(repeated["count_t1"], 3)
            self.assertEqual(repeated["count_t2"], 2)
            self.assertEqual(repeated["historical_total"], 5)
            self.assertEqual(repeated["years_observed"], 2)
            self.assertEqual(result["total_rows"], 5)

    def test_rejects_target_without_history(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self._write_year(root, 2024, [(2024, "CHEMBL1", "CHEMBL2", "r", 1)])
            with self.assertRaisesRegex(ValueError, "no history"):
                temporal.build_features([2024], [2024], root, root / "out", root / "r.json", False)


if __name__ == "__main__":
    unittest.main()
