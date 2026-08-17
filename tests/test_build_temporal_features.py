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

    def test_stage4_three_complete_lags(self):
        """Stage 4 worked example: support of 2, 4, 7 in 2018-2020 must map
        to count_t3=2, count_t2=4, count_t1=7 for target year 2021 -- the
        first year with three complete lags once 2013-2022 is available."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            gold = root / "gold"
            self._write_year(gold, 2018, [(2018, "CHEMBL1", "CHEMBL2", "r", 2)])
            self._write_year(gold, 2019, [(2019, "CHEMBL1", "CHEMBL2", "r", 4)])
            self._write_year(gold, 2020, [(2020, "CHEMBL1", "CHEMBL2", "r", 7)])
            self._write_year(gold, 2021, [(2021, "CHEMBL1", "CHEMBL2", "r", 9)])
            output = root / "features"
            temporal.build_features(
                [2018, 2019, 2020, 2021], [2021], gold, output,
                root / "report.json", False, "512MB", 1,
            )
            row = pd.read_parquet(output / "target_year=2021" / "part-00000.parquet").iloc[0]
            self.assertEqual(row["count_t1"], 7)
            self.assertEqual(row["count_t2"], 4)
            self.assertEqual(row["count_t3"], 2)
            self.assertEqual(row["historical_total"], 13)
            self.assertEqual(row["years_observed"], 3)
            self.assertEqual(row["target"], 1)

    def test_target_year_support_does_not_leak_into_lag_features(self):
        """Changing the target year's own support must not change count_t1,
        count_t2, count_t3, or historical_total computed for that target."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)

            def build(target_support):
                gold = root / f"gold_{target_support}"
                self._write_year(gold, 2018, [(2018, "CHEMBL1", "CHEMBL2", "r", 2)])
                self._write_year(gold, 2019, [(2019, "CHEMBL1", "CHEMBL2", "r", 4)])
                self._write_year(gold, 2020, [(2020, "CHEMBL1", "CHEMBL2", "r", 7)])
                self._write_year(gold, 2021, [(2021, "CHEMBL1", "CHEMBL2", "r", target_support)])
                output = root / f"features_{target_support}"
                temporal.build_features(
                    [2018, 2019, 2020, 2021], [2021], gold, output,
                    root / f"report_{target_support}.json", False, "512MB", 1,
                )
                return pd.read_parquet(output / "target_year=2021" / "part-00000.parquet").iloc[0]

            low = build(1)
            high = build(500)
            for column in ("count_t1", "count_t2", "count_t3", "historical_total", "years_observed"):
                self.assertEqual(low[column], high[column], column)
            self.assertEqual(low["target"], 1)
            self.assertEqual(high["target"], 1)


if __name__ == "__main__":
    unittest.main()
