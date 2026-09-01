import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import run_historical_baseline as driver  # noqa: E402
import run_ms7c_models as ms7c  # noqa: E402


def synthetic_partition(year: int) -> list[tuple]:
    columns = ["target_year", "drug_a", "drug_b", "reaction_pt", *ms7c.FEATURES, "target"]
    rows = [
        (year, "CHEMBL1", "CHEMBL2", "a", 5, 3, 1, 1, 1, 1, 9, 3, 1, 1),
        (year, "CHEMBL1", "CHEMBL3", "b", 1, 0, None, 1, 0, None, 1, 2, 1, 0),
        (year, "CHEMBL2", "CHEMBL3", "c", 10, 8, 6, 1, 1, 1, 24, 3, 1, 1),
        (year, "CHEMBL2", "CHEMBL4", "d", 0, 2, 0, 0, 1, 0, 2, 2, 2, 0),
    ]
    return columns, rows


class RunHistoricalBaselineTests(unittest.TestCase):
    def test_rolling_years_and_ablation_wiring(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            features = root / "features"
            for year in (2015, 2016, 2017):
                columns, rows = synthetic_partition(year)
                partition = features / f"target_year={year}"
                partition.mkdir(parents=True)
                pd.DataFrame(rows, columns=columns).to_parquet(partition / "part-00000.parquet", index=False)

            manifest = driver.run_all(
                features, root / "output", root / "reports",
                target_years=[2016, 2017], batch_size=2, overwrite=False,
            )

            self.assertEqual(
                [row["train_year"] for row in manifest["metrics_by_target_year"]], [2015, 2016]
            )
            self.assertEqual(
                [row["target_year"] for row in manifest["metrics_by_target_year"]], [2016, 2017]
            )
            # M1/M2/M3 must produce distinct coefficient-vector lengths (1, 3, 9 raw features).
            for feature_set, expected_len in (("m1", 1), ("m2", 3), ("m3", 9)):
                predictions = root / "output" / feature_set / "target_year=2016" / "predictions" / "target_year=2016" / "part-00000.parquet"
                self.assertTrue(predictions.exists())
                _, model_features = ms7c.model_features_for(feature_set)
                self.assertEqual(len(ms7c.FEATURE_SETS[feature_set]), expected_len)
                self.assertTrue(len(model_features) >= expected_len)

            self.assertIn("m1_mean_f1", manifest["ablation"])
            self.assertIn("m2_mean_f1", manifest["ablation"])
            self.assertIn("m3_mean_f1", manifest["ablation"])
            self.assertEqual(manifest["recency_consistency"]["years_checked"], 2)

            report_root = root / "reports"
            self.assertTrue((report_root / "metrics_by_target_year.csv").exists())
            self.assertTrue((report_root / "coefficients_by_target_year.csv").exists())
            self.assertTrue((report_root / "ablation_comparison.csv").exists())
            self.assertTrue((report_root / "historical_baseline_v1_manifest.json").exists())


if __name__ == "__main__":
    unittest.main()
