import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import run_ms7c_models as modeling  # noqa: E402


class Ms7cModelTests(unittest.TestCase):
    def test_binary_metrics(self):
        metrics = modeling.BinaryMetrics()
        metrics.update(np.array([1, 1, 0, 0]), np.array([1, 0, 1, 0]))
        result = metrics.result()
        self.assertEqual((result["tp"], result["fp"], result["tn"], result["fn"]), (1, 1, 1, 1))
        self.assertEqual(result["f1"], 0.5)

    def test_end_to_end_small_pilot(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            features = root / "features"
            columns = ["target_year", "drug_a", "drug_b", "reaction_pt", *modeling.FEATURES, "target"]
            for year, rows in {
                2024: [
                    (2024, "CHEMBL1", "CHEMBL2", "a", 5, None, None, 1, None, None, 5, 1, 1, 1),
                    (2024, "CHEMBL1", "CHEMBL3", "b", 1, None, None, 1, None, None, 1, 1, 1, 0),
                    (2024, "CHEMBL2", "CHEMBL3", "c", 10, None, None, 1, None, None, 10, 1, 1, 1),
                    (2024, "CHEMBL2", "CHEMBL4", "d", 2, None, None, 1, None, None, 2, 1, 1, 0),
                ],
                2025: [
                    (2025, "CHEMBL1", "CHEMBL2", "a", 4, 5, None, 1, 1, None, 9, 2, 1, 1),
                    (2025, "CHEMBL1", "CHEMBL3", "b", 0, 1, None, 0, 1, None, 1, 1, 2, 0),
                ],
            }.items():
                partition = features / f"target_year={year}"
                partition.mkdir(parents=True)
                pd.DataFrame(rows, columns=columns).to_parquet(partition / "part-00000.parquet", index=False)
            result = modeling.run_ms7c(features, root / "output", root / "reports", 2024, 2025, 2, False)
            self.assertEqual(result["model_order"][0], "always_negative")
            self.assertTrue((root / "output" / "_SUCCESS").exists())
            prediction = pd.read_parquet(root / "output" / "predictions" / "target_year=2025" / "part-00000.parquet")
            self.assertEqual(len(prediction), 2)
            self.assertTrue(prediction["logistic_probability"].between(0, 1).all())


if __name__ == "__main__":
    unittest.main()
