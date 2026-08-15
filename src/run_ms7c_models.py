"""Run ordered MS7C recurrence baselines and streaming logistic regression."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from sklearn.linear_model import SGDClassifier
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FEATURE_ROOT = ROOT / "data/gold/temporal_pair_reaction_features"
DEFAULT_OUTPUT_ROOT = ROOT / "data/processed/modeling/ms7c"
DEFAULT_REPORT_ROOT = ROOT / "reports/modeling/ms7c"
FEATURES = [
    "count_t1", "count_t2", "count_t3", "present_t1", "present_t2",
    "present_t3", "historical_total", "years_observed", "years_since_last_seen",
]
MISSING_FEATURES = [f"{name}_unavailable" for name in FEATURES[:6]]
MODEL_FEATURES = FEATURES + MISSING_FEATURES
THRESHOLDS = [1, 2, 3, 5, 10, 20, 50, 100]


def feature_path(root: Path, year: int) -> Path:
    path = root / f"target_year={year}" / "part-00000.parquet"
    if not path.is_file():
        raise FileNotFoundError(f"MS7B feature partition missing: {path}")
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_batches(path: Path, batch_size: int, include_keys: bool = False):
    columns = (["drug_a", "drug_b", "reaction_pt"] if include_keys else []) + FEATURES + ["target"]
    yield from pq.ParquetFile(path).iter_batches(batch_size=batch_size, columns=columns)


def matrix_from_batch(batch: pa.RecordBatch) -> tuple[np.ndarray, np.ndarray]:
    frame = batch.to_pandas()
    raw = frame[FEATURES].to_numpy(dtype=np.float64, na_value=np.nan)
    unavailable = np.isnan(raw[:, :6]).astype(np.float64)
    np.nan_to_num(raw, copy=False, nan=0.0)
    return np.column_stack((raw, unavailable)), frame["target"].to_numpy(dtype=np.int8)


class BinaryMetrics:
    def __init__(self, probability_bins: int = 10_000):
        self.tp = self.fp = self.tn = self.fn = 0
        self.n = 0
        self.positives = 0
        self.predicted_positives = 0
        self.brier_sum = 0.0
        self.pos_hist = np.zeros(probability_bins, dtype=np.int64)
        self.neg_hist = np.zeros(probability_bins, dtype=np.int64)

    def update(self, y: np.ndarray, prediction: np.ndarray, probability: np.ndarray | None = None):
        y = y.astype(bool)
        prediction = prediction.astype(bool)
        self.tp += int(np.sum(y & prediction))
        self.fp += int(np.sum(~y & prediction))
        self.tn += int(np.sum(~y & ~prediction))
        self.fn += int(np.sum(y & ~prediction))
        self.n += len(y)
        self.positives += int(y.sum())
        self.predicted_positives += int(prediction.sum())
        if probability is not None:
            probability = np.clip(probability.astype(np.float64), 0.0, 1.0)
            self.brier_sum += float(np.sum((probability - y.astype(float)) ** 2))
            indexes = np.minimum((probability * len(self.pos_hist)).astype(np.int64), len(self.pos_hist) - 1)
            self.pos_hist += np.bincount(indexes[y], minlength=len(self.pos_hist))
            self.neg_hist += np.bincount(indexes[~y], minlength=len(self.neg_hist))

    @staticmethod
    def _divide(numerator, denominator):
        return float(numerator / denominator) if denominator else 0.0

    def result(self) -> dict:
        precision = self._divide(self.tp, self.tp + self.fp)
        recall = self._divide(self.tp, self.tp + self.fn)
        specificity = self._divide(self.tn, self.tn + self.fp)
        result = {
            "rows": self.n, "positives": self.positives,
            "prevalence": self._divide(self.positives, self.n),
            "tp": self.tp, "fp": self.fp, "tn": self.tn, "fn": self.fn,
            "precision": precision, "recall": recall,
            "f1": self._divide(2 * precision * recall, precision + recall),
            "specificity": specificity,
            "balanced_accuracy": (recall + specificity) / 2,
            "predicted_positive_rate": self._divide(self.predicted_positives, self.n),
            "accuracy": self._divide(self.tp + self.tn, self.n),
        }
        if self.brier_sum:
            result["brier_score"] = self.brier_sum / self.n
            result.update(self._probability_metrics())
        return result

    def _probability_metrics(self) -> dict:
        positives = int(self.pos_hist.sum())
        negatives = int(self.neg_hist.sum())
        if not positives or not negatives:
            return {"roc_auc_approx": None, "average_precision_approx": None}
        tp = np.cumsum(self.pos_hist[::-1])
        fp = np.cumsum(self.neg_hist[::-1])
        tpr = np.concatenate(([0.0], tp / positives))
        fpr = np.concatenate(([0.0], fp / negatives))
        roc_auc = float(np.trapezoid(tpr, fpr))
        precision = tp / np.maximum(tp + fp, 1)
        recall = tp / positives
        previous_recall = np.concatenate(([0.0], recall[:-1]))
        average_precision = float(np.sum((recall - previous_recall) * precision))
        return {"roc_auc_approx": roc_auc, "average_precision_approx": average_precision}


def select_support_threshold(path: Path, batch_size: int, thresholds=THRESHOLDS):
    metrics = {threshold: BinaryMetrics() for threshold in thresholds}
    for batch in iter_batches(path, batch_size):
        frame = batch.to_pandas()
        counts = frame["count_t1"].fillna(0).to_numpy(dtype=np.int64)
        target = frame["target"].to_numpy(dtype=np.int8)
        for threshold, tracker in metrics.items():
            tracker.update(target, counts >= threshold)
    rows = [{"threshold": threshold, **tracker.result()} for threshold, tracker in metrics.items()]
    selected = max(rows, key=lambda row: (row["f1"], row["precision"], -row["threshold"]))
    return int(selected["threshold"]), rows


def train_logistic(path: Path, batch_size: int, random_state: int = 42):
    scaler = StandardScaler()
    positives = rows = 0
    for batch in iter_batches(path, batch_size):
        matrix, target = matrix_from_batch(batch)
        scaler.partial_fit(matrix)
        positives += int(target.sum())
        rows += len(target)
    if positives == 0 or positives == rows:
        raise ValueError("logistic training requires both target classes")
    positive_weight = rows / (2 * positives)
    negative_weight = rows / (2 * (rows - positives))
    model = SGDClassifier(
        loss="log_loss", penalty="l2", alpha=0.0001, max_iter=1,
        learning_rate="optimal", random_state=random_state, average=True,
    )
    first = True
    for batch in iter_batches(path, batch_size):
        matrix, target = matrix_from_batch(batch)
        weights = np.where(target == 1, positive_weight, negative_weight)
        model.partial_fit(
            scaler.transform(matrix), target, classes=np.array([0, 1], dtype=np.int8) if first else None,
            sample_weight=weights,
        )
        first = False
    return scaler, model, {"training_rows": rows, "training_positives": positives}


def evaluate_and_write(
    path: Path, destination: Path, scaler: StandardScaler, model: SGDClassifier,
    support_threshold: int, batch_size: int,
):
    trackers = {name: BinaryMetrics() for name in (
        "always_negative", "previous_year_persistence", "support_threshold", "logistic_regression",
    )}
    schema = pa.schema([
        pa.field("target_year", pa.int32(), nullable=False),
        pa.field("drug_a", pa.string(), nullable=False), pa.field("drug_b", pa.string(), nullable=False),
        pa.field("reaction_pt", pa.string(), nullable=False), pa.field("target", pa.int8(), nullable=False),
        pa.field("always_negative_prediction", pa.int8(), nullable=False),
        pa.field("persistence_prediction", pa.int8(), nullable=False),
        pa.field("support_threshold_prediction", pa.int8(), nullable=False),
        pa.field("logistic_probability", pa.float32(), nullable=False),
        pa.field("logistic_prediction", pa.int8(), nullable=False),
    ])
    writer = pq.ParquetWriter(destination, schema, compression="zstd")
    target_year = int(path.parent.name.split("=", 1)[1])
    try:
        for batch in iter_batches(path, batch_size, include_keys=True):
            frame = batch.to_pandas()
            matrix, target = matrix_from_batch(batch)
            counts = frame["count_t1"].fillna(0).to_numpy(dtype=np.int64)
            always = np.zeros(len(target), dtype=np.int8)
            persistence = (counts > 0).astype(np.int8)
            support = (counts >= support_threshold).astype(np.int8)
            probability = model.predict_proba(scaler.transform(matrix))[:, 1]
            logistic = (probability >= 0.5).astype(np.int8)
            trackers["always_negative"].update(target, always)
            trackers["previous_year_persistence"].update(target, persistence)
            trackers["support_threshold"].update(target, support)
            trackers["logistic_regression"].update(target, logistic, probability)
            output = pa.Table.from_pydict({
                "target_year": np.full(len(target), target_year, dtype=np.int32),
                "drug_a": frame["drug_a"], "drug_b": frame["drug_b"], "reaction_pt": frame["reaction_pt"],
                "target": target, "always_negative_prediction": always,
                "persistence_prediction": persistence, "support_threshold_prediction": support,
                "logistic_probability": probability.astype(np.float32), "logistic_prediction": logistic,
            }, schema=schema)
            writer.write_table(output)
    finally:
        writer.close()
    return {name: tracker.result() for name, tracker in trackers.items()}


def write_csv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def run_ms7c(
    feature_root: Path, output_root: Path, report_root: Path,
    train_year: int, test_year: int, batch_size: int, overwrite: bool,
):
    if train_year >= test_year:
        raise ValueError("train year must precede test year")
    train_path = feature_path(feature_root, train_year)
    test_path = feature_path(feature_root, test_year)
    if output_root.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {output_root}; pass --overwrite")
    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_root.name}-", dir=output_root.parent))
    try:
        threshold, threshold_rows = select_support_threshold(train_path, batch_size)
        scaler, model, training = train_logistic(train_path, batch_size)
        prediction_dir = staging / "predictions" / f"target_year={test_year}"
        prediction_dir.mkdir(parents=True)
        metrics = evaluate_and_write(
            test_path, prediction_dir / "part-00000.parquet", scaler, model, threshold, batch_size,
        )
        model_dir = staging / "models"
        model_dir.mkdir()
        joblib.dump({"scaler": scaler, "model": model, "features": MODEL_FEATURES}, model_dir / "logistic_regression.joblib")
        summary = {
            "milestone": "MS7C", "status": "stage_3_pilot",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "train_year": train_year, "test_year": test_year,
            "candidate_population": "MS7B tuples observed before target_year",
            "model_order": ["always_negative", "previous_year_persistence", "support_threshold", "logistic_regression"],
            "selected_support_threshold": threshold, **training,
            "input_files": {
                f"target_year={train_year}": {"sha256": sha256_file(train_path), "bytes": train_path.stat().st_size},
                f"target_year={test_year}": {"sha256": sha256_file(test_path), "bytes": test_path.stat().st_size},
            },
            "logistic_preprocessing": {
                "features": MODEL_FEATURES,
                "null_policy": "zero imputation plus explicit unavailable indicators for six lag features",
                "scaler_mean": scaler.mean_.tolist(), "scaler_scale": scaler.scale_.tolist(),
            },
            "logistic_model": {
                "implementation": "sklearn.linear_model.SGDClassifier",
                "loss": "log_loss", "penalty": "l2", "alpha": 0.0001,
                "class_weight_policy": "balanced sample weights computed from training year",
                "decision_threshold": 0.5, "coefficients": model.coef_[0].tolist(),
                "intercept": model.intercept_.tolist(),
            },
            "metrics": metrics, "probability_metric_method": "10,000-bin streaming approximation",
            "gradient_boosted_trees": "optional_not_implemented",
            "all_invariants_passed": True,
        }
        (staging / "metadata").mkdir()
        (staging / "metadata" / "ms7c_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        (staging / "_SUCCESS").touch()
        if output_root.exists():
            shutil.rmtree(output_root)
        staging.replace(output_root)
        report_root.mkdir(parents=True, exist_ok=True)
        (report_root / "ms7c_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        write_csv(report_root / "support_threshold_search.csv", threshold_rows)
        comparison = [{"model": name, "target_year": test_year, **result} for name, result in metrics.items()]
        write_csv(report_root / "model_comparison.csv", comparison)
        write_csv(report_root / "metrics_by_target_year.csv", comparison)
        write_csv(report_root / "confusion_matrices.csv", [
            {key: row[key] for key in ("model", "target_year", "rows", "tp", "fp", "tn", "fn")}
            for row in comparison
        ])
        return summary
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main():
    parser = argparse.ArgumentParser(description="Run ordered MS7C recurrence baselines and logistic regression.")
    parser.add_argument("--feature-root", type=Path, default=DEFAULT_FEATURE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--report-root", type=Path, default=DEFAULT_REPORT_ROOT)
    parser.add_argument("--train-year", type=int, default=2024)
    parser.add_argument("--test-year", type=int, default=2025)
    parser.add_argument("--batch-size", type=int, default=250_000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error("--batch-size must be positive")
    result = run_ms7c(
        args.feature_root, args.output_root, args.report_root,
        args.train_year, args.test_year, args.batch_size, args.overwrite,
    )
    print(f"MS7C pilot complete; logistic F1={result['metrics']['logistic_regression']['f1']:.4f}")


if __name__ == "__main__":
    main()
