"""Roll the existing MS7C evaluation across every target year with three real
lags (2016-2025) and every ablation feature set (M1/M2/M3), then freeze the
result as Historical Baseline v1.

This is a thin driver over run_ms7c_models.run_ms7c -- it does not
reimplement modeling, preprocessing, or metrics.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import run_ms7c_models as ms7c
from integrate_global_case_assignments import repository_path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FEATURE_ROOT = ROOT / "data/gold/temporal_pair_reaction_features_historical_v1"
DEFAULT_OUTPUT_ROOT = ROOT / "data/processed/modeling/historical_baseline_v1"
DEFAULT_REPORT_ROOT = ROOT / "reports/modeling/historical_baseline_v1"
TARGET_YEARS = list(range(2016, 2026))
PRIMARY_FEATURE_SET = "m3"
COEFFICIENT_FEATURE_SETS = ("m2", "m3")
AGGREGATE_METRIC_KEYS = [
    ("support_threshold", "precision"), ("support_threshold", "recall"), ("support_threshold", "f1"),
    ("logistic_regression", "precision"), ("logistic_regression", "recall"), ("logistic_regression", "f1"),
    ("logistic_regression", "roc_auc_approx"), ("logistic_regression", "average_precision_approx"),
]
DETAILED_METRIC_FIELDS = [
    "rows", "positives", "prevalence", "tp", "fp", "tn", "fn", "precision", "recall", "f1",
    "specificity", "balanced_accuracy", "predicted_positive_rate", "accuracy",
    "brier_score", "roc_auc_approx", "average_precision_approx",
]
ABLATION_METRICS = ("f1", "roc_auc_approx", "average_precision_approx")


def _git_metadata() -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout.strip()
    dirty = bool(subprocess.run(
        ["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout)
    return commit, dirty


def write_csv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def run_all(
    feature_root: Path, output_root: Path, report_root: Path,
    target_years: list[int], batch_size: int, overwrite: bool,
) -> dict:
    if len(target_years) != len(set(target_years)):
        raise ValueError("target years must be unique")
    runs = {}
    for feature_set in sorted(ms7c.FEATURE_SETS):
        for test_year in target_years:
            train_year = test_year - 1
            run_output = output_root / feature_set / f"target_year={test_year}"
            run_report = report_root / feature_set / f"target_year={test_year}"
            summary = ms7c.run_ms7c(
                feature_root, run_output, run_report,
                train_year, test_year, batch_size, overwrite, feature_set,
            )
            runs[(feature_set, test_year)] = summary

    metrics_by_year = []
    coefficients_by_year = []
    detailed_metrics = []
    for feature_set in sorted(ms7c.FEATURE_SETS):
        for test_year in target_years:
            summary = runs[(feature_set, test_year)]
            for model_name, metrics in summary["metrics"].items():
                row = {"feature_set": feature_set, "model": model_name, "target_year": test_year}
                for field in DETAILED_METRIC_FIELDS:
                    row[field] = metrics.get(field)
                detailed_metrics.append(row)

    for test_year in target_years:
        summary = runs[(PRIMARY_FEATURE_SET, test_year)]
        metrics = summary["metrics"]
        row = {
            "target_year": test_year, "train_year": summary["train_year"],
            "n_samples": metrics["logistic_regression"]["rows"],
            "positive_count": metrics["logistic_regression"]["positives"],
            "positive_rate": metrics["logistic_regression"]["prevalence"],
            "selected_support_threshold": summary["selected_support_threshold"],
        }
        for model, key in AGGREGATE_METRIC_KEYS:
            row[f"{model}_{key}"] = metrics[model][key]
        metrics_by_year.append(row)

    for feature_set in COEFFICIENT_FEATURE_SETS:
        for test_year in target_years:
            summary = runs[(feature_set, test_year)]
            model_features = summary["logistic_preprocessing"]["features"]
            coefficients = summary["logistic_model"]["coefficients"]
            coefficients_by_year.append({
                "feature_set": feature_set, "target_year": test_year,
                **{f"coef_{name}": value for name, value in zip(model_features, coefficients)},
                "intercept": summary["logistic_model"]["intercept"][0],
            })

    # Coefficient dominance: rank each model_feature by |coefficient| within its own
    # (feature_set, target_year) model. Rank 1 = largest |coefficient| that year.
    dominance_rows = []
    for feature_set in COEFFICIENT_FEATURE_SETS:
        for test_year in target_years:
            summary = runs[(feature_set, test_year)]
            model_features = summary["logistic_preprocessing"]["features"]
            coefficients = summary["logistic_model"]["coefficients"]
            ranked = sorted(
                zip(model_features, coefficients), key=lambda item: abs(item[1]), reverse=True,
            )
            for rank, (name, value) in enumerate(ranked, start=1):
                dominance_rows.append({
                    "feature_set": feature_set, "target_year": test_year, "feature": name,
                    "coefficient": value, "abs_coefficient": abs(value), "rank_within_model": rank,
                })

    dominance_summary_rows = []
    for feature_set in COEFFICIENT_FEATURE_SETS:
        by_feature: dict[str, list[dict]] = {}
        for row in dominance_rows:
            if row["feature_set"] == feature_set:
                by_feature.setdefault(row["feature"], []).append(row)
        n_years = len(target_years)
        for feature, rows in by_feature.items():
            coefs = np.array([r["coefficient"] for r in rows], dtype=np.float64)
            abs_coefs = np.abs(coefs)
            rank1_count = sum(r["rank_within_model"] == 1 for r in rows)
            dominance_summary_rows.append({
                "feature_set": feature_set, "feature": feature,
                "mean_coefficient": float(coefs.mean()),
                "mean_abs_coefficient": float(abs_coefs.mean()),
                "std_coefficient": float(coefs.std(ddof=1)) if len(coefs) > 1 else 0.0,
                "years_rank1_count": rank1_count, "years_rank1_pct": rank1_count / n_years,
            })

    # Recency check: does |coef(count_t1)| > |coef(count_t2)| > |coef(count_t3)| hold in the M3 model?
    m3_coef_rows = {row["target_year"]: row for row in coefficients_by_year if row["feature_set"] == "m3"}
    recency_consistent_years = sum(
        abs(row["coef_count_t1"]) > abs(row["coef_count_t2"]) > abs(row["coef_count_t3"])
        for row in m3_coef_rows.values()
    )

    ablation_rows = []
    for feature_set in sorted(ms7c.FEATURE_SETS):
        for test_year in target_years:
            metrics = runs[(feature_set, test_year)]["metrics"]["logistic_regression"]
            ablation_rows.append({
                "feature_set": feature_set, "target_year": test_year,
                "features": ",".join(ms7c.FEATURE_SETS[feature_set]),
                "precision": metrics["precision"], "recall": metrics["recall"], "f1": metrics["f1"],
                "roc_auc_approx": metrics.get("roc_auc_approx"),
                "average_precision_approx": metrics.get("average_precision_approx"),
            })

    # Wide per-year table: target_year, then {metric}_{feature_set} for each ablation metric.
    ablation_by_year = {
        feature_set: {row["target_year"]: row for row in ablation_rows if row["feature_set"] == feature_set}
        for feature_set in sorted(ms7c.FEATURE_SETS)
    }
    ablation_per_year_wide = []
    for test_year in target_years:
        wide_row = {"target_year": test_year}
        for feature_set in sorted(ms7c.FEATURE_SETS):
            row = ablation_by_year[feature_set][test_year]
            for metric in ABLATION_METRICS:
                wide_row[f"{feature_set}_{metric}"] = row[metric]
        ablation_per_year_wide.append(wide_row)

    ablation_aggregate_rows = []
    for feature_set in sorted(ms7c.FEATURE_SETS):
        fs_rows = [row for row in ablation_rows if row["feature_set"] == feature_set]
        agg_row = {"feature_set": feature_set}
        for metric in ABLATION_METRICS:
            values = np.array([row[metric] for row in fs_rows], dtype=np.float64)
            agg_row[f"mean_{metric}"] = float(values.mean())
            agg_row[f"std_{metric}"] = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        ablation_aggregate_rows.append(agg_row)

    aggregate = {}
    for model, key in AGGREGATE_METRIC_KEYS:
        values = np.array([row[f"{model}_{key}"] for row in metrics_by_year], dtype=np.float64)
        aggregate[f"{model}_{key}"] = {
            "mean": float(values.mean()), "std": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
            "median": float(np.median(values)), "min": float(values.min()), "max": float(values.max()),
        }

    m1_f1 = np.array([row["f1"] for row in ablation_rows if row["feature_set"] == "m1"])
    m2_f1 = np.array([row["f1"] for row in ablation_rows if row["feature_set"] == "m2"])
    m3_f1 = np.array([row["f1"] for row in ablation_rows if row["feature_set"] == "m3"])

    commit, dirty = _git_metadata()
    manifest = {
        "name": "Historical Baseline v1", "milestone": "post-MS7C historical rolling evaluation",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_commit": commit, "source_tree_dirty": dirty,
        "feature_root": repository_path(feature_root),
        "target_years": target_years,
        "feature_sets": {name: features for name, features in ms7c.FEATURE_SETS.items()},
        "model": {
            "rule_baselines": ["always_negative", "previous_year_persistence", "support_threshold"],
            "logistic": "sklearn.linear_model.SGDClassifier (streaming, same methodology as the Stage 3 pilot)",
            "train_test_split": "train_year = target_year - 1, per target year (no random split)",
        },
        "metrics_by_target_year": metrics_by_year,
        "aggregate_metrics": aggregate,
        "coefficients_by_target_year": coefficients_by_year,
        "feature_scaling": (
            "sklearn.StandardScaler (zero mean, unit variance) is fit per training year "
            "before logistic regression (run_ms7c_models.train_logistic); coefficients are "
            "therefore in standardized-feature space, so magnitudes are directly comparable "
            "across features within a model, and mean/SD comparisons across years are meaningful. "
            "They are NOT comparable to raw, unscaled feature units."
        ),
        "coefficient_dominance": dominance_summary_rows,
        "recency_consistency": {
            "years_checked": len(m3_coef_rows),
            "years_with_count_t1_gt_t2_gt_t3": recency_consistent_years,
        },
        "ablation": {
            "m1_mean_f1": float(m1_f1.mean()), "m2_mean_f1": float(m2_f1.mean()), "m3_mean_f1": float(m3_f1.mean()),
            "m1_vs_m3_mean_f1_delta": float(m3_f1.mean() - m1_f1.mean()),
            "m2_vs_m3_mean_f1_delta": float(m3_f1.mean() - m2_f1.mean()),
            "aggregate_by_feature_set": ablation_aggregate_rows,
        },
        "stage_3_pilot_comparison": "reports/modeling/ms7c/ (untouched, kept for direct comparison)",
    }

    report_root.mkdir(parents=True, exist_ok=True)
    write_csv(report_root / "metrics_by_target_year.csv", metrics_by_year)
    write_csv(report_root / "all_metrics_detailed.csv", detailed_metrics)
    write_csv(report_root / "coefficients_by_target_year.csv", coefficients_by_year)
    write_csv(report_root / "coefficient_dominance.csv", dominance_rows)
    write_csv(report_root / "coefficient_dominance_summary.csv", dominance_summary_rows)
    write_csv(report_root / "ablation_comparison.csv", ablation_rows)
    write_csv(report_root / "ablation_per_year_wide.csv", ablation_per_year_wide)
    write_csv(report_root / "ablation_aggregate.csv", ablation_aggregate_rows)
    (report_root / "historical_baseline_v1_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main():
    parser = argparse.ArgumentParser(description="Roll MS7C across every 3-lag target year and ablation feature set.")
    parser.add_argument("--feature-root", type=Path, default=DEFAULT_FEATURE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--report-root", type=Path, default=DEFAULT_REPORT_ROOT)
    parser.add_argument("--target-years", type=int, nargs="+", default=TARGET_YEARS)
    parser.add_argument("--batch-size", type=int, default=250_000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    manifest = run_all(
        args.feature_root, args.output_root, args.report_root,
        sorted(args.target_years), args.batch_size, args.overwrite,
    )
    print(f"Historical Baseline v1 complete; {len(manifest['metrics_by_target_year'])} target years evaluated.")


if __name__ == "__main__":
    main()
