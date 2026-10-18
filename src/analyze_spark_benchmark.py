"""MS8C: speedup/efficiency figures and summary tables from run_spark_benchmark.py output.

Outputs (figures via plot_paper_figures.save, PNG + PDF in reports/figures/paper/):
  fig_spark_speedup             speedup and efficiency vs cores, ideal line, Amdahl fit
  fig_spark_runtime_by_year     runtime vs join size per year, one line per core count
  reports/spark/benchmark_summary.json   per-core totals, speedup, efficiency, serial fraction
  reports/spark/benchmark_by_year.csv    local pipeline vs Spark 1 core vs best core count
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_paper_figures import INK, MUTED, SERIES, line, save

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports" / "spark"
MS6_REPORTS = [ROOT / "reports/data_quality/historical_2013_2022", ROOT / "reports/data_quality"]
# Ordinal blue ramp (light -> dark = few -> many cores); lightest step clears 2:1 on white.
CORE_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]


def amdahl_serial_fraction(cores: np.ndarray, speedup: np.ndarray) -> float:
    """Least-squares fit of 1/S = s + (1 - s)/p; returns the serial fraction s."""
    slope, intercept = np.polyfit(1 / cores, 1 / speedup, 1)
    return float(intercept / (intercept + slope))


def summarize(report_dir: Path) -> dict:
    raw = pd.read_csv(report_dir / "benchmark_scaling.csv")
    raw["repeat"] = raw.groupby(["year", "cores"]).cumcount() + 1
    # Repeats are averaged per (year, cores) before anything is summed.
    scaling = raw.groupby(["year", "cores"], as_index=False).agg(
        seconds=("seconds", "mean"), output_rows=("output_rows", "first"))
    totals = scaling.groupby("cores")["seconds"].sum().sort_index()
    repeat_totals = raw.groupby(["cores", "repeat"])["seconds"].agg(["sum", "count"])
    complete = repeat_totals[repeat_totals["count"] == raw["year"].nunique()]["sum"]
    total_sd = complete.groupby("cores").std()
    speedup = totals.iloc[0] / totals
    cores = totals.index.to_numpy(dtype=float)
    summary = {
        "years": f"{scaling['year'].min()}-{scaling['year'].max()}",
        "output_rows": int(scaling.groupby("year")["output_rows"].first().sum()),
        "a_repeats": int(raw["repeat"].max()),
        "by_cores": [{"cores": int(c), "seconds": round(float(t), 1),
                      "sd_seconds": None if pd.isna(total_sd.get(c)) else round(float(total_sd[c]), 1),
                      "speedup": round(float(s), 3), "efficiency": round(float(s / c), 3)}
                     for c, t, s in zip(totals.index, totals, speedup)],
        "serial_fraction_amdahl": round(amdahl_serial_fraction(cores, speedup.to_numpy()), 3)
        if len(cores) > 2 else None,
    }
    repeats_path = report_dir / "benchmark_repeats.csv"
    if repeats_path.exists():
        repeats = pd.read_csv(repeats_path).groupby("cores")["seconds"].agg(["mean", "std", "count"])
        summary["repeats"] = [{"cores": int(c), "mean_seconds": round(r["mean"], 1),
                               "sd_seconds": round(r["std"], 1) if r["count"] > 1 else None,
                               "n": int(r["count"])} for c, r in repeats.iterrows()]
    noop_path = report_dir / "benchmark_noop.csv"
    if noop_path.exists() and repeats_path.exists():
        noop = pd.read_csv(noop_path).groupby("cores")["seconds"].mean()
        summary["write_share"] = {int(c): round(1 - noop[c] / repeats.loc[c, "mean"], 3)
                                  for c in noop.index if c in repeats.index}
    return summary, scaling, totals, speedup


def join_rows() -> pd.Series:
    """Case-pair-reaction rows per year (the join that dominates the work), from the MS6 reports."""
    frames = [pd.read_csv(d / "ms6_annual_pair_reaction_validation.csv") for d in MS6_REPORTS]
    return pd.concat(frames).groupby("year")["case_pair_reaction_rows"].first()


def by_year_table(scaling: pd.DataFrame) -> pd.DataFrame:
    local = pd.concat(pd.read_csv(d / "ms6_stage_performance_by_year.csv") for d in MS6_REPORTS)
    local = local.groupby("year")["duration_seconds"].sum().rename("local_seconds")
    spark = scaling.pivot_table(index="year", columns="cores", values="seconds")
    best = spark.min(axis=1)
    table = pd.DataFrame({
        "output_rows": scaling.groupby("year")["output_rows"].first(),
        "join_rows": join_rows(),
        "local_seconds": local,
        "spark_1core_seconds": spark[spark.columns.min()],
        "spark_best_seconds": best,
        "spark_best_cores": spark.idxmin(axis=1),
    }).dropna(subset=["spark_1core_seconds"])
    table["local_vs_spark_1core"] = table["local_seconds"] / table["spark_1core_seconds"]
    table["spark_1core_vs_best"] = table["spark_1core_seconds"] / table["spark_best_seconds"]
    return table.round(3)


def speedup_figure(totals: pd.Series, speedup: pd.Series, serial_fraction: float | None) -> None:
    cores = totals.index.to_numpy(dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), layout="constrained")
    axes[0].plot(cores, cores, color=MUTED, linestyle="--", linewidth=1, label="Ideal (linear)")
    line(axes[0], cores, speedup, 0, "Measured, 2013-2025")
    if serial_fraction is not None:
        grid = np.linspace(1, cores.max(), 100)
        axes[0].plot(grid, 1 / (serial_fraction + (1 - serial_fraction) / grid), color=SERIES[1][0],
                     linewidth=1.2, label=f"Amdahl fit (serial {serial_fraction:.0%})")
    line(axes[1], cores, speedup / cores, 0, "Measured")
    axes[1].axhline(1, color=MUTED, linestyle="--", linewidth=1)
    axes[0].set_title("Speedup", fontsize=9, color=INK, loc="left")
    axes[1].set_title("Parallel efficiency", fontsize=9, color=INK, loc="left")
    axes[0].set_ylim(0, cores.max() * 1.05)
    axes[1].set_ylim(0, 1.1)
    for ax in axes:
        ax.set_xscale("log", base=2)
        ax.set_xticks(cores, [str(int(c)) for c in cores])
        ax.set_xlabel("Cores (Spark local[N])")
    if 12 < cores.max():
        for ax in axes:
            ax.axvline(12, color=MUTED, linewidth=0.6, linestyle=":")
        axes[1].annotate("12 physical cores", (12, 1.03), xytext=(3, 0), textcoords="offset points",
                         fontsize=7, color=MUTED)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3)
    save(fig, "fig_spark_speedup")


def runtime_by_year_figure(scaling: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(4.8, 3.2), layout="constrained")
    sizes = join_rows()
    core_counts = sorted(scaling["cores"].unique())
    ramp = CORE_RAMP[-len(core_counts):] if len(core_counts) <= len(CORE_RAMP) else CORE_RAMP
    for color, cores in zip(ramp, core_counts):
        rows = scaling[scaling["cores"] == cores].assign(join=lambda d: d["year"].map(sizes)).sort_values("join")
        ax.plot(rows["join"] / 1e6, rows["seconds"], color=color, marker="o", markersize=3.5,
                linewidth=1.4, label=f"{cores} core{'s' if cores > 1 else ''}",
                markeredgecolor="white", markeredgewidth=0.5)
    ax.set_xlabel("Case-pair-reaction join rows per year (millions)")
    ax.set_ylabel("Seconds per year")
    ax.legend(title="Spark local[N]", fontsize=7, title_fontsize=7)
    save(fig, "fig_spark_runtime_by_year")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--report-dir", type=Path, default=REPORTS)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()

    summary, scaling, totals, speedup = summarize(args.report_dir)
    table = by_year_table(scaling)
    table.to_csv(args.report_dir / "benchmark_by_year.csv")
    (args.report_dir / "benchmark_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if not args.no_figures:
        speedup_figure(totals, speedup, summary["serial_fraction_amdahl"])
        runtime_by_year_figure(scaling)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
