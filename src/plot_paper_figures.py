"""Paper figures from the Historical Baseline v1 reports (PNG 300 dpi + PDF)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports" / "modeling" / "historical_baseline_v1"
OUT = ROOT / "reports" / "figures" / "paper"

# Validated categorical slots 1-3 (all-pairs CVD safe); markers carry identity in grayscale print.
SERIES = [("#2a78d6", "o"), ("#eb6834", "s"), ("#1baf7a", "^")]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def line(ax, years, values, slot, label):
    color, marker = SERIES[slot]
    ax.plot(years, values, color=color, marker=marker, markersize=4.5, linewidth=1.6, label=label,
            markeredgecolor="white", markeredgewidth=0.6)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)


def ablation_figure():
    wide = pd.read_csv(REPORTS / "ablation_per_year_wide.csv")
    years = wide["target_year"]
    panels = [("f1", "F1 (threshold 0.5)"), ("roc_auc_approx", "AUROC"), ("average_precision_approx", "AUPRC")]
    sets = [("m1", "M1: count_t1"), ("m2", "M2: 3 annual lags"), ("m3", "M3: 9 temporal features")]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), sharex=True, layout="constrained")
    for ax, (metric, title) in zip(axes, panels):
        for slot, (fs, label) in enumerate(sets):
            line(ax, years, wide[f"{fs}_{metric}"], slot, label)
        ax.set_title(title, fontsize=9, color=INK, loc="left")
        ax.set_xticks(years[::3])
    axes[0].set_ylabel("Score")
    axes[1].set_xlabel("Target year")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3)
    save(fig, "fig_ablation_m1_m2_m3")


def precision_recall_figure():
    detail = pd.read_csv(REPORTS / "all_metrics_detailed.csv")
    m3 = detail[detail["feature_set"] == "m3"]
    models = [("logistic_regression", "Logistic regression (M3)"),
              ("support_threshold", "Support-threshold rule"),
              ("previous_year_persistence", "Previous-year persistence")]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharex=True, sharey=True, layout="constrained")
    for slot, (model, label) in enumerate(models):
        rows = m3[m3["model"] == model].sort_values("target_year")
        line(axes[0], rows["target_year"], rows["precision"], slot, label)
        line(axes[1], rows["target_year"], rows["recall"], slot, label)
    base = m3[m3["model"] == "logistic_regression"].sort_values("target_year")
    axes[0].plot(base["target_year"], base["prevalence"], color=MUTED, linestyle="--", linewidth=1)
    axes[0].annotate("chance level (prevalence)", (base["target_year"].iloc[-1], base["prevalence"].iloc[-1]),
                     xytext=(-4, -12), textcoords="offset points", ha="right", fontsize=8, color=MUTED)
    axes[0].set_title("Precision", fontsize=9, color=INK, loc="left")
    axes[1].set_title("Recall", fontsize=9, color=INK, loc="left")
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("Score")
    for ax in axes:
        ax.set_xlabel("Target year")
        ax.set_xticks(base["target_year"][::3])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3)
    save(fig, "fig_precision_recall_by_year")


if __name__ == "__main__":
    ablation_figure()
    precision_recall_figure()
    print(f"Figures written to {OUT}")
