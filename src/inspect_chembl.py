"""
Inspect ChEMBL processed artifacts used by build_decagon_graph.py.

Reads from data/processed/chembl  (mechanism-only / strict)
      and data/processed/chembl_all (union with activities / broad)

Run after process_chembl.py has written its outputs.
"""

import json
import os
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")

MODES = {
    "strict (drug_mechanism)": "chembl",
    "broad  (mechanism + activities)": "chembl_all",
}


# ---------- Loaders ----------

def load_pairs(artifact_dir: str) -> pd.DataFrame:
    path = os.path.join(artifact_dir, "drug_target_pairs.tsv")
    return pd.read_csv(path, sep="\t", dtype=str)


def load_json(artifact_dir: str, filename: str) -> dict:
    with open(os.path.join(artifact_dir, filename), encoding="utf-8") as f:
        return json.load(f)


# ---------- Display helpers ----------

def sep(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def print_stats(stats: dict) -> None:
    print(f"  Unique drugs:        {stats['n_unique_drugs']:>8,}")
    print(f"  Unique targets:      {stats['n_unique_targets']:>8,}")
    print(f"  Unique pairs:        {stats['n_unique_pairs']:>8,}")
    print(f"  Avg targets/drug:    {stats['avg_targets_per_drug']:>8.2f}")
    print(f"  Avg drugs/target:    {stats['avg_drugs_per_target']:>8.2f}")

    dist = stats.get("action_type_distribution", {})
    if dist:
        total = sum(dist.values())
        print(f"\n  Action / source type breakdown (n={total:,}):")
        for action, cnt in sorted(dist.items(), key=lambda x: -x[1]):
            print(f"    {action:30s}: {cnt:6,}  ({100*cnt/total:.1f}%)")


def top_n_by_degree(
    pairs: pd.DataFrame,
    group_col: str,
    count_col: str,
    label: str,
    n: int = 20,
) -> pd.DataFrame:
    degree = (
        pairs.groupby(group_col)[count_col]
        .nunique()
        .rename("degree")
        .sort_values(ascending=False)
        .head(n)
        .reset_index()
    )
    degree.columns = [group_col, f"n_{label}"]
    return degree


def save_csv(df: pd.DataFrame, name: str) -> None:
    os.makedirs(os.path.join(DATA_PROCESSED), exist_ok=True)
    out = os.path.join(DATA_PROCESSED, f"{name}.csv")
    df.to_csv(out, index=False)
    print(f"  Saved: {out}")


# ---------- Main inspection ----------

def inspect_mode(label: str, subdir: str) -> None:
    artifact_dir = os.path.join(DATA_PROCESSED, subdir)

    if not os.path.isdir(artifact_dir):
        print(f"  [SKIP] Directory not found: {artifact_dir}")
        return

    pairs = load_pairs(artifact_dir)
    stats = load_json(artifact_dir, "stats.json")
    meta  = load_json(artifact_dir, "meta.json")
    drug_to_idx   = load_json(artifact_dir, "drug_to_idx.json")
    target_to_idx = load_json(artifact_dir, "target_to_idx.json")

    # ── 1. Coverage summary ───────────────────────────────────────────
    sep(f"ChEMBL [{label}]  —  coverage summary")
    print_stats(stats)
    print(f"\n  Source: {meta.get('source', 'n/a')}")
    print(f"  Index sizes  →  drugs: {len(drug_to_idx):,}   targets: {len(target_to_idx):,}")

    # ── 2. Sample pairs ───────────────────────────────────────────────
    sep(f"ChEMBL [{label}]  —  drug_target_pairs.tsv  (first 100 rows)")
    print(pairs.head(100).to_string(index=False))
    save_csv(pairs.head(100), f"chembl_{subdir}_pairs_top100")

    # ── 3. Top drugs by target count ─────────────────────────────────
    sep(f"ChEMBL [{label}]  —  top 20 drugs by number of targets")
    top_drugs = top_n_by_degree(pairs, "molecule_chembl_id", "uniprot_accession", "targets")
    print(top_drugs.to_string(index=False))
    save_csv(top_drugs, f"chembl_{subdir}_top_drugs")

    # ── 4. Top targets by drug count ─────────────────────────────────
    sep(f"ChEMBL [{label}]  —  top 20 targets by number of drugs")
    top_targets = top_n_by_degree(pairs, "uniprot_accession", "molecule_chembl_id", "drugs")
    print(top_targets.to_string(index=False))
    save_csv(top_targets, f"chembl_{subdir}_top_targets")

    # ── 5. Degree distribution summary ───────────────────────────────
    sep(f"ChEMBL [{label}]  —  degree distribution (percentiles)")
    drug_deg  = pairs.groupby("molecule_chembl_id")["uniprot_accession"].nunique()
    tgt_deg   = pairs.groupby("uniprot_accession")["molecule_chembl_id"].nunique()

    for name, series in [("targets per drug", drug_deg), ("drugs per target", tgt_deg)]:
        desc = series.describe(percentiles=[0.25, 0.5, 0.75, 0.90, 0.99])
        print(f"\n  {name}:")
        for stat, val in desc.items():
            print(f"    {stat:6s}: {val:.1f}")


if __name__ == "__main__":
    os.makedirs(DATA_PROCESSED, exist_ok=True)

    for label, subdir in MODES.items():
        inspect_mode(label, subdir)

    print("\n" + "=" * 80)
    print("Done.")
    print("=" * 80)
