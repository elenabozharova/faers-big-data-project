"""
build_decagon_inputs.py — Assemble Decagon-format input tensors.

Converts the flat count tables from process_faers.py into the sparse data
structures the Decagon GCN expects, alongside the already-built PPI and
drug-target matrices.

Decagon expects (per the original Zitnik et al. 2018 implementation):
  adj_mats[(0,0)] = [ppi_adj]                     shape: (n_prot × n_prot)
  adj_mats[(0,1)] = [drug_target_adj.T]            shape: (n_prot × n_drug)
  adj_mats[(1,0)] = [drug_target_adj]              shape: (n_drug × n_prot)
  adj_mats[(1,1)] = [se_0_adj, se_1_adj, ...]      shape: (n_drug × n_drug) each

This script produces all four blocks and saves them in a single directory.

Combo storage format
--------------------
Per-SE sparse matrices for (n_drugs x n_drugs) are only practical for the
strict drug set (3,471 drugs).  For the broad set (672 k drugs) they are
too large to materialise all at once.

  strict mode  →  per-SE .npz files under combo_adjs/  (feasible)
  all mode     →  COO flat format only; per-SE materialisation is deferred
                  to the training script via load_combo_adjacencies()

Pre-requisites
--------------
  process_faers.py must be run first (produces combo_counts.parquet etc.)

Usage
-----
  python src/build_decagon_inputs.py
  python src/build_decagon_inputs.py --mode strict
  python src/build_decagon_inputs.py --min-se-pairs 50
  python src/build_decagon_inputs.py --top-se 964   # Decagon paper used 964 SEs

Outputs  (data/processed/decagon_inputs[_all]/)
-----------------------------------------------
  combo_adj_coo.npz    COO arrays: drug_i, drug_j, se_idx  (after SE filtering)
  mono_adj.npz         sparse (n_drugs × n_se_mono)
  ppi_adj.npz          symlink/copy from decagon[_all]/
  drug_target_adj.npz  symlink/copy from decagon[_all]/
  se_names.json        filtered se_idx (0-based) -> MedDRA PT string
  mono_se_names.json   filtered mono se_idx -> MedDRA PT string
  meta.json            shape/count statistics
  combo_adjs/          (strict mode only) per-SE sparse .npz files
"""

import argparse
import json
import os
import shutil
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(_SRC_DIR)
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")

FAERS_DIR = os.path.join(DATA_PROCESSED, "faers")


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_faers_outputs(faers_dir: str) -> tuple:
    """Load combo_counts, mono_counts and side_effect_to_idx from process_faers."""
    combo_path = os.path.join(faers_dir, "combo_counts.parquet")
    mono_path = os.path.join(faers_dir, "mono_counts.parquet")
    se_path = os.path.join(faers_dir, "side_effect_to_idx.json")

    for p in (combo_path, mono_path, se_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"Missing FAERS output: {p}\n"
                "Run process_faers.py first."
            )

    combo_df = pd.read_parquet(combo_path)
    mono_df = pd.read_parquet(mono_path)
    with open(se_path, encoding="utf-8") as f:
        se_to_idx = json.load(f)          # pt -> original idx (from process_faers)
    idx_to_se = {v: k for k, v in se_to_idx.items()}

    print(f"  Combo triples:      {len(combo_df):,}")
    print(f"  Mono pairs:         {len(mono_df):,}")
    print(f"  Raw SE vocabulary:  {len(se_to_idx):,}")
    return combo_df, mono_df, idx_to_se


def load_decagon_artifacts(mode: str) -> tuple:
    """Load ppi_adj.npz, drug_target_adj.npz, drug_to_idx, uniprot_to_idx."""
    subdir = "decagon" if mode == "strict" else "decagon_all"
    d = os.path.join(DATA_PROCESSED, subdir)

    ppi = sp.load_npz(os.path.join(d, "ppi_adj.npz"))
    dt = sp.load_npz(os.path.join(d, "drug_target_adj.npz"))

    with open(os.path.join(d, "drug_to_idx.json"), encoding="utf-8") as f:
        drug_to_idx = json.load(f)
    with open(os.path.join(d, "uniprot_to_idx.json"), encoding="utf-8") as f:
        prot_to_idx = json.load(f)

    n_drugs = len(drug_to_idx)
    n_prots = len(prot_to_idx)
    print(f"  PPI adj shape:          {ppi.shape}  nnz={ppi.nnz:,}")
    print(f"  Drug-target adj shape:  {dt.shape}   nnz={dt.nnz:,}")
    print(f"  Drugs: {n_drugs:,}   Proteins: {n_prots:,}")
    return ppi, dt, drug_to_idx, prot_to_idx, d


# ---------------------------------------------------------------------------
# SE filtering — keep only side effects with enough supporting drug pairs
# ---------------------------------------------------------------------------

def filter_combo_se(
    combo_df: pd.DataFrame,
    idx_to_se: dict,
    min_se_pairs: int,
    top_se: int,
) -> tuple:
    """
    Filter combo side effects to those with >= min_se_pairs unique drug pairs.
    Then optionally truncate to the top_se most-connected SEs.

    Returns:
      filtered_combo_df   — combo_df restricted to kept SEs, se_idx re-indexed 0..n
      se_names            — list of MedDRA PT strings in new index order
    """
    # Count unique drug pairs per SE
    pair_counts = (
        combo_df
        .groupby("se_idx")[["drug_i_idx", "drug_j_idx"]]
        .apply(lambda g: g.drop_duplicates().shape[0])
        .rename("n_pairs")
        .reset_index()
    )
    kept = pair_counts[pair_counts["n_pairs"] >= min_se_pairs].copy()
    kept = kept.sort_values("n_pairs", ascending=False)

    if top_se and len(kept) > top_se:
        kept = kept.head(top_se)

    print(f"  SE after min_pairs>={min_se_pairs}: {len(kept):,}")
    if top_se:
        print(f"  SE after top_se={top_se}:          {len(kept):,}")

    kept_old_idx = set(kept["se_idx"])
    # Re-index SEs to dense 0..n-1
    old_to_new = {old: new for new, old in enumerate(kept["se_idx"].tolist())}

    filtered = combo_df[combo_df["se_idx"].isin(kept_old_idx)].copy()
    filtered["se_idx"] = filtered["se_idx"].map(old_to_new)

    se_names = [idx_to_se[old] for old in kept["se_idx"].tolist()]
    return filtered, se_names


def filter_mono_se(
    mono_df: pd.DataFrame,
    idx_to_se: dict,
    min_mono_drugs: int,
) -> tuple:
    """
    Keep mono side effects seen in >= min_mono_drugs drugs.
    Re-index to dense 0..n-1.
    """
    drug_counts = mono_df.groupby("se_idx")["drug_idx"].nunique().reset_index(name="n_drugs")
    kept = drug_counts[drug_counts["n_drugs"] >= min_mono_drugs]
    print(f"  Mono SE after min_drugs>={min_mono_drugs}: {len(kept):,}")

    kept_old_idx = set(kept["se_idx"])
    old_to_new = {old: new for new, old in enumerate(sorted(kept_old_idx))}

    filtered = mono_df[mono_df["se_idx"].isin(kept_old_idx)].copy()
    filtered["se_idx"] = filtered["se_idx"].map(old_to_new)

    se_names = [idx_to_se[old] for old in sorted(kept_old_idx)]
    return filtered, se_names


# ---------------------------------------------------------------------------
# Build adjacency structures
# ---------------------------------------------------------------------------

def build_combo_coo(combo_df: pd.DataFrame, n_drugs: int, n_se: int) -> dict:
    """
    Store combo data in COO form: three integer arrays + shape.

    This representation is compact and can reconstruct any per-SE sparse matrix
    on demand via load_combo_adjacencies().
    """
    drug_i = combo_df["drug_i_idx"].to_numpy(dtype=np.int32)
    drug_j = combo_df["drug_j_idx"].to_numpy(dtype=np.int32)
    se_idx = combo_df["se_idx"].to_numpy(dtype=np.int32)

    print(f"  Combo COO entries:  {len(drug_i):,}")
    print(f"  Unique drug pairs:  {len(set(zip(drug_i.tolist(), drug_j.tolist()))):,}")
    print(f"  Side effects:       {n_se:,}")
    return {"drug_i": drug_i, "drug_j": drug_j, "se_idx": se_idx,
            "n_drugs": np.array(n_drugs), "n_se": np.array(n_se)}


def build_mono_adj(mono_df: pd.DataFrame, n_drugs: int, n_se: int) -> sp.csr_matrix:
    """Build sparse (n_drugs × n_se_mono) indicator matrix."""
    rows = mono_df["drug_idx"].to_numpy(dtype=np.int32)
    cols = mono_df["se_idx"].to_numpy(dtype=np.int32)
    data = np.ones(len(rows), dtype=np.int8)
    M = sp.csr_matrix((data, (rows, cols)), shape=(n_drugs, n_se), dtype=np.int8)
    M.sum_duplicates()
    M.data[:] = 1
    print(f"  Mono adj shape: {M.shape}  nnz={M.nnz:,}")
    return M


def materialise_per_se_adjs(
    combo_df: pd.DataFrame,
    n_drugs: int,
    n_se: int,
    out_dir: str,
) -> None:
    """
    Write one sparse .npz per side effect under out_dir/combo_adjs/.
    Only called in strict mode where n_drugs is small (~3,471).
    """
    combo_adjs_dir = os.path.join(out_dir, "combo_adjs")
    os.makedirs(combo_adjs_dir, exist_ok=True)

    for se_i in range(n_se):
        sub = combo_df[combo_df["se_idx"] == se_i]
        rows = np.concatenate([sub["drug_i_idx"].values, sub["drug_j_idx"].values])
        cols = np.concatenate([sub["drug_j_idx"].values, sub["drug_i_idx"].values])
        data = np.ones(len(rows), dtype=np.int8)
        A = sp.csr_matrix((data, (rows, cols)), shape=(n_drugs, n_drugs), dtype=np.int8)
        A.sum_duplicates()
        A.data[:] = 1
        sp.save_npz(os.path.join(combo_adjs_dir, f"se_{se_i:05d}.npz"), A)

        if (se_i + 1) % 100 == 0 or se_i == n_se - 1:
            print(f"  Materialised {se_i+1:,}/{n_se:,} per-SE matrices", end="\r")
    print()
    print(f"  Saved {n_se:,} per-SE adjacency matrices to {combo_adjs_dir}/")


# ---------------------------------------------------------------------------
# Save all outputs
# ---------------------------------------------------------------------------

def save_outputs(
    out_dir: str,
    combo_coo: dict,
    mono_adj: sp.csr_matrix,
    ppi_adj: sp.csr_matrix,
    drug_target_adj: sp.csr_matrix,
    se_names: list,
    mono_se_names: list,
    meta: dict,
    drug_to_idx: dict,
    prot_to_idx: dict,
) -> None:
    os.makedirs(out_dir, exist_ok=True)

    # COO combo data
    np.savez_compressed(
        os.path.join(out_dir, "combo_adj_coo.npz"), **combo_coo
    )
    print(f"Saved: combo_adj_coo.npz")

    # Mono sparse matrix
    sp.save_npz(os.path.join(out_dir, "mono_adj.npz"), mono_adj)
    print(f"Saved: mono_adj.npz  shape={mono_adj.shape}")

    # PPI and drug-target (copy so this dir is self-contained)
    sp.save_npz(os.path.join(out_dir, "ppi_adj.npz"), ppi_adj)
    print(f"Saved: ppi_adj.npz   shape={ppi_adj.shape}")
    sp.save_npz(os.path.join(out_dir, "drug_target_adj.npz"), drug_target_adj)
    print(f"Saved: drug_target_adj.npz  shape={drug_target_adj.shape}")

    # ID mappings (copied so this dir is fully self-contained)
    with open(os.path.join(out_dir, "drug_to_idx.json"), "w", encoding="utf-8") as f:
        json.dump(drug_to_idx, f)
    print(f"Saved: drug_to_idx.json  ({len(drug_to_idx):,} entries)")

    with open(os.path.join(out_dir, "uniprot_to_idx.json"), "w", encoding="utf-8") as f:
        json.dump(prot_to_idx, f)
    print(f"Saved: uniprot_to_idx.json  ({len(prot_to_idx):,} entries)")

    # SE name lists
    with open(os.path.join(out_dir, "se_names.json"), "w", encoding="utf-8") as f:
        json.dump(se_names, f, indent=2)
    print(f"Saved: se_names.json  ({len(se_names):,} entries)")

    with open(os.path.join(out_dir, "mono_se_names.json"), "w", encoding="utf-8") as f:
        json.dump(mono_se_names, f, indent=2)
    print(f"Saved: mono_se_names.json  ({len(mono_se_names):,} entries)")

    with open(os.path.join(out_dir, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    print(f"Saved: meta.json")


# ---------------------------------------------------------------------------
# Public loader — use this from a training script
# ---------------------------------------------------------------------------

def load_combo_adjacencies(inputs_dir: str) -> list:
    """
    Reconstruct the list of per-SE sparse (n_drugs × n_drugs) matrices from
    the compact COO format saved by this script.

    Usage in a training script:
        from build_decagon_inputs import load_combo_adjacencies
        adj_list = load_combo_adjacencies("data/processed/decagon_inputs_all")
        # adj_list[i] is a (n_drugs × n_drugs) scipy csr_matrix

    For the broad drug set (672k drugs) each matrix is very sparse but
    still large; iterate rather than stacking.
    """
    coo = np.load(os.path.join(inputs_dir, "combo_adj_coo.npz"))
    drug_i = coo["drug_i"]
    drug_j = coo["drug_j"]
    se_idx = coo["se_idx"]
    n_drugs = int(coo["n_drugs"])
    n_se = int(coo["n_se"])

    adj_list = []
    ones = np.ones(2, dtype=np.int8)  # reused per-entry (avoid alloc in loop)

    for s in range(n_se):
        mask = se_idx == s
        rows = np.concatenate([drug_i[mask], drug_j[mask]])
        cols = np.concatenate([drug_j[mask], drug_i[mask]])
        data = np.ones(len(rows), dtype=np.int8)
        A = sp.csr_matrix((data, (rows, cols)), shape=(n_drugs, n_drugs), dtype=np.int8)
        A.sum_duplicates()
        A.data[:] = 1
        adj_list.append(A)

    return adj_list


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Assemble Decagon-format sparse tensors from processed artifacts."
    )
    parser.add_argument(
        "--mode", choices=["strict", "all"], default="all",
        help="Drug universe: 'strict' (3,471 drugs) or 'all' (672,337 drugs). Default: all.",
    )
    parser.add_argument(
        "--min-se-pairs", type=int, default=10,
        help=(
            "Min unique drug pairs a side effect must appear in to be included. "
            "Decagon paper used 500; start lower to assess your data. Default: 10."
        ),
    )
    parser.add_argument(
        "--top-se", type=int, default=None,
        help="Keep only the top N side effects by drug-pair count (optional cap).",
    )
    parser.add_argument(
        "--min-mono-drugs", type=int, default=3,
        help="Min drugs a mono side effect must appear in to be included. Default: 3.",
    )
    parser.add_argument(
        "--materialise", action="store_true",
        help=(
            "Write one .npz per side effect to combo_adjs/. "
            "Only practical for --mode strict (small drug count)."
        ),
    )
    args = parser.parse_args()

    out_subdir = "decagon_inputs" if args.mode == "strict" else "decagon_inputs_all"
    OUT_DIR = os.path.join(DATA_PROCESSED, out_subdir)

    print("=" * 60)
    print("Build Decagon Inputs")
    print("=" * 60)
    print(f"  Mode:           {args.mode}")
    print(f"  Min SE pairs:   {args.min_se_pairs}")
    print(f"  Top SE cap:     {args.top_se or 'none'}")
    print(f"  Min mono drugs: {args.min_mono_drugs}")
    print(f"  Materialise:    {args.materialise}")

    # ── 1. Load process_faers outputs ───────────────────────────────────────
    print("\n--- Step 1: Load FAERS side-effect outputs ---")
    combo_df, mono_df, idx_to_se = load_faers_outputs(FAERS_DIR)

    # ── 2. Load Decagon graph matrices ──────────────────────────────────────
    print("\n--- Step 2: Load Decagon graph artifacts ---")
    ppi_adj, drug_target_adj, drug_to_idx, prot_to_idx, _ = load_decagon_artifacts(args.mode)
    n_drugs = len(drug_to_idx)
    n_prots = len(prot_to_idx)

    # Guard: drug indices in combo/mono must be within the chosen drug universe
    if len(combo_df) > 0:
        max_drug_idx = max(combo_df["drug_i_idx"].max(), combo_df["drug_j_idx"].max())
        if max_drug_idx >= n_drugs:
            raise ValueError(
                f"combo_counts contains drug_idx={max_drug_idx} but "
                f"decagon drug universe only has {n_drugs} drugs. "
                f"Re-run process_faers.py with --mode {args.mode}."
            )

    # ── 3. Filter side effects ───────────────────────────────────────────────
    print("\n--- Step 3: Filter side effects ---")
    combo_df, se_names = filter_combo_se(
        combo_df, idx_to_se, args.min_se_pairs, args.top_se
    )
    mono_df, mono_se_names = filter_mono_se(mono_df, idx_to_se, args.min_mono_drugs)
    n_se = len(se_names)
    n_se_mono = len(mono_se_names)

    # ── 4. Build combo COO ──────────────────────────────────────────────────
    print("\n--- Step 4: Build combo COO structure ---")
    combo_coo = build_combo_coo(combo_df, n_drugs, n_se)

    # ── 5. Build mono adj ────────────────────────────────────────────────────
    print("\n--- Step 5: Build mono adjacency matrix ---")
    mono_adj = build_mono_adj(mono_df, n_drugs, n_se_mono)

    # ── 6. Optionally materialise per-SE matrices (strict only) ─────────────
    if args.materialise:
        if args.mode == "all":
            print(
                "\nWARNING: --materialise with --mode all requires storing "
                f"{n_se:,} sparse matrices of shape ({n_drugs:,} × {n_drugs:,}). "
                "This may use tens of GB. Proceeding anyway..."
            )
        print(f"\n--- Step 6: Materialise {n_se:,} per-SE adjacency matrices ---")
        materialise_per_se_adjs(combo_df, n_drugs, n_se, OUT_DIR)

    # ── 7. Save all outputs ──────────────────────────────────────────────────
    print("\n--- Step 7: Save outputs ---")
    meta = {
        "mode": args.mode,
        "min_se_pairs": args.min_se_pairs,
        "top_se": args.top_se,
        "min_mono_drugs": args.min_mono_drugs,
        "n_drugs": n_drugs,
        "n_proteins": n_prots,
        "n_combo_se": n_se,
        "n_mono_se": n_se_mono,
        "n_combo_triples": len(combo_df),
        "n_mono_pairs": len(mono_df),
        "ppi_shape": list(ppi_adj.shape),
        "ppi_nnz": int(ppi_adj.nnz),
        "drug_target_shape": list(drug_target_adj.shape),
        "drug_target_nnz": int(drug_target_adj.nnz),
        "adj_mats_layout": {
            "(0,0)": "ppi_adj.npz — protein×protein",
            "(0,1)": "drug_target_adj.npz transposed — protein×drug",
            "(1,0)": "drug_target_adj.npz — drug×protein",
            "(1,1)": "combo_adj_coo.npz — n_combo_se drug×drug matrices (COO)",
        },
    }
    save_outputs(
        OUT_DIR, combo_coo, mono_adj,
        ppi_adj, drug_target_adj,
        se_names, mono_se_names, meta,
        drug_to_idx, prot_to_idx,
    )

    print("\n" + "=" * 60)
    print("Decagon inputs assembled!")
    print("=" * 60)
    print(f"  Proteins:           {n_prots:,}")
    print(f"  Drugs:              {n_drugs:,}")
    print(f"  Combo side effects: {n_se:,}")
    print(f"  Mono  side effects: {n_se_mono:,}")
    print(f"  Output dir:         {OUT_DIR}")
    print()
    print("To reconstruct per-SE adjacency matrices in a training script:")
    print("  from src.build_decagon_inputs import load_combo_adjacencies")
    print(f"  adj_list = load_combo_adjacencies('{OUT_DIR}')")
    print("  # adj_list[i] -> scipy csr_matrix (n_drugs × n_drugs)")
