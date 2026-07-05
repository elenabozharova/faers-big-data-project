"""
Build unified Decagon-ready graph from independently processed pipelines.

Loads STRING PPI and ChEMBL drug-target outputs, unions the protein
ID spaces, and produces re-indexed adjacency matrices that share a single
consistent protein index.

Uses UNION (not intersection) of protein sets so the full PPI graph is
preserved for GCN message passing. ChEMBL-only targets become isolated
protein nodes (no PPI neighbors, but reachable via drug-target edges).

Run AFTER process_string.py and process_chembl.py have produced their outputs.
"""

import json
import os
from typing import Dict, List, Set, Tuple

import numpy as np
import scipy.sparse as sp

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")


# ── loaders ──────────────────────────────────────────────────────────

def load_string_artifacts(string_dir: str) -> Tuple[sp.csr_matrix, Dict[str, int]]:
    A = sp.load_npz(os.path.join(string_dir, "ppi_uniprot_csr.npz"))
    with open(os.path.join(string_dir, "uniprot_to_idx.json"), encoding="utf-8") as f:
        uniprot_to_idx = json.load(f)
    return A, uniprot_to_idx


def load_chembl_artifacts(chembl_dir: str) -> Tuple[List[Tuple[str, str]], Dict[str, int]]:
    pairs: List[Tuple[str, str]] = []
    with open(os.path.join(chembl_dir, "drug_target_pairs.tsv"), encoding="utf-8") as f:
        f.readline()  # skip header
        for line in f:
            drug, target = line.strip().split("\t")
            pairs.append((drug, target))

    with open(os.path.join(chembl_dir, "drug_to_idx.json"), encoding="utf-8") as f:
        drug_to_idx = json.load(f)

    return pairs, drug_to_idx


# ── union ────────────────────────────────────────────────────────────

def union_proteins(
    ppi_proteins: Set[str],
    chembl_targets: Set[str],
) -> Dict[str, int]:
    all_proteins = sorted(ppi_proteins | chembl_targets)
    unified = {acc: i for i, acc in enumerate(all_proteins)}

    intersection = ppi_proteins & chembl_targets
    ppi_only = ppi_proteins - chembl_targets
    chembl_only = chembl_targets - ppi_proteins

    print(f"PPI proteins:       {len(ppi_proteins):,}")
    print(f"ChEMBL targets:     {len(chembl_targets):,}")
    print(f"Union (total):      {len(unified):,}")
    print(f"  intersection:     {len(intersection):,}")
    print(f"  PPI-only:         {len(ppi_only):,}")
    print(f"  ChEMBL-only:      {len(chembl_only):,} (isolated in PPI)")

    return unified


# ── re-index PPI ─────────────────────────────────────────────────────

def reindex_ppi(
    A_orig: sp.csr_matrix,
    old_idx: Dict[str, int],
    unified_idx: Dict[str, int],
) -> sp.csr_matrix:
    old_to_new = {}
    for acc, new_i in unified_idx.items():
        if acc in old_idx:
            old_to_new[old_idx[acc]] = new_i

    A_coo = A_orig.tocoo()
    rows, cols, data = [], [], []

    for r, c, v in zip(A_coo.row, A_coo.col, A_coo.data):
        rows.append(old_to_new[r])
        cols.append(old_to_new[c])
        data.append(v)

    n = len(unified_idx)
    A_new = sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=np.int8)
    A_new.sum_duplicates()
    A_new.data[:] = 1

    print(f"\nPPI re-indexing:")
    print(f"  original nnz:  {A_orig.nnz:,}")
    print(f"  new nnz:       {A_new.nnz:,}")
    print(f"  new shape:     {A_new.shape}")

    return A_new


# ── build drug-target matrix ────────────────────────────────────────

def build_drug_target_adj(
    pairs: List[Tuple[str, str]],
    drug_to_idx: Dict[str, int],
    unified_protein_idx: Dict[str, int],
) -> sp.csr_matrix:
    rows, cols, data = [], [], []

    for drug, target in pairs:
        di = drug_to_idx[drug]
        pi = unified_protein_idx[target]  # all targets are in the union
        rows.append(di)
        cols.append(pi)
        data.append(1)

    n_drugs = len(drug_to_idx)
    n_proteins = len(unified_protein_idx)
    DT = sp.csr_matrix((data, (rows, cols)), shape=(n_drugs, n_proteins), dtype=np.int8)
    DT.sum_duplicates()
    DT.data[:] = 1

    print(f"\nDrug-target matrix:")
    print(f"  pairs:   {len(pairs):,}")
    print(f"  nnz:     {DT.nnz:,}")
    print(f"  shape:   {DT.shape}")
    print(f"  drugs:   {n_drugs:,}")

    return DT


# ── save ─────────────────────────────────────────────────────────────

def save_unified(
    out_dir: str,
    A_ppi: sp.csr_matrix,
    DT: sp.csr_matrix,
    unified_protein_idx: Dict[str, int],
    drug_to_idx: Dict[str, int],
    meta: dict,
) -> None:
    os.makedirs(out_dir, exist_ok=True)

    sp.save_npz(os.path.join(out_dir, "ppi_adj.npz"), A_ppi)
    sp.save_npz(os.path.join(out_dir, "drug_target_adj.npz"), DT)

    with open(os.path.join(out_dir, "uniprot_to_idx.json"), "w", encoding="utf-8") as f:
        json.dump(unified_protein_idx, f)

    with open(os.path.join(out_dir, "drug_to_idx.json"), "w", encoding="utf-8") as f:
        json.dump(drug_to_idx, f)

    with open(os.path.join(out_dir, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print(f"\nSaved unified artifacts to {out_dir}/")


# ── main ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    # Default to chembl_all (broad); pass "strict" as arg to use mechanism-only
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode == "strict":
        chembl_subdir = "chembl"
        out_subdir = "decagon"
    else:
        chembl_subdir = "chembl_all"
        out_subdir = "decagon_all"

    STRING_DIR = os.path.join(DATA_PROCESSED, "string")
    CHEMBL_DIR = os.path.join(DATA_PROCESSED, chembl_subdir)
    OUT_DIR = os.path.join(DATA_PROCESSED, out_subdir)

    print("=" * 60)
    print("Building unified Decagon graph")
    print("=" * 60)

    # 1. Load both pipeline outputs
    A_ppi_orig, ppi_uniprot_to_idx = load_string_artifacts(STRING_DIR)
    chembl_pairs, drug_to_idx = load_chembl_artifacts(CHEMBL_DIR)

    ppi_proteins = set(ppi_uniprot_to_idx.keys())
    chembl_targets = {target for _, target in chembl_pairs}

    # 2. Union protein sets (keep full PPI for message passing)
    unified_protein_idx = union_proteins(ppi_proteins, chembl_targets)

    # 3. Re-index PPI adjacency into unified space
    A_ppi = reindex_ppi(A_ppi_orig, ppi_uniprot_to_idx, unified_protein_idx)

    # 4. Build drug-target matrix
    DT = build_drug_target_adj(chembl_pairs, drug_to_idx, unified_protein_idx)

    # 5. Save
    intersection = ppi_proteins & chembl_targets
    meta = {
        "ppi_proteins_original": len(ppi_proteins),
        "chembl_targets_original": len(chembl_targets),
        "proteins_union": len(unified_protein_idx),
        "proteins_intersection": len(intersection),
        "proteins_ppi_only": len(ppi_proteins - chembl_targets),
        "proteins_chembl_only": len(chembl_targets - ppi_proteins),
        "ppi_nnz": int(A_ppi.nnz),
        "ppi_shape": list(A_ppi.shape),
        "drug_target_shape": list(DT.shape),
        "drug_target_nnz": int(DT.nnz),
        "n_drugs": int(DT.shape[0]),
        "n_proteins": len(unified_protein_idx),
    }
    save_unified(OUT_DIR, A_ppi, DT, unified_protein_idx, drug_to_idx, meta)

    print("\n" + "=" * 60)
    print("Done. Unified graph ready for Decagon.")
    print("=" * 60)
