"""
Process STRING physical PPI into a UniProt-indexed sparse adjacency matrix.

Pipeline
--------
  1. load_mapping()    Build ENSP -> UniProt mapping from UniProt's bulk FTP
                       idmapping file (HUMAN_9606_idmapping.dat.gz).
                       Downloads the file on first run (~80 MB, one-time).
                       Prefers reviewed (Swiss-Prot) entries over TrEMBL.
  2. stream_edges()    Stream STRING physical-links file; map, filter, canonicalise
                       (i < j), and optionally deduplicate on the fly.
  3. build_sparse()    Expand canonical edges to both directions; build symmetric
                       CSR (binary int8 or weighted float32).
  4. save_artifacts()  Write npz, json, parquet, stats.json, meta.json.

Outputs  (data/processed/string/)
----------------------------------
  ppi_uniprot_csr.npz          Symmetric CSR adjacency in UniProt-index space.
  uniprot_to_idx.json          UniProt accession -> integer index (0-based).
  idx_to_uniprot.json          Reverse: list where [i] = UniProt accession.
  ppi_edges_idx.parquet        Canonical edges (i, j, weight), i < j, unique.
  ppi_edges_uniprot.parquet    Same edges with UniProt string IDs.
  reviewed_human_acs.txt       Cached list of Swiss-Prot human accessions.
  stats.json                   Full processing counters for reproducibility.
  meta.json                    Run configuration (paths, thresholds, flags).

ID-space note
-------------
  STRING uses Ensembl protein IDs (9606.ENSP...).
  ChEMBL uses UniProt accessions (e.g. P04637).
  We unify to UniProt so PPI nodes and drug-target nodes share one index space
  in build_decagon_graph.py.

  Mapping is done via UniProt's bulk FTP idmapping file which stores every
  Ensembl_PRO <-> UniProt cross-reference directly.  This avoids the UniProt
  async job API which proved unreliable for large batches.
"""

import gzip
import json
import os
import urllib.request
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# UniProt data sources
# ---------------------------------------------------------------------------

_IDMAPPING_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release"
    "/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"
)
_REVIEWED_ACS_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=organism_id%3A9606+AND+reviewed%3Atrue&format=list"
)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_RAW = os.path.join(PROJECT_ROOT, "data", "raw")
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


# ---------------------------------------------------------------------------
# Step 1: Load ENSP -> UniProt mapping via UniProt FTP idmapping file
# ---------------------------------------------------------------------------

def collect_ensp_ids(ppi_path: str, score_min: int) -> Set[str]:
    """First pass through the PPI file — collect ENSPs from score-passing edges."""
    ensp_set: Set[str] = set()
    with _open_text(ppi_path) as f:
        f.readline()  # header
        for line in f:
            parts = line.split()
            if len(parts) >= 3 and int(parts[2]) >= score_min:
                ensp_set.add(parts[0])
                ensp_set.add(parts[1])
    print(f"  Unique ENSPs with score >= {score_min}: {len(ensp_set):,}")
    return ensp_set


def fetch_reviewed_acs(cache_path: str) -> Set[str]:
    """
    Return the set of reviewed (Swiss-Prot) human UniProt accessions.

    Downloads once from the UniProt search API and caches to cache_path.
    The request is a simple synchronous GET — no async job needed.
    """
    if os.path.exists(cache_path):
        with open(cache_path, encoding="utf-8") as f:
            acs = {line.strip() for line in f if line.strip()}
        print(f"  Loaded {len(acs):,} reviewed ACs from cache")
        return acs

    print("  Fetching reviewed human UniProt accessions (Swiss-Prot)...")
    req = urllib.request.Request(_REVIEWED_ACS_URL, headers={"Accept": "text/plain"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        content = resp.read().decode("utf-8")

    acs = {line.strip() for line in content.splitlines() if line.strip()}
    print(f"  Found {len(acs):,} reviewed human accessions")

    with open(cache_path, "w", encoding="utf-8") as f:
        f.write("\n".join(sorted(acs)))
    print(f"  Cached: {cache_path}")
    return acs


def _download_idmapping(dest_path: str) -> None:
    """Download HUMAN_9606_idmapping.dat.gz with progress reporting."""
    print(f"  Source: {_IDMAPPING_URL}")
    print(f"  Dest:   {dest_path}")
    print("  Size:   ~80-100 MB (one-time download)")

    last_pct = [-1]

    def _progress(count, block_size, total_size):
        if total_size > 0:
            pct = min(100, count * block_size * 100 // total_size)
            if pct != last_pct[0] and pct % 10 == 0:
                print(f"    {pct}%")
                last_pct[0] = pct

    urllib.request.urlretrieve(_IDMAPPING_URL, dest_path, _progress)
    print("  Download complete.")


def load_mapping(
    idmapping_path: str,
    ensp_set: Set[str],
    reviewed_acs: Set[str],
) -> Tuple[Dict[str, str], dict]:
    """
    Build ENSP -> canonical UniProt mapping from UniProt's bulk idmapping file.

    File format (tab-separated, no header):
        UniProtAC  ID_type  ID
    We extract lines where ID_type == 'Ensembl_PRO'.  The ID is a bare ENSP
    (no '9606.' prefix); we match against ensp_set after stripping that prefix.

    Resolution rule per ENSP:
      1. Prefer reviewed (Swiss-Prot) UniProt ACs (present in reviewed_acs).
      2. Among reviewed candidates, keep the first seen in file order.
      3. If no reviewed entry exists, fall back to the first TrEMBL entry seen.

    Returns
    -------
    ensp_to_uniprot  {9606.ENSP... -> canonical UniProt accession}
    stats            Counters for stats.json.
    """
    bare_ensp_set = {e.removeprefix("9606.") for e in ensp_set}

    reviewed_map: Dict[str, str] = {}  # bare_ensp -> first reviewed AC
    trembl_map:   Dict[str, str] = {}  # bare_ensp -> first TrEMBL AC

    print(f"  Streaming: {idmapping_path}")
    n_lines = 0
    with gzip.open(idmapping_path, "rt", encoding="utf-8") as f:
        for line in f:
            n_lines += 1
            if n_lines % 5_000_000 == 0:
                print(f"    {n_lines:,} lines parsed...")

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            ac, id_type, bare_ensp = parts

            if id_type != "Ensembl_PRO":
                continue
            bare_ensp = bare_ensp.split(".")[0]  # strip version suffix e.g. ENSP00000269305.4 -> ENSP00000269305
            if bare_ensp not in bare_ensp_set:
                continue

            if ac in reviewed_acs:
                reviewed_map.setdefault(bare_ensp, ac)
            else:
                trembl_map.setdefault(bare_ensp, ac)

    print(f"  Parsed {n_lines:,} lines total")

    ensp_to_uniprot: Dict[str, str] = {}
    n_swiss = n_trembl = n_unmapped = 0

    for bare in bare_ensp_set:
        full_ensp = "9606." + bare
        if bare in reviewed_map:
            ensp_to_uniprot[full_ensp] = reviewed_map[bare]
            n_swiss += 1
        elif bare in trembl_map:
            ensp_to_uniprot[full_ensp] = trembl_map[bare]
            n_trembl += 1
        else:
            n_unmapped += 1

    print(f"  ENSP -> UniProt resolved: {len(ensp_to_uniprot):,}")
    print(f"    Swiss-Prot (reviewed): {n_swiss:,}")
    print(f"    TrEMBL (unreviewed):   {n_trembl:,}")
    print(f"    Unmapped:              {n_unmapped:,}")

    stats: dict = {
        "mapping_source":        "UniProt FTP idmapping (Ensembl_PRO)",
        "n_mappings":            len(ensp_to_uniprot),
        "n_mappings_swiss_prot": n_swiss,
        "n_mappings_trembl":     n_trembl,
        "n_mappings_unmapped":   n_unmapped,
    }
    return ensp_to_uniprot, stats


# ---------------------------------------------------------------------------
# Step 2: Stream edges from STRING PPI file
# ---------------------------------------------------------------------------

def stream_edges(
    ppi_path: str,
    ensp_to_uniprot: Dict[str, str],
    score_min: int = 700,
    max_edges: Optional[int] = None,
    drop_self_loops: bool = True,
    dedup_in_memory: bool = True,
    print_every_seen: int = 2_000_000,
) -> Tuple[List[Tuple[int, int, float, str, str]], Dict[str, int], dict]:
    """
    Stream STRING physical-links, map ENSP -> UniProt, filter, and collect
    canonical undirected edges (i < j).
    """
    uniprot_to_idx: Dict[str, int] = {}

    def _idx(uniprot: str) -> int:
        if uniprot not in uniprot_to_idx:
            uniprot_to_idx[uniprot] = len(uniprot_to_idx)
        return uniprot_to_idx[uniprot]

    seen_edges: Optional[Set[Tuple[int, int]]] = set() if dedup_in_memory else None
    edges: List[Tuple[int, int, float, str, str]] = []

    seen = kept = 0
    dropped_score = dropped_no_map = dropped_self_loop = dedup_skipped = 0

    ensp_observed: Set[str] = set()
    ensp_mapped: Set[str] = set()

    with _open_text(ppi_path) as f:
        header = f.readline().strip()
        print(f"PPI header: {header}")

        for line in f:
            seen += 1
            if print_every_seen and seen % print_every_seen == 0:
                print(
                    f"  seen={seen:,}  kept={kept:,}  "
                    f"drop_score={dropped_score:,}  drop_no_map={dropped_no_map:,}  "
                    f"drop_self={dropped_self_loop:,}  dedup_skip={dedup_skipped:,}  "
                    f"proteins={len(uniprot_to_idx):,}"
                )

            p1, p2, score_s = line.split()
            score = int(score_s)

            if score < score_min:
                dropped_score += 1
                continue

            ensp_observed.add(p1)
            ensp_observed.add(p2)

            u1 = ensp_to_uniprot.get(p1)
            u2 = ensp_to_uniprot.get(p2)
            if u1 is not None:
                ensp_mapped.add(p1)
            if u2 is not None:
                ensp_mapped.add(p2)
            if u1 is None or u2 is None:
                dropped_no_map += 1
                continue

            if drop_self_loops and u1 == u2:
                dropped_self_loop += 1
                continue

            i, j = _idx(u1), _idx(u2)

            if drop_self_loops and i == j:
                dropped_self_loop += 1
                continue

            if i > j:
                i, j = j, i
                u1, u2 = u2, u1

            if seen_edges is not None:
                if (i, j) in seen_edges:
                    dedup_skipped += 1
                    continue
                seen_edges.add((i, j))

            edges.append((i, j, score / 1000.0, u1, u2))
            kept += 1

            if max_edges is not None and kept >= max_edges:
                break

    print(
        f"PPI stream done: seen={seen:,}  kept={kept:,}  "
        f"drop_score={dropped_score:,}  drop_no_map={dropped_no_map:,}  "
        f"drop_self={dropped_self_loop:,}  dedup_skip={dedup_skipped:,}  "
        f"proteins={len(uniprot_to_idx):,}"
    )

    n_ensp_observed = len(ensp_observed)
    n_ensp_mapped   = len(ensp_mapped)

    stats: dict = {
        "ppi_seen":               seen,
        "ppi_pass_score":         seen - dropped_score,
        "ppi_pass_mapping":       seen - dropped_score - dropped_no_map - dropped_self_loop,
        "dropped_score":          dropped_score,
        "dropped_no_map":         dropped_no_map,
        "dropped_self_loop":      dropped_self_loop,
        "dedup_skipped":          dedup_skipped,
        "kept_edges":             kept,
        "n_ensp_observed_in_ppi": n_ensp_observed,
        "n_ensp_mapped_in_ppi":   n_ensp_mapped,
        "mapping_coverage_ppi":   n_ensp_mapped / n_ensp_observed if n_ensp_observed else 0.0,
        "n_uniprot_in_graph":     len(uniprot_to_idx),
    }
    return edges, uniprot_to_idx, stats


# ---------------------------------------------------------------------------
# Step 3: Build symmetric sparse adjacency
# ---------------------------------------------------------------------------

def build_sparse(
    edges: List[Tuple[int, int, float, str, str]],
    uniprot_to_idx: Dict[str, int],
    use_weights: bool = False,
) -> sp.csr_matrix:
    """
    Build a symmetric CSR adjacency from canonical (i < j) edges.

    use_weights=False  binary int8
    use_weights=True   float32 STRING score/1000, max-coalesced on duplicates
    """
    n = len(uniprot_to_idx)

    if not edges:
        dtype = np.float32 if use_weights else np.int8
        return sp.csr_matrix((n, n), dtype=dtype)

    if use_weights:
        weight_dict: Dict[Tuple[int, int], np.float32] = {}
        for i, j, norm_score, *_ in edges:
            val = np.float32(norm_score)
            if val > weight_dict.get((i, j), np.float32(0.0)):
                weight_dict[(i, j)] = val

        ci_arr = np.array([k[0] for k in weight_dict], dtype=np.int32)
        cj_arr = np.array([k[1] for k in weight_dict], dtype=np.int32)
        w_arr  = np.array(list(weight_dict.values()), dtype=np.float32)

        rows_arr = np.concatenate([ci_arr, cj_arr])
        cols_arr = np.concatenate([cj_arr, ci_arr])
        data_arr = np.concatenate([w_arr,  w_arr])

        A = sp.csr_matrix((data_arr, (rows_arr, cols_arr)), shape=(n, n))
        A.sum_duplicates()
        np.clip(A.data, 0.0, 1.0, out=A.data)

    else:
        ne = len(edges)
        rows_arr = np.empty(2 * ne, dtype=np.int32)
        cols_arr = np.empty(2 * ne, dtype=np.int32)
        data_arr = np.ones(2 * ne, dtype=np.int8)

        for k, (i, j, *_) in enumerate(edges):
            rows_arr[2 * k],     cols_arr[2 * k]     = i, j
            rows_arr[2 * k + 1], cols_arr[2 * k + 1] = j, i

        A = sp.csr_matrix((data_arr, (rows_arr, cols_arr)), shape=(n, n))
        A.sum_duplicates()
        A.data[:] = 1

    return A


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def degrees_from_adj(A: sp.csr_matrix) -> np.ndarray:
    return np.asarray(A.sum(axis=1)).ravel()


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_artifacts(
    edges: List[Tuple[int, int, float, str, str]],
    A: sp.csr_matrix,
    uniprot_to_idx: Dict[str, int],
    use_weights: bool = False,
    drop_self_loops: bool = True,
) -> None:
    n = len(uniprot_to_idx)
    print("Validating artifacts...")

    for k, (i, j, w, *_) in enumerate(edges):
        assert 0 <= i < n and 0 <= j < n, f"Edge #{k}: index out of range i={i} j={j} n={n}"
        assert i < j, f"Edge #{k}: not canonical i={i} j={j}"
        assert np.isfinite(w), f"Edge #{k}: non-finite weight w={w}"
        if use_weights:
            assert 0.0 < w <= 1.0, f"Edge #{k}: weight outside (0,1] w={w}"

    diff = A - A.T
    diff.eliminate_zeros()
    assert diff.nnz == 0, f"Adjacency not symmetric: (A-A.T).nnz={diff.nnz}"

    if drop_self_loops:
        diag_nnz = int((A.diagonal() != 0).sum())
        assert diag_nnz == 0, f"Adjacency has {diag_nnz} self-loop entries"

    print(
        f"  OK: {len(edges):,} edges  n={n:,}  symmetric  "
        + ("no diagonal  " if drop_self_loops else "")
        + ("weights in (0,1]" if use_weights else "binary")
    )


# ---------------------------------------------------------------------------
# Step 4: Save all artifacts
# ---------------------------------------------------------------------------

def save_artifacts(
    out_dir: str,
    A: sp.csr_matrix,
    uniprot_to_idx: Dict[str, int],
    edges: List[Tuple[int, int, float, str, str]],
    stats: dict,
    meta: dict,
) -> None:
    os.makedirs(out_dir, exist_ok=True)

    npz_path = os.path.join(out_dir, "ppi_uniprot_csr.npz")
    sp.save_npz(npz_path, A)
    print(f"Saved: {npz_path}  shape={A.shape}  nnz={A.nnz:,}")

    idx_path = os.path.join(out_dir, "uniprot_to_idx.json")
    with open(idx_path, "w", encoding="utf-8") as f:
        json.dump(uniprot_to_idx, f)
    print(f"Saved: {idx_path}  ({len(uniprot_to_idx):,} entries)")

    idx_to_uniprot: List[str] = [""] * len(uniprot_to_idx)
    for acc, idx in uniprot_to_idx.items():
        idx_to_uniprot[idx] = acc
    itu_path = os.path.join(out_dir, "idx_to_uniprot.json")
    with open(itu_path, "w", encoding="utf-8") as f:
        json.dump(idx_to_uniprot, f)
    print(f"Saved: {itu_path}  ({len(idx_to_uniprot):,} entries)")

    edges_df = (
        pd.DataFrame(
            [(i, j, w) for i, j, w, *_ in edges],
            columns=["i", "j", "weight"],
        )
        .groupby(["i", "j"], as_index=False)["weight"]
        .max()
    )
    edges_df["uniprot_1"] = [idx_to_uniprot[i] for i in edges_df["i"]]
    edges_df["uniprot_2"] = [idx_to_uniprot[j] for j in edges_df["j"]]
    n_canonical = len(edges_df)

    edges_idx_path = os.path.join(out_dir, "ppi_edges_idx.parquet")
    (
        edges_df[["i", "j", "weight"]]
        .astype({"i": "int32", "j": "int32", "weight": "float32"})
        .to_parquet(edges_idx_path, index=False)
    )
    print(f"Saved: {edges_idx_path}  ({n_canonical:,} canonical edges)")

    edges_uniprot_path = os.path.join(out_dir, "ppi_edges_uniprot.parquet")
    edges_df[["uniprot_1", "uniprot_2", "weight"]].to_parquet(
        edges_uniprot_path, index=False
    )
    print(f"Saved: {edges_uniprot_path}  ({n_canonical:,} canonical edges)")

    deg = degrees_from_adj(A)
    stats_full: dict = {
        **stats,
        "n_proteins":        int(A.shape[0]),
        "nnz":               int(A.nnz),
        "n_canonical_edges": n_canonical,
        "degree_min":        float(deg.min()),
        "degree_mean":       float(deg.mean()),
        "degree_max":        float(deg.max()),
    }
    stats_path = os.path.join(out_dir, "stats.json")
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats_full, f, indent=2)
    print(f"Saved: {stats_path}")

    meta_path = os.path.join(out_dir, "meta.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    print(f"Saved: {meta_path}")


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys as _sys
    _src = os.path.dirname(os.path.abspath(__file__))
    if _src not in _sys.path:
        _sys.path.insert(0, _src)
    import config as _cfg

    PPI_PATH       = _cfg.STRING_PPI_PATH
    IDMAPPING_PATH = _cfg.UNIPROT_IDMAPPING_PATH
    OUT_DIR        = os.path.join(DATA_PROCESSED, "string")
    REVIEWED_CACHE = os.path.join(OUT_DIR, "reviewed_human_acs.txt")

    # ── Configuration ────────────────────────────────────────────────
    SCORE_MIN       = 700    # high-confidence physical interactions (≈10k proteins, ≈80k edges)
    USE_WEIGHTS     = False
    DEDUP_IN_MEMORY = True
    MAX_EDGES       = None

    print("=" * 60)
    print("Step 1: Build ENSP -> UniProt mapping")
    print("=" * 60)
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Collecting ENSPs from PPI file (first pass)...")
    ensp_set = collect_ensp_ids(PPI_PATH, SCORE_MIN)

    print("Loading reviewed (Swiss-Prot) human accessions...")
    reviewed_acs = fetch_reviewed_acs(REVIEWED_CACHE)

    if not os.path.exists(IDMAPPING_PATH):
        os.makedirs(os.path.dirname(IDMAPPING_PATH), exist_ok=True)
        print("Downloading UniProt idmapping file (one-time, ~80 MB)...")
        _download_idmapping(IDMAPPING_PATH)

    ensp_to_uniprot, mapping_stats = load_mapping(
        IDMAPPING_PATH, ensp_set, reviewed_acs
    )

    print("\n" + "=" * 60)
    print("Step 2: Stream STRING PPI edges")
    print("=" * 60)
    edges, uniprot_to_idx, edge_stats = stream_edges(
        PPI_PATH,
        ensp_to_uniprot=ensp_to_uniprot,
        score_min=SCORE_MIN,
        max_edges=MAX_EDGES,
        drop_self_loops=True,
        dedup_in_memory=DEDUP_IN_MEMORY,
        print_every_seen=2_000_000,
    )

    print("\n" + "=" * 60)
    print("Step 3: Build sparse adjacency")
    print("=" * 60)
    A = build_sparse(edges, uniprot_to_idx, use_weights=USE_WEIGHTS)
    print(f"Adjacency: shape={A.shape}  nnz={A.nnz:,}  dtype={A.dtype}")
    deg = degrees_from_adj(A)
    print(f"Degree: min={deg.min():.0f}  mean={deg.mean():.2f}  max={deg.max():.0f}")

    print("\n" + "=" * 60)
    print("Step 3b: Validate artifacts")
    print("=" * 60)
    validate_artifacts(edges, A, uniprot_to_idx, use_weights=USE_WEIGHTS, drop_self_loops=True)

    print("\n" + "=" * 60)
    print("Step 4: Save artifacts")
    print("=" * 60)
    meta = {
        "ppi_path":        PPI_PATH,
        "idmapping_path":  IDMAPPING_PATH,
        "mapping_source":  "UniProt FTP idmapping (Ensembl_PRO -> UniProtKB)",
        "score_min":       SCORE_MIN,
        "use_weights":     USE_WEIGHTS,
        "dedup_in_memory": DEDUP_IN_MEMORY,
        "max_edges":       MAX_EDGES,
        "drop_self_loops": True,
    }
    combined_stats = {**mapping_stats, **edge_stats}
    save_artifacts(OUT_DIR, A, uniprot_to_idx, edges, combined_stats, meta)

    print("\nDone.")
