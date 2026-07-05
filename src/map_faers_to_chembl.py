"""
Map FAERS drug names to ChEMBL molecule IDs.

Tiered matching strategy
------------------------
Tier 0  Field selection   prefer prod_ai (active ingredient); fall back to drugname.
Tier 1  Normalization     via normalize_drug_names.normalize_drug_string.
Tier 2a exact_pref        normalized name == normalized ChEMBL pref_name.
Tier 2b exact_synonym     normalized name == normalized ChEMBL synonym.
Tier 2c combo_all         every "|"-separated component individually matches
                          (exact_pref or exact_synonym), and none is ambiguous.
Tier 3  fuzzy             token-overlap prefilter + SequenceMatcher ratio >= threshold.
                          Applied only to names still unmatched after Tier 2.
                          Every fuzzy match is logged with match_score and
                          match_method="fuzzy" so it can be excluded from strict runs.

Ambiguity rule
--------------
  candidate_count == 1  -> chembl_id is set; is_ambiguous = False
  candidate_count  > 1  -> chembl_id = None;  is_ambiguous = True

Outputs (data/processed/faers/)
--------------------------------
  chembl_name_index.parquet   cached normalized ChEMBL name -> ID lookup (rebuild with --force-rebuild-index)
  faers_drug_to_chembl.parquet  one row per unique normalized FAERS drug name
  unmatched_drugs.csv           unmatched names sorted by FAERS frequency
  ambiguous_drugs.csv           ambiguous names sorted by FAERS frequency
  mapping_stats.json            validation summary (counts, rates, top-20 lists)

Usage
-----
  PGPASSWORD=secret python src/map_faers_to_chembl.py
  PGPASSWORD=secret python src/map_faers_to_chembl.py --no-fuzzy
  PGPASSWORD=secret python src/map_faers_to_chembl.py --max-rows 50000   # test mode
  PGPASSWORD=secret python src/map_faers_to_chembl.py --force-rebuild-index

Run AFTER process_string.py and process_chembl.py (ChEMBL DB must be up).
"""

import json
import os
import sys
import zipfile
from collections import defaultdict
from difflib import SequenceMatcher
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import psycopg2
from psycopg2.extras import RealDictCursor

# Resolve the src/ directory so that normalize_drug_names is importable
# regardless of the working directory the user runs from.
_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from normalize_drug_names import is_combo, normalize_drug_string, split_combo
import config as _cfg

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(_SRC_DIR)
DATA_RAW = os.path.join(PROJECT_ROOT, "data", "raw")
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")

FAERS_ZIP        = _cfg.FAERS_ZIP_PATH
FAERS_DRUG_INNER = _cfg.FAERS_DRUG_INNER

OUT_DIR = os.path.join(DATA_PROCESSED, "faers")


# ---------------------------------------------------------------------------
# Tier 0 + 1: Load FAERS and build unique normalized name table
# ---------------------------------------------------------------------------

def load_faers_drug_table(zip_path: str, inner_path: str) -> pd.DataFrame:
    """Load the FAERS DRUG table from the quarterly ZIP archive."""
    with zipfile.ZipFile(zip_path) as z:
        with z.open(inner_path) as f:
            df = pd.read_csv(f, sep="$", dtype=str, low_memory=False)
    print(f"FAERS DRUG table loaded: {len(df):,} rows, {df.shape[1]} columns")
    print(f"  Columns: {list(df.columns)}")
    return df


def build_unique_name_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply Tier 0 (field selection) and Tier 1 (normalization), then
    deduplicate by normalized form so each unique drug name is matched
    only once.

    Field selection (Tier 0):
      Use prod_ai when it is present, non-empty, has >2 characters, and is
      not purely numeric. Otherwise fall back to drugname.

    Returns a DataFrame with columns:
      raw_name      selected raw string (prod_ai or drugname)
      raw_ai        prod_ai value (may be empty)
      field_used    "prod_ai" or "drugname"
      normalized    output of normalize_drug_string(raw_name)
      faers_count   number of DRUG-table rows that share this normalized form
    """
    # Vectorized field selection
    prod_ai = (
        df["prod_ai"].fillna("").astype(str).str.strip()
        if "prod_ai" in df.columns
        else pd.Series([""] * len(df), dtype=str)
    )
    drugname = (
        df["drugname"].fillna("").astype(str).str.strip()
        if "drugname" in df.columns
        else pd.Series([""] * len(df), dtype=str)
    )

    # prod_ai is "meaningful" when it has >2 chars and is not purely numeric
    is_meaningful = (prod_ai.str.len() > 2) & (
        ~prod_ai.str.replace(".", "", regex=False).str.isnumeric()
    )

    raw_name = prod_ai.where(is_meaningful, drugname)
    field_used = pd.Series("prod_ai", index=df.index).where(is_meaningful, "drugname")

    names_df = pd.DataFrame({
        "raw_name": raw_name,
        "raw_ai": prod_ai,
        "field_used": field_used,
    })

    # Normalize (Tier 1)
    names_df["normalized"] = names_df["raw_name"].apply(normalize_drug_string)

    # Frequency count: how many DRUG-table rows map to each normalized form
    freq = names_df.groupby("normalized").size().rename("faers_count")

    # Deduplicate by normalized form; keep the first occurrence's raw info
    unique_df = (
        names_df.drop_duplicates(subset=["normalized"])
        .join(freq, on="normalized")
        .reset_index(drop=True)
    )

    n_prod_ai = is_meaningful.sum()
    n_empty = (unique_df["normalized"] == "").sum()

    print(f"\nFAERS unique name table (after Tier 0 + 1):")
    print(f"  Total DRUG rows:          {len(df):,}")
    print(f"  Rows using prod_ai:       {n_prod_ai:,} ({100 * n_prod_ai / len(df):.1f}%)")
    print(f"  Rows using drugname:      {len(df) - n_prod_ai:,}")
    print(f"  Unique raw_name values:   {names_df['raw_name'].nunique():,}")
    print(f"  Unique normalized forms:  {unique_df['normalized'].nunique():,}")
    print(f"  Empty after normalization:{n_empty:,}")

    return unique_df


# ---------------------------------------------------------------------------
# Build ChEMBL name index (Tier 2)
# ---------------------------------------------------------------------------

_CHEMBL_INDEX_QUERY = """
SELECT
    md.chembl_id,
    md.pref_name,
    md.molecule_type,
    ms.synonyms   AS synonym,
    ms.syn_type
FROM molecule_dictionary md
LEFT JOIN molecule_synonyms ms ON md.molregno = ms.molregno
WHERE md.molecule_type NOT IN ('Unknown')
  AND (md.pref_name IS NOT NULL OR ms.synonyms IS NOT NULL)
ORDER BY md.chembl_id
"""


def build_chembl_name_index(
    conn: psycopg2.extensions.connection,
    out_dir: str,
    force_rebuild: bool = False,
) -> Tuple[Dict[str, List[dict]], Dict[str, List[dict]]]:
    """
    Build two normalized lookup dicts from ChEMBL molecule names and synonyms.

    pref_index  normalized_pref_name  -> [{"chembl_id": ..., "original": ...}, ...]
    syn_index   normalized_synonym    -> [{"chembl_id": ..., "original": ...,
                                           "syn_type": ...}, ...]

    A parquet cache is written to out_dir/chembl_name_index.parquet so that
    the database does not need to be queried on every run.  Pass
    force_rebuild=True to re-query and overwrite the cache.
    """
    cache_path = os.path.join(out_dir, "chembl_name_index.parquet")

    if not force_rebuild and os.path.exists(cache_path):
        print(f"Loading cached ChEMBL name index from {cache_path}")
        raw = pd.read_parquet(cache_path)
        print(f"  Loaded {len(raw):,} rows from cache")
    else:
        print("Querying ChEMBL for molecule names and synonyms...")
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(_CHEMBL_INDEX_QUERY)
            rows = cur.fetchall()

        raw = pd.DataFrame(rows)
        print(f"  Fetched {len(raw):,} rows from ChEMBL")

        os.makedirs(out_dir, exist_ok=True)
        raw.to_parquet(cache_path, index=False)
        print(f"  Cached to {cache_path}")

    print("Building normalized name lookup dicts...")

    pref_index: Dict[str, List[dict]] = defaultdict(list)
    syn_index: Dict[str, List[dict]] = defaultdict(list)

    # Track ChEMBL IDs already added to pref_index to avoid duplicates
    # when the LEFT JOIN produces multiple synonym rows per molecule.
    pref_seen: Set[str] = set()

    for _, row in raw.iterrows():
        chembl_id: str = row["chembl_id"]
        pref_name: str = row.get("pref_name") or ""
        synonym: str = row.get("synonym") or ""
        syn_type: str = row.get("syn_type") or ""

        if pref_name and chembl_id not in pref_seen:
            norm_pref = normalize_drug_string(pref_name)
            if norm_pref:
                pref_index[norm_pref].append({
                    "chembl_id": chembl_id,
                    "original": pref_name,
                })
            pref_seen.add(chembl_id)

        if synonym:
            norm_syn = normalize_drug_string(synonym)
            if norm_syn:
                syn_index[norm_syn].append({
                    "chembl_id": chembl_id,
                    "original": synonym,
                    "syn_type": syn_type,
                })

    print(f"  pref_index unique keys: {len(pref_index):,}")
    print(f"  syn_index  unique keys: {len(syn_index):,}")

    return dict(pref_index), dict(syn_index)


# ---------------------------------------------------------------------------
# Token inverted index for fuzzy prefiltering (Tier 3)
# ---------------------------------------------------------------------------

def build_token_index(
    pref_index: Dict[str, List[dict]],
    syn_index: Dict[str, List[dict]],
) -> Dict[str, Set[str]]:
    """
    Build an inverted index:  token -> set of normalized ChEMBL names
    that contain that token.

    This dramatically reduces the candidate pool for fuzzy matching;
    only normalized names sharing at least one token with the query
    are scored with SequenceMatcher.
    """
    token_index: Dict[str, Set[str]] = defaultdict(set)
    for norm_name in list(pref_index.keys()) + list(syn_index.keys()):
        for token in norm_name.split():
            token_index[token].add(norm_name)
    print(f"  token_index unique tokens: {len(token_index):,}")
    return dict(token_index)


# ---------------------------------------------------------------------------
# Tier 2: exact matching helpers
# ---------------------------------------------------------------------------

def _resolve_candidates(
    candidates: List[dict],
) -> Tuple[Optional[str], bool, List[str]]:
    """
    Given a list of candidate dicts (each with a 'chembl_id' key),
    return (chembl_id, is_ambiguous, all_candidate_ids).

    If exactly one unique ChEMBL ID is found -> chembl_id set, is_ambiguous=False.
    If multiple unique ChEMBL IDs -> chembl_id=None, is_ambiguous=True.
    """
    chembl_ids = list(dict.fromkeys(c["chembl_id"] for c in candidates))  # preserve order, dedup
    if len(chembl_ids) == 1:
        return chembl_ids[0], False, chembl_ids
    return None, True, chembl_ids


def match_single_name(
    normalized: str,
    pref_index: Dict[str, List[dict]],
    syn_index: Dict[str, List[dict]],
) -> Optional[dict]:
    """
    Attempt Tier 2a (exact_pref) then Tier 2b (exact_synonym) for a
    single normalized drug name.  Returns a result dict or None.
    """
    # Tier 2a: preferred name match (highest confidence)
    if normalized in pref_index:
        chembl_id, is_ambiguous, candidates = _resolve_candidates(pref_index[normalized])
        return {
            "chembl_id": chembl_id,
            "match_method": "exact_pref",
            "match_score": 1.0,
            "candidate_count": len(candidates),
            "is_ambiguous": is_ambiguous,
            "candidates": candidates,
        }

    # Tier 2b: synonym match
    if normalized in syn_index:
        chembl_id, is_ambiguous, candidates = _resolve_candidates(syn_index[normalized])
        return {
            "chembl_id": chembl_id,
            "match_method": "exact_synonym",
            "match_score": 1.0,
            "candidate_count": len(candidates),
            "is_ambiguous": is_ambiguous,
            "candidates": candidates,
        }

    return None


def match_combo_name(
    normalized: str,
    pref_index: Dict[str, List[dict]],
    syn_index: Dict[str, List[dict]],
) -> Optional[dict]:
    """
    Tier 2c: match a combo drug by matching each "|"-separated component
    individually.  Only succeeds when EVERY component matches exactly and
    non-ambiguously.  Returns a result dict or None.

    The chembl_id field is set to None for combos (a combo product is
    represented by its component IDs in the 'candidates' list).
    Component ChEMBL IDs are stored in candidates and match_detail.
    """
    components = split_combo(normalized)
    if len(components) < 2:
        return None

    component_results = []
    for comp in components:
        m = match_single_name(comp, pref_index, syn_index)
        if m is None or m["is_ambiguous"]:
            return None  # any unmatched or ambiguous component fails the whole combo
        component_results.append((comp, m["chembl_id"], m["match_method"]))

    component_ids = [cid for _, cid, _ in component_results]
    component_methods = list({method for _, _, method in component_results})

    return {
        "chembl_id": None,  # no single ID for a multi-component combo
        "match_method": "combo_all",
        "match_score": 1.0,
        "candidate_count": len(component_ids),
        "is_ambiguous": False,
        "candidates": component_ids,
        "match_detail": json.dumps({
            comp: cid for comp, cid, _ in component_results
        }),
        "component_methods": component_methods,
    }


# ---------------------------------------------------------------------------
# Tier 3: fuzzy matching
# ---------------------------------------------------------------------------

def match_fuzzy(
    normalized: str,
    pref_index: Dict[str, List[dict]],
    syn_index: Dict[str, List[dict]],
    token_index: Dict[str, Set[str]],
    min_score: float = 0.85,
) -> Optional[dict]:
    """
    Tier 3: fuzzy match using token-overlap prefilter + SequenceMatcher.
    Applied only to names still unmatched after Tier 2.

    Steps:
      1. Find ChEMBL normalized names sharing >= 1 token with the query.
      2. Score each candidate with SequenceMatcher.ratio().
      3. Accept the best candidate above min_score; resolve ambiguity.

    Every fuzzy match is tagged match_method="fuzzy" with the similarity
    score recorded, so it can be excluded from strict downstream pipelines.
    """
    query_tokens = set(normalized.split())
    if not query_tokens:
        return None

    # Prefilter: only names that share at least one token
    candidate_names: Set[str] = set()
    for token in query_tokens:
        candidate_names.update(token_index.get(token, set()))

    if not candidate_names:
        return None

    # Score and rank
    scored: List[Tuple[float, str]] = []
    for cand_name in candidate_names:
        score = SequenceMatcher(None, normalized, cand_name).ratio()
        if score >= min_score:
            scored.append((score, cand_name))

    if not scored:
        return None

    scored.sort(key=lambda x: -x[0])
    best_score, best_name = scored[0]

    # Resolve ChEMBL ID(s) for the best-matching normalized name
    candidates_data = pref_index.get(best_name) or syn_index.get(best_name) or []
    chembl_id, is_ambiguous, candidate_ids = _resolve_candidates(candidates_data)

    return {
        "chembl_id": chembl_id,
        "match_method": "fuzzy",
        "match_score": round(best_score, 4),
        "candidate_count": len(candidate_ids),
        "is_ambiguous": is_ambiguous,
        "candidates": candidate_ids,
        "match_detail": json.dumps({"matched_chembl_name": best_name}),
    }


# ---------------------------------------------------------------------------
# Main mapping loop
# ---------------------------------------------------------------------------

def run_mapping(
    unique_names_df: pd.DataFrame,
    pref_index: Dict[str, List[dict]],
    syn_index: Dict[str, List[dict]],
    fuzzy: bool = True,
    fuzzy_min_score: float = 0.85,
) -> pd.DataFrame:
    """
    Apply the full tiered matching strategy to every unique normalized FAERS
    drug name.  Returns a DataFrame with one row per unique name and columns
    describing the match outcome.
    """
    token_index: Dict[str, Set[str]] = {}
    if fuzzy:
        print("Building token index for fuzzy prefiltering...")
        token_index = build_token_index(pref_index, syn_index)

    results: List[dict] = []
    n = len(unique_names_df)

    for i, row in enumerate(unique_names_df.itertuples(index=False), 1):
        if i % 5_000 == 0 or i == n:
            print(f"  Matching {i:,}/{n:,}...", end="\r")

        normalized: str = row.normalized
        raw_name: str = row.raw_name
        raw_ai: str = row.raw_ai
        field_used: str = row.field_used
        faers_count: int = int(getattr(row, "faers_count", 1))

        # Base result: defaults for unmatched
        result: dict = {
            "raw_name": raw_name,
            "raw_ai": raw_ai,
            "field_used": field_used,
            "normalized": normalized,
            "faers_count": faers_count,
            "chembl_id": None,
            "match_method": "unmatched",
            "match_score": 0.0,
            "candidate_count": 0,
            "is_ambiguous": False,
            "candidates": [],
            "match_detail": None,
        }

        if not normalized:
            result["match_method"] = "empty_normalized"
            results.append(result)
            continue

        # Tier 2a/2b: exact match on full normalized string
        match = match_single_name(normalized, pref_index, syn_index)

        # Tier 2c: if unmatched AND it's a combo, try component-wise match
        if match is None and is_combo(normalized):
            match = match_combo_name(normalized, pref_index, syn_index)

        # Tier 3: fuzzy fallback (only for still-unmatched names)
        if match is None and fuzzy:
            match = match_fuzzy(
                normalized, pref_index, syn_index, token_index,
                min_score=fuzzy_min_score,
            )

        if match:
            result.update(match)

        results.append(result)

    print()  # newline after progress indicator
    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Validation statistics
# ---------------------------------------------------------------------------

def compute_mapping_stats(results_df: pd.DataFrame) -> dict:
    """
    Compute and print the validation numbers your professor will ask about.
    Returns a dict that is also serialized to mapping_stats.json.
    """
    total = len(results_df)
    total_faers_rows = int(results_df["faers_count"].sum())

    method_counts: dict = results_df["match_method"].value_counts().to_dict()

    unmatched = results_df[results_df["match_method"] == "unmatched"]
    ambiguous = results_df[results_df["is_ambiguous"]]
    exact = results_df[results_df["match_method"].isin(
        ["exact_pref", "exact_synonym", "combo_all"]
    )]
    fuzzy_matched = results_df[results_df["match_method"] == "fuzzy"]
    matched = results_df[~results_df["match_method"].isin(
        ["unmatched", "empty_normalized"]
    )]

    matched_faers_rows = int(matched["faers_count"].sum())

    stats = {
        "total_unique_normalized_names":      total,
        "total_faers_drug_rows_represented":  total_faers_rows,
        "method_counts":                      method_counts,
        "n_matched":                          len(matched),
        "n_unmatched":                        len(unmatched),
        "n_ambiguous":                        len(ambiguous),
        "n_exact":                            len(exact),
        "n_fuzzy":                            len(fuzzy_matched),
        "exact_match_pct":   round(100 * len(exact)        / total, 2) if total else 0,
        "fuzzy_match_pct":   round(100 * len(fuzzy_matched) / total, 2) if total else 0,
        "unmatched_pct":     round(100 * len(unmatched)    / total, 2) if total else 0,
        "ambiguity_rate_pct":round(100 * len(ambiguous)    / total, 2) if total else 0,
        # Weighted by FAERS row frequency (more meaningful for coverage)
        "weighted_match_pct":round(100 * matched_faers_rows / total_faers_rows, 2)
                             if total_faers_rows else 0,
    }

    # Top-20 lists for inspection
    top_unmatched = (
        unmatched
        .nlargest(20, "faers_count")
        [["raw_name", "normalized", "faers_count"]]
        .to_dict("records")
    )
    top_ambiguous = (
        ambiguous
        .nlargest(20, "faers_count")
        [["raw_name", "normalized", "candidate_count", "faers_count"]]
        .to_dict("records")
    )
    stats["top_20_unmatched"] = top_unmatched
    stats["top_20_ambiguous"] = top_ambiguous

    # Print summary
    print("\n" + "=" * 60)
    print("FAERS -> ChEMBL Mapping Statistics")
    print("=" * 60)
    print(f"  Unique normalized names:      {total:,}")
    print(f"  FAERS DRUG rows represented:  {total_faers_rows:,}")
    print()
    print(f"  exact_pref    matches:        {method_counts.get('exact_pref', 0):,}")
    print(f"  exact_synonym matches:        {method_counts.get('exact_synonym', 0):,}")
    print(f"  combo_all     matches:        {method_counts.get('combo_all', 0):,}")
    print(f"  fuzzy         matches:        {method_counts.get('fuzzy', 0):,}")
    print(f"  empty_normalized:             {method_counts.get('empty_normalized', 0):,}")
    print(f"  unmatched:                    {method_counts.get('unmatched', 0):,}")
    print()
    print(f"  Exact match rate (by name):   {stats['exact_match_pct']:.1f}%")
    print(f"  Fuzzy match rate (by name):   {stats['fuzzy_match_pct']:.1f}%")
    print(f"  Unmatched rate (by name):     {stats['unmatched_pct']:.1f}%")
    print(f"  Ambiguity rate (by name):     {stats['ambiguity_rate_pct']:.1f}%")
    print(f"  Weighted match % (by FAERS):  {stats['weighted_match_pct']:.1f}%")

    print(f"\n  Top unmatched (by FAERS frequency):")
    for r in top_unmatched[:5]:
        print(f"    {r['raw_name']!r:40s}  norm={r['normalized']!r:30s}  n={r['faers_count']}")
    if len(top_unmatched) > 5:
        print(f"    ... ({len(top_unmatched) - 5} more in mapping_stats.json)")

    print(f"\n  Top ambiguous (by FAERS frequency):")
    for r in top_ambiguous[:5]:
        print(f"    {r['raw_name']!r:40s}  norm={r['normalized']!r:25s}  candidates={r['candidate_count']}")
    if len(top_ambiguous) > 5:
        print(f"    ... ({len(top_ambiguous) - 5} more in mapping_stats.json)")

    return stats


# ---------------------------------------------------------------------------
# Save artifacts
# ---------------------------------------------------------------------------

def save_artifacts(
    out_dir: str,
    results_df: pd.DataFrame,
    stats: dict,
) -> None:
    """Write all output files to out_dir."""
    os.makedirs(out_dir, exist_ok=True)

    # Full mapping table (parquet; requires pyarrow)
    mapping_path = os.path.join(out_dir, "faers_drug_to_chembl.parquet")
    results_df.to_parquet(mapping_path, index=False)
    print(f"\nSaved: {mapping_path}  ({len(results_df):,} rows)")

    # Unmatched names
    unmatched = results_df[results_df["match_method"] == "unmatched"].copy()
    unmatched_path = os.path.join(out_dir, "unmatched_drugs.csv")
    (
        unmatched[["raw_name", "raw_ai", "normalized", "faers_count"]]
        .sort_values("faers_count", ascending=False)
        .to_csv(unmatched_path, index=False)
    )
    print(f"Saved: {unmatched_path}  ({len(unmatched):,} rows)")

    # Ambiguous names
    ambiguous = results_df[results_df["is_ambiguous"]].copy()
    ambiguous_path = os.path.join(out_dir, "ambiguous_drugs.csv")
    (
        ambiguous[["raw_name", "raw_ai", "normalized", "candidate_count", "candidates", "faers_count"]]
        .sort_values("faers_count", ascending=False)
        .to_csv(ambiguous_path, index=False)
    )
    print(f"Saved: {ambiguous_path}  ({len(ambiguous):,} rows)")

    # Validation stats
    stats_path = os.path.join(out_dir, "mapping_stats.json")
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2, default=str)
    print(f"Saved: {stats_path}")


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Map FAERS drug names to ChEMBL molecule IDs (tiered strategy)."
    )
    parser.add_argument(
        "--no-fuzzy", action="store_true",
        help="Skip Tier 3 fuzzy matching (faster; strict exact-only output).",
    )
    parser.add_argument(
        "--fuzzy-threshold", type=float, default=0.85,
        help="Minimum SequenceMatcher.ratio() for a fuzzy match (default: 0.85).",
    )
    parser.add_argument(
        "--force-rebuild-index", action="store_true",
        help="Re-query ChEMBL and overwrite the cached chembl_name_index.parquet.",
    )
    parser.add_argument(
        "--max-rows", type=int, default=None,
        help="Limit FAERS DRUG rows for a quick test run.",
    )
    parser.add_argument(
        "--db-password", default=os.environ.get("PGPASSWORD", ""),
        help="PostgreSQL password (default: $PGPASSWORD env var).",
    )
    args = parser.parse_args()

    # ── 1. Connect to ChEMBL and build/load name index ────────────────
    print("=" * 60)
    print("Step 1: ChEMBL name index")
    print("=" * 60)
    conn = psycopg2.connect(
        dbname=_cfg.CHEMBL_DB_NAME,
        user=_cfg.CHEMBL_DB_USER,
        password=args.db_password,
        host=_cfg.CHEMBL_DB_HOST,
        port=_cfg.CHEMBL_DB_PORT,
    )
    print(f"Connected to ChEMBL ({_cfg.CHEMBL_DB_NAME}@{_cfg.CHEMBL_DB_HOST}:{_cfg.CHEMBL_DB_PORT})")

    try:
        pref_index, syn_index = build_chembl_name_index(
            conn, OUT_DIR, force_rebuild=args.force_rebuild_index,
        )
    finally:
        conn.close()
        print("Database connection closed")

    # ── 2. Load FAERS and build unique name table ─────────────────────
    print("\n" + "=" * 60)
    print("Step 2: FAERS drug names")
    print("=" * 60)
    faers_df = load_faers_drug_table(FAERS_ZIP, FAERS_DRUG_INNER)
    if args.max_rows:
        faers_df = faers_df.head(args.max_rows)
        print(f"(Test mode: limited to first {args.max_rows:,} rows)")

    unique_names_df = build_unique_name_table(faers_df)

    # ── 3. Run tiered matching ────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Step 3: Tiered matching")
    print("=" * 60)
    results_df = run_mapping(
        unique_names_df,
        pref_index,
        syn_index,
        fuzzy=not args.no_fuzzy,
        fuzzy_min_score=args.fuzzy_threshold,
    )

    # ── 4. Validation stats ───────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Step 4: Validation")
    print("=" * 60)
    stats = compute_mapping_stats(results_df)

    # ── 5. Save artifacts ─────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Step 5: Saving artifacts")
    print("=" * 60)
    save_artifacts(OUT_DIR, results_df, stats)

    print("\nDone. Artifacts written to:", OUT_DIR)
