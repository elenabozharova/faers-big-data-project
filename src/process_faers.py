"""
process_faers.py — Extract polypharmacy side-effect signals from FAERS.

Builds three Decagon inputs from raw FAERS adverse-event reports:
  combo  — (drug_i, drug_j, side_effect) co-occurrence triples
  mono   — (drug, side_effect) single-drug association pairs
  se_idx — MedDRA Preferred Term string → integer index

Pre-requisites (run in order before this script)
-------------------------------------------------
  1. map_faers_to_chembl.py  →  data/processed/faers/faers_drug_to_chembl.parquet
  2. build_decagon_graph.py  →  data/processed/decagon[_all]/drug_to_idx.json

Signal design
-------------
Mono:
  Reports where exactly ONE known drug (a drug in the Decagon universe) appears
  as PRIMARY SUSPECT (role_cod='PS').  This minimises confounding from
  co-administered drugs.

Combo:
  Reports with TWO OR MORE known drugs in any role (PS / SS / C / I).
  All (i, j) pairs with i < j are enumerated; side effects shared by the report
  are attributed to every pair.

Minimum-support filtering removes rare (drug, SE) combinations that would
introduce noise into model training.

Outputs  (data/processed/faers/)
---------------------------------
  side_effect_to_idx.json        MedDRA PT string -> integer index (sorted)
  combo_counts.parquet           (drug_i_idx int, drug_j_idx int, se_idx int,
                                   n_reports int)  — i < j always
  mono_counts.parquet            (drug_idx int, se_idx int, n_reports int)
  meta.json                      statistics and run configuration

Usage
-----
  python src/process_faers.py
  python src/process_faers.py --mode strict
  python src/process_faers.py --min-combo-reports 5 --min-mono-reports 3
  python src/process_faers.py --max-rows 200000   # test mode (first N DRUG rows)
"""

import argparse
import json
import os
import sys
import zipfile
from itertools import combinations
from typing import Dict, List, Tuple

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(_SRC_DIR)
DATA_RAW = os.path.join(PROJECT_ROOT, "data", "raw")
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")

# Ensure src/ is importable when run from project root
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

import config as _cfg  # noqa: E402
from normalize_drug_names import normalize_drug_string  # noqa: E402

FAERS_ZIP        = _cfg.FAERS_ZIP_PATH
FAERS_DRUG_INNER = _cfg.FAERS_DRUG_INNER
FAERS_REAC_INNER = _cfg.FAERS_REAC_INNER

OUT_DIR = os.path.join(DATA_PROCESSED, "faers")

# ---------------------------------------------------------------------------
# MedDRA PTs that are NOT adverse events
# These are administrative / medication-error / effectiveness terms that would
# pollute the side-effect signal.  Lowercase, as returned after normalization.
# ---------------------------------------------------------------------------
_SE_EXCLUDE: frozenset = frozenset({
    "drug ineffective",
    "off label use",
    "product use in unapproved indication",
    "intentional overdose",
    "accidental overdose",
    "drug interaction",
    "drug dose omission",
    "wrong drug administered",
    "drug dispensing error",
    "expired product administered",
    "inappropriate schedule of product administration",
    "product administered to patient of inappropriate age",
    "medication error",
    "no adverse event",
    "therapeutic response unexpected",
    "incorrect dose administered",
    "drug administration error",
    "product quality issue",
    "underdose",
    "overdose",
    "drug level increased",
    "drug level decreased",
    "drug level above therapeutic",
    "drug level below therapeutic",
    "drug resistance",
    "lack of efficacy",
})


# ---------------------------------------------------------------------------
# Step 1 — load pre-built mapping and decagon drug index
# ---------------------------------------------------------------------------

def load_drug_mapping(faers_dir: str) -> pd.DataFrame:
    """
    Load the FAERS drug name → ChEMBL ID mapping produced by map_faers_to_chembl.py.

    Columns used downstream: normalized, chembl_id, is_ambiguous, match_method.
    """
    path = os.path.join(faers_dir, "faers_drug_to_chembl.parquet")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Drug mapping not found at:\n  {path}\n"
            "Run map_faers_to_chembl.py first to produce this file."
        )
    df = pd.read_parquet(path)
    print(f"  Loaded drug mapping: {len(df):,} unique normalized names")
    return df


def load_decagon_drug_idx(mode: str) -> Dict[str, int]:
    """
    Load drug_to_idx.json for the chosen Decagon graph variant.

    mode='strict'  →  data/processed/decagon/drug_to_idx.json   (3,471 drugs)
    mode='all'     →  data/processed/decagon_all/drug_to_idx.json (672,337 drugs)
    """
    subdir = "decagon" if mode == "strict" else "decagon_all"
    path = os.path.join(DATA_PROCESSED, subdir, "drug_to_idx.json")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Decagon drug index not found at:\n  {path}\n"
            f"Run build_decagon_graph.py {'strict' if mode == 'strict' else 'all'} first."
        )
    with open(path, encoding="utf-8") as f:
        idx = json.load(f)
    print(f"  Decagon drug universe ({mode}): {len(idx):,} drugs")
    return idx


def build_name_to_decagon_lookup(
    mapping_df: pd.DataFrame,
    decagon_drug_idx: Dict[str, int],
    include_fuzzy: bool = False,
) -> Dict[str, int]:
    """
    Combine the FAERS→ChEMBL mapping table with the Decagon drug index to
    produce a direct lookup:  normalized_faers_name -> decagon_drug_idx.

    Rules applied:
      - Skip rows where chembl_id is None (unmatched, ambiguous, or combo_all).
      - Skip rows where is_ambiguous is True.
      - Skip rows whose chembl_id is not in the Decagon drug universe.
      - By default (include_fuzzy=False), skip rows where match_method='fuzzy'.
        Pass include_fuzzy=True (--include-fuzzy flag) to include them.
    """
    lookup: Dict[str, int] = {}
    n_no_id = 0
    n_not_in_decagon = 0
    n_fuzzy_skipped = 0

    for row in mapping_df.itertuples(index=False):
        cid = row.chembl_id
        if pd.isna(cid) or cid is None or row.is_ambiguous:
            n_no_id += 1
            continue
        if not include_fuzzy and row.match_method == "fuzzy":
            n_fuzzy_skipped += 1
            continue
        if cid not in decagon_drug_idx:
            n_not_in_decagon += 1
            continue
        lookup[row.normalized] = decagon_drug_idx[cid]

    print(f"  normalized->decagon entries: {len(lookup):,}")
    print(f"  Skipped (no/ambiguous ChEMBL ID): {n_no_id:,}")
    if not include_fuzzy:
        print(f"  Skipped (fuzzy match, strict mode): {n_fuzzy_skipped:,}")
    print(f"  Skipped (ChEMBL ID not in Decagon): {n_not_in_decagon:,}")
    return lookup


# ---------------------------------------------------------------------------
# Step 2 — load FAERS tables
# ---------------------------------------------------------------------------

def load_faers_table(zip_path: str, inner_path: str, max_rows: int = None) -> pd.DataFrame:
    """Load a FAERS ASCII table (dollar-delimited) from the quarterly ZIP."""
    with zipfile.ZipFile(zip_path) as z:
        with z.open(inner_path) as f:
            df = pd.read_csv(
                f, sep="$", dtype=str, low_memory=False,
                nrows=max_rows,
            )
    df.columns = df.columns.str.lower().str.strip()
    print(f"  {inner_path}: {len(df):,} rows, columns={list(df.columns)}")
    return df


# ---------------------------------------------------------------------------
# Step 3 — map FAERS DRUG rows to Decagon drug indices
# ---------------------------------------------------------------------------

def map_drug_rows_to_decagon(
    drug_df: pd.DataFrame,
    name_to_decagon: Dict[str, int],
) -> pd.DataFrame:
    """
    Apply Tier 0 field selection and Tier 1 normalization to every FAERS DRUG
    row, then look up the resulting normalized name in name_to_decagon.

    Tier 0: prefer prod_ai when it is present, non-empty (>2 chars), and not
            purely numeric; otherwise fall back to drugname.
    Tier 1: normalize_drug_string() from normalize_drug_names.py.

    Returns a DataFrame with columns:
      primaryid  (str)
      role_cod   (str: PS / SS / C / I / DN / '')
      drug_idx   (int or NaN for unmapped drugs)
    """
    prod_ai = (
        drug_df["prod_ai"].fillna("").astype(str).str.strip()
        if "prod_ai" in drug_df.columns
        else pd.Series([""] * len(drug_df), dtype=str)
    )
    drugname = (
        drug_df["drugname"].fillna("").astype(str).str.strip()
        if "drugname" in drug_df.columns
        else pd.Series([""] * len(drug_df), dtype=str)
    )

    # Tier 0: prefer prod_ai when meaningful
    is_meaningful = (prod_ai.str.len() > 2) & (
        ~prod_ai.str.replace(".", "", regex=False).str.isnumeric()
    )
    raw_name = prod_ai.where(is_meaningful, drugname)

    # Tier 1: normalize (row-by-row; ~30 s for 1.8 M rows)
    print("  Normalizing drug names (this may take ~30 s for the full table)...")
    normalized = raw_name.apply(normalize_drug_string)

    mapped_df = pd.DataFrame({
        "primaryid": drug_df["primaryid"].astype(str).str.strip(),
        "role_cod": (
            drug_df["role_cod"].fillna("").astype(str).str.strip()
            if "role_cod" in drug_df.columns
            else pd.Series([""] * len(drug_df), dtype=str)
        ),
        "drug_idx": normalized.map(name_to_decagon),
    })

    n_mapped = mapped_df["drug_idx"].notna().sum()
    print(
        f"  DRUG rows mapped to Decagon: {n_mapped:,} / {len(mapped_df):,} "
        f"({100 * n_mapped / len(mapped_df):.1f}%)"
    )
    return mapped_df


# ---------------------------------------------------------------------------
# Step 4 — filter REAC to genuine adverse events
# ---------------------------------------------------------------------------

def filter_reactions(
    reac_df: pd.DataFrame,
    se_exclude: frozenset = _SE_EXCLUDE,
) -> pd.DataFrame:
    """
    Lowercase MedDRA PT strings and remove non-adverse-event terms.

    Returns a DataFrame with columns:  primaryid (str), pt (str, lowercase).
    """
    pt_col = "pt" if "pt" in reac_df.columns else reac_df.columns[2]
    out = pd.DataFrame({
        "primaryid": reac_df["primaryid"].astype(str).str.strip(),
        "pt": reac_df[pt_col].fillna("").astype(str).str.strip().str.lower(),
    })
    before = len(out)
    out = out[(out["pt"] != "") & (~out["pt"].isin(se_exclude))]
    print(
        f"  REAC rows after filtering: {len(out):,} "
        f"(removed {before - len(out):,} noise rows)"
    )
    return out


# ---------------------------------------------------------------------------
# Step 5 — classify reports and aggregate signals
# ---------------------------------------------------------------------------

def extract_signals(
    mapped_drug_df: pd.DataFrame,
    filtered_reac_df: pd.DataFrame,
    min_combo_reports: int,
    min_mono_reports: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, int]]:
    """
    Build mono and combo side-effect count tables.

    Mono logic:
      Keep only primaryids where the primary-suspect (PS) drug set contains
      exactly ONE Decagon-known drug.  Map that drug with every reaction in
      the report.

    Combo logic:
      Keep primaryids where ANY role yields >= 2 Decagon-known drugs.
      Enumerate all ordered pairs (i < j) and attribute each reaction to
      every pair.

    Returns
    -------
    combo_df  : DataFrame (drug_i_idx, drug_j_idx, se_idx, n_reports), i < j
    mono_df   : DataFrame (drug_idx, se_idx, n_reports)
    se_to_idx : Dict[str, int] sorted MedDRA PT -> integer index
    """
    known = mapped_drug_df[mapped_drug_df["drug_idx"].notna()].copy()
    known["drug_idx"] = known["drug_idx"].astype(int)

    # ── Mono: exactly one PS drug known ────────────────────────────────────
    print("  Identifying mono reports (1 PS drug)...")
    ps = known[known["role_cod"] == "PS"]
    ps_per_report = ps.groupby("primaryid")["drug_idx"].nunique()
    mono_pids = ps_per_report[ps_per_report == 1].index

    mono_drug_map = (
        ps[ps["primaryid"].isin(mono_pids)]
        .groupby("primaryid")["drug_idx"]
        .first()          # single PS drug
        .reset_index()
    )
    print(f"    Mono reports: {len(mono_pids):,}")

    # ── Combo: >= 2 distinct known drugs (any role) ─────────────────────────
    print("  Identifying combo reports (>=2 known drugs)...")
    all_per_report = known.groupby("primaryid")["drug_idx"].nunique()
    combo_pids = all_per_report[all_per_report >= 2].index

    combo_drug_sets = (
        known[known["primaryid"].isin(combo_pids)]
        .groupby("primaryid")["drug_idx"]
        .apply(lambda x: sorted(set(x)))
        .reset_index()
        .rename(columns={"drug_idx": "drug_set"})
    )
    print(f"    Combo reports: {len(combo_pids):,}")

    # ── Join with REAC ───────────────────────────────────────────────────────
    print("  Joining with REAC table...")
    mono_joined = mono_drug_map.merge(filtered_reac_df, on="primaryid", how="inner")
    combo_joined = combo_drug_sets.merge(filtered_reac_df, on="primaryid", how="inner")
    print(f"    Mono join:  {len(mono_joined):,} rows")
    print(f"    Combo join: {len(combo_joined):,} rows")

    # ── Side effect index (union of both sets, sorted) ───────────────────────
    all_pts = set(mono_joined["pt"].unique()) | set(combo_joined["pt"].unique())
    se_to_idx: Dict[str, int] = {pt: i for i, pt in enumerate(sorted(all_pts))}
    print(f"  Unique side effects (MedDRA PTs): {len(se_to_idx):,}")

    mono_joined["se_idx"] = mono_joined["pt"].map(se_to_idx)
    combo_joined["se_idx"] = combo_joined["pt"].map(se_to_idx)

    # ── Mono counts ──────────────────────────────────────────────────────────
    print("  Aggregating mono counts...")
    mono_counts = (
        mono_joined
        .groupby(["drug_idx", "se_idx"])
        .size()
        .reset_index(name="n_reports")
    )
    mono_counts = mono_counts[mono_counts["n_reports"] >= min_mono_reports]
    mono_counts["drug_idx"] = mono_counts["drug_idx"].astype(int)
    mono_counts["se_idx"] = mono_counts["se_idx"].astype(int)
    print(
        f"  Mono (drug, SE) pairs >= {min_mono_reports} reports: "
        f"{len(mono_counts):,}"
    )

    # ── Combo pair expansion ─────────────────────────────────────────────────
    print("  Expanding combo drug pairs...")
    combo_rows: List[Tuple[int, int, int]] = []
    for row in combo_joined.itertuples(index=False):
        se_idx = row.se_idx
        drug_set: List[int] = row.drug_set
        for di, dj in combinations(drug_set, 2):
            combo_rows.append((min(di, dj), max(di, dj), se_idx))

    print(f"  Expanded {len(combo_rows):,} (pair, SE) observations")

    combo_pairs_df = pd.DataFrame(
        combo_rows, columns=["drug_i_idx", "drug_j_idx", "se_idx"]
    )
    combo_counts = (
        combo_pairs_df
        .groupby(["drug_i_idx", "drug_j_idx", "se_idx"])
        .size()
        .reset_index(name="n_reports")
    )
    combo_counts = combo_counts[combo_counts["n_reports"] >= min_combo_reports]
    for col in ["drug_i_idx", "drug_j_idx", "se_idx"]:
        combo_counts[col] = combo_counts[col].astype(int)
    print(
        f"  Combo (di, dj, SE) triples >= {min_combo_reports} reports: "
        f"{len(combo_counts):,}"
    )

    return combo_counts, mono_counts, se_to_idx


# ---------------------------------------------------------------------------
# Step 6 — save artifacts
# ---------------------------------------------------------------------------

def save_artifacts(
    out_dir: str,
    combo_df: pd.DataFrame,
    mono_df: pd.DataFrame,
    se_to_idx: Dict[str, int],
    meta: dict,
) -> None:
    os.makedirs(out_dir, exist_ok=True)

    se_path = os.path.join(out_dir, "side_effect_to_idx.json")
    with open(se_path, "w", encoding="utf-8") as f:
        json.dump(se_to_idx, f, indent=2)
    print(f"Saved: {se_path}  ({len(se_to_idx):,} side effects)")

    combo_path = os.path.join(out_dir, "combo_counts.parquet")
    combo_df.to_parquet(combo_path, index=False)
    print(f"Saved: {combo_path}  ({len(combo_df):,} rows)")

    mono_path = os.path.join(out_dir, "mono_counts.parquet")
    mono_df.to_parquet(mono_path, index=False)
    print(f"Saved: {mono_path}  ({len(mono_df):,} rows)")

    meta_path = os.path.join(out_dir, "meta.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    print(f"Saved: {meta_path}")


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Extract polypharmacy side-effect signals from FAERS for the "
            "Decagon model.  Requires map_faers_to_chembl.py and "
            "build_decagon_graph.py to have been run first."
        )
    )
    parser.add_argument(
        "--mode", choices=["strict", "all"], default="all",
        help=(
            "Drug universe to use.  'strict' = mechanism-based ChEMBL "
            "(3,471 drugs); 'all' = broad activity-based ChEMBL "
            "(672,337 drugs).  Default: all."
        ),
    )
    parser.add_argument(
        "--min-combo-reports", type=int, default=10,
        help=(
            "Minimum number of reports sharing a (drug_i, drug_j, SE) triple "
            "for it to be included in the combo output.  Default: 10."
        ),
    )
    parser.add_argument(
        "--min-mono-reports", type=int, default=5,
        help=(
            "Minimum number of reports for a (drug, SE) pair to be included "
            "in the mono output.  Default: 5."
        ),
    )
    parser.add_argument(
        "--include-fuzzy", action="store_true",
        help=(
            "Include fuzzy-matched drug names (match_method='fuzzy') from "
            "map_faers_to_chembl.py.  By default these are excluded so that "
            "only exact and synonym matches contribute to the signal."
        ),
    )
    parser.add_argument(
        "--max-rows", type=int, default=None,
        help=(
            "Limit DRUG and REAC tables to the first N rows each.  "
            "Intended for quick smoke-testing only."
        ),
    )
    args = parser.parse_args()

    print("=" * 60)
    print("FAERS Side-Effect Signal Extraction")
    print("=" * 60)
    print(f"  Mode:              {args.mode}")
    print(f"  Fuzzy matching:    {'included' if args.include_fuzzy else 'excluded (strict)'}")
    print(f"  Min combo reports: {args.min_combo_reports}")
    print(f"  Min mono reports:  {args.min_mono_reports}")
    if args.max_rows:
        print(f"  Max rows:          {args.max_rows:,}  (TEST MODE)")

    # ── 1. Load mappings ────────────────────────────────────────────────────
    print("\n--- Step 1: Load drug mapping and Decagon drug universe ---")
    mapping_df = load_drug_mapping(OUT_DIR)
    decagon_drug_idx = load_decagon_drug_idx(args.mode)
    name_to_decagon = build_name_to_decagon_lookup(
        mapping_df, decagon_drug_idx, include_fuzzy=args.include_fuzzy
    )

    # ── 2. Load FAERS tables ────────────────────────────────────────────────
    print("\n--- Step 2: Load FAERS DRUG and REAC tables ---")
    drug_df = load_faers_table(FAERS_ZIP, FAERS_DRUG_INNER, max_rows=args.max_rows)
    reac_df = load_faers_table(FAERS_ZIP, FAERS_REAC_INNER, max_rows=args.max_rows)

    # ── 3. Map drug rows to Decagon indices ─────────────────────────────────
    print("\n--- Step 3: Map FAERS drug names to Decagon drug indices ---")
    mapped_drug_df = map_drug_rows_to_decagon(drug_df, name_to_decagon)

    # ── 4. Filter reactions ─────────────────────────────────────────────────
    print("\n--- Step 4: Filter REAC to genuine adverse events ---")
    filtered_reac_df = filter_reactions(reac_df)

    # ── 5. Extract signals ───────────────────────────────────────────────────
    print("\n--- Step 5: Classify reports and aggregate signals ---")
    combo_df, mono_df, se_to_idx = extract_signals(
        mapped_drug_df,
        filtered_reac_df,
        min_combo_reports=args.min_combo_reports,
        min_mono_reports=args.min_mono_reports,
    )

    # ── 6. Save ─────────────────────────────────────────────────────────────
    print("\n--- Step 6: Save artifacts ---")
    meta = {
        "mode": args.mode,
        "include_fuzzy": args.include_fuzzy,
        "min_combo_reports": args.min_combo_reports,
        "min_mono_reports": args.min_mono_reports,
        "faers_zip": os.path.basename(FAERS_ZIP),
        "n_decagon_drugs_total": len(decagon_drug_idx),
        "n_decagon_drugs_with_faers_name": len(name_to_decagon),
        "n_side_effects": len(se_to_idx),
        "n_combo_triples": len(combo_df),
        "n_combo_drug_pairs": (
            int(combo_df[["drug_i_idx", "drug_j_idx"]].drop_duplicates().__len__())
            if len(combo_df) > 0 else 0
        ),
        "n_mono_pairs": len(mono_df),
        "n_mono_drugs_covered": (
            int(mono_df["drug_idx"].nunique()) if len(mono_df) > 0 else 0
        ),
    }
    save_artifacts(OUT_DIR, combo_df, mono_df, se_to_idx, meta)

    print("\n" + "=" * 60)
    print("FAERS processing complete!")
    print("=" * 60)
    print(f"  Side effects indexed:  {len(se_to_idx):,}")
    print(f"  Combo (di,dj,SE):      {len(combo_df):,}")
    print(f"  Mono  (drug,SE):       {len(mono_df):,}")
    print(f"  Artifacts in:          {OUT_DIR}")
