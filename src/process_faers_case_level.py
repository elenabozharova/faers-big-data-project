"""
process_faers_case_level.py -- Corrected CASE-LEVEL FAERS preprocessing (NEW pipeline).

This is a new, independent entry point. It does not modify, import the
internals of, or change the behavior of process_faers.py. That script keeps
running exactly as before.

Why a new pipeline
-------------------
process_faers.py treats every `primaryid` as an independent case. FAERS
resubmits cases over time: the DEMO table's `caseid` + `caseversion` identify
the same real-world case reported/amended multiple times under different
`primaryid`s. Counting every primaryid inflates report counts and distorts
co-occurrence signal strength. This pipeline resolves CASEID/CASEVERSION so
each real-world case is represented exactly once before any drug-pair or
side-effect signal is derived.

Milestone status
-----------------
  [x] Milestone 1 -- CLI skeleton, path wiring, placeholders
  [x] Milestone 2 -- load raw DEMO/DRUG/REAC; raw-table diagnostics           (this file)
  [ ] Milestone 3 -- resolve CASEID/CASEVERSION; filter to retained primaryids; dedupe reactions/drugs per case
  [ ] Milestone 4 -- canonical (drug_i, drug_j) pair generation per case
  [ ] Milestone 5 -- drug_pair_reaction_counts aggregation (PRR comes later, not in Milestone 5 either)
  [ ] Milestone 6 -- validation checks (distinct CASEID counts, no duplicate caseids, etc.)

Milestone 2 scope (this file)
------------------------------
Loads DEMO, DRUG, and REAC as-is and reports raw-table diagnostics only:
  - no CASEID/CASEVERSION resolution
  - no drug-name mapping/normalization
  - no reaction deduplication
  - no drug-pair generation
`caseversion` is coerced to a nullable numeric column; rows that fail
coercion are counted and flagged, not dropped, so the raw row count is
never silently changed by this step.

This file does NOT:
  - implement full preprocessing
  - generate drug pairs
  - compute PRR
  - touch R-GCN, ChEMBL, STRING, or Decagon adjacency / training code

Usage
-----
  python src/process_faers_case_level.py
  python src/process_faers_case_level.py --quarter 2025Q4 --max-rows 50000
  python src/process_faers_case_level.py --faers-dir data/raw/faers/extracted_2025Q4
"""

import argparse
import json
import os
import re
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(_SRC_DIR)
DATA_RAW = os.path.join(PROJECT_ROOT, "data", "raw")
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")

if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

import config as _cfg  # noqa: E402

# Reused as-is (read-only helpers; not modified, not copied):
from process_faers import load_faers_table  # noqa: E402  -- used in Milestone 2
from normalize_drug_names import normalize_drug_string  # noqa: E402  -- used in Milestone 3/4

DEFAULT_OUT_DIR = os.path.join(DATA_PROCESSED, "faers_case_level")
DEFAULT_MAPPING_FILE = os.path.join(DATA_PROCESSED, "faers", "faers_drug_to_chembl.parquet")


def _infer_quarter_from_zip(zip_path: str) -> str:
    """Best-effort extraction of a '2025Q4'-style label from a FAERS zip filename."""
    match = re.search(r"(\d{4}Q[1-4])", os.path.basename(zip_path), flags=re.IGNORECASE)
    return match.group(1).upper() if match else "unknown"


# ---------------------------------------------------------------------------
# Milestone 2 -- load raw DEMO / DRUG / REAC and report raw-table diagnostics
# ---------------------------------------------------------------------------

def _load_raw_faers_table(
    faers_zip: str,
    faers_dir: str,
    inner_path: str,
    max_rows: int = None,
) -> pd.DataFrame:
    """
    Load a dollar-delimited FAERS ASCII table from either an extracted
    directory (--faers-dir) or a quarterly zip (--faers-zip).

    Zip mode reuses process_faers.load_faers_table() unmodified. Directory
    mode is new (the old pipeline never reads from an extracted directory).
    """
    if faers_dir:
        candidates = [
            os.path.join(faers_dir, inner_path),
            os.path.join(faers_dir, os.path.basename(inner_path)),
        ]
        path = next((p for p in candidates if os.path.exists(p)), None)
        if path is None:
            raise FileNotFoundError(
                f"Could not find {inner_path!r} under {faers_dir!r}. Tried: {candidates}"
            )
        df = pd.read_csv(path, sep="$", dtype=str, low_memory=False, nrows=max_rows)
        df.columns = df.columns.str.lower().str.strip()
        print(f"  {inner_path} (dir): {len(df):,} rows, columns={list(df.columns)}")
        return df
    return load_faers_table(faers_zip, inner_path, max_rows=max_rows)


def _require_columns(df: pd.DataFrame, required: list, table_name: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{table_name} table is missing required column(s): {missing}. "
            f"Available columns: {list(df.columns)}"
        )


def _warn_missing_columns(df: pd.DataFrame, expected: list, table_name: str) -> None:
    missing = [c for c in expected if c not in df.columns]
    if missing:
        print(f"  [WARN] {table_name}: expected column(s) not found: {missing}")


def _is_missing(series: pd.Series) -> pd.Series:
    """True where a raw string-typed FAERS field is null or blank."""
    as_str = series.astype(str)
    return series.isna() | (as_str.str.strip() == "") | (as_str.str.strip().str.lower() == "nan")


def load_demo_table(
    faers_zip: str,
    faers_dir: str,
    inner_path: str,
    max_rows: int = None,
) -> pd.DataFrame:
    """
    Load the raw FAERS DEMO table.

    Guarantees primaryid, caseid, caseversion are present (raises otherwise).
    Warns (does not fail) if fda_dt is absent. Keeps primaryid/caseid as
    strings. Adds caseversion_numeric (nullable Int64, safely coerced) and
    caseversion_malformed (bool) without dropping or fixing bad rows --
    CASEID/CASEVERSION resolution itself is Milestone 3.
    """
    print("  Loading DEMO table...")
    df = _load_raw_faers_table(faers_zip, faers_dir, inner_path, max_rows=max_rows)
    _require_columns(df, ["primaryid", "caseid", "caseversion"], "DEMO")
    if "fda_dt" not in df.columns:
        print("  [WARN] DEMO: 'fda_dt' column not found (optional); continuing without it")

    df["primaryid"] = df["primaryid"].astype(str).str.strip()
    df["caseid"] = df["caseid"].astype(str).str.strip()

    caseversion_str = df["caseversion"].fillna("").astype(str).str.strip()
    is_missing = caseversion_str == ""
    numeric = pd.to_numeric(caseversion_str, errors="coerce")
    is_malformed = (~is_missing) & numeric.isna()

    df["caseversion_numeric"] = numeric.astype("Int64")
    df["caseversion_malformed"] = is_malformed

    n_malformed = int(is_malformed.sum())
    if n_malformed:
        print(f"  [WARN] DEMO: {n_malformed:,} rows have non-numeric caseversion values (kept, flagged)")

    print(f"  DEMO loaded: {len(df):,} rows")
    return df


def load_drug_table(
    faers_zip: str,
    faers_dir: str,
    inner_path: str,
    max_rows: int = None,
) -> pd.DataFrame:
    """Load the raw FAERS DRUG table. Guarantees primaryid; warns on missing drugname/prod_ai/role_cod."""
    print("  Loading DRUG table...")
    df = _load_raw_faers_table(faers_zip, faers_dir, inner_path, max_rows=max_rows)
    _require_columns(df, ["primaryid"], "DRUG")
    _warn_missing_columns(df, ["drugname", "prod_ai", "role_cod"], "DRUG")
    df["primaryid"] = df["primaryid"].astype(str).str.strip()
    print(f"  DRUG loaded: {len(df):,} rows")
    return df


def load_reac_table(
    faers_zip: str,
    faers_dir: str,
    inner_path: str,
    max_rows: int = None,
) -> pd.DataFrame:
    """Load the raw FAERS REAC table. Guarantees primaryid; falls back positionally to find 'pt' if unnamed."""
    print("  Loading REAC table...")
    df = _load_raw_faers_table(faers_zip, faers_dir, inner_path, max_rows=max_rows)
    _require_columns(df, ["primaryid"], "REAC")
    if "pt" not in df.columns:
        fallback_col = df.columns[2] if len(df.columns) > 2 else None
        if fallback_col:
            print(f"  [WARN] REAC: 'pt' column not found; using '{fallback_col}' as pt (positional fallback)")
            df["pt"] = df[fallback_col]
        else:
            print("  [WARN] REAC: 'pt' column not found and no fallback column available")
    df["primaryid"] = df["primaryid"].astype(str).str.strip()
    print(f"  REAC loaded: {len(df):,} rows")
    return df


def compute_raw_diagnostics(
    demo_df: pd.DataFrame,
    drug_df: pd.DataFrame,
    reac_df: pd.DataFrame,
) -> dict:
    """Compute the raw-table diagnostics required by Milestone 2 (counts, uniques, top values)."""
    demo_primaryid = demo_df["primaryid"]
    demo_caseid = demo_df["caseid"]
    caseversion_numeric = demo_df["caseversion_numeric"]

    demo_diag = {
        "n_rows": int(len(demo_df)),
        "n_unique_primaryid": int(demo_primaryid.nunique()),
        "n_unique_caseid": int(demo_caseid.nunique()),
        "n_missing_primaryid": int(_is_missing(demo_primaryid).sum()),
        "n_missing_caseid": int(_is_missing(demo_caseid).sum()),
        "n_missing_caseversion": int((demo_df["caseversion"].fillna("").astype(str).str.strip() == "").sum()),
        "n_malformed_caseversion": int(demo_df["caseversion_malformed"].sum()),
        "caseversion_min": int(caseversion_numeric.min()) if caseversion_numeric.notna().any() else None,
        "caseversion_max": int(caseversion_numeric.max()) if caseversion_numeric.notna().any() else None,
    }

    drug_primaryid = drug_df["primaryid"]
    top_role_cod = (
        drug_df["role_cod"].fillna("(missing)").astype(str).str.strip().replace("", "(missing)")
        .value_counts().head(10)
        if "role_cod" in drug_df.columns else pd.Series(dtype=int)
    )
    drug_diag = {
        "n_rows": int(len(drug_df)),
        "n_unique_primaryid": int(drug_primaryid.nunique()),
        "n_missing_primaryid": int(_is_missing(drug_primaryid).sum()),
        "top_role_cod": {str(k): int(v) for k, v in top_role_cod.items()},
    }

    reac_primaryid = reac_df["primaryid"]
    top_pt = (
        reac_df["pt"].fillna("(missing)").astype(str).str.strip().replace("", "(missing)")
        .value_counts().head(10)
        if "pt" in reac_df.columns else pd.Series(dtype=int)
    )
    reac_diag = {
        "n_rows": int(len(reac_df)),
        "n_unique_primaryid": int(reac_primaryid.nunique()),
        "n_missing_primaryid": int(_is_missing(reac_primaryid).sum()),
        "top_pt": {str(k): int(v) for k, v in top_pt.items()},
    }

    return {"demo": demo_diag, "drug": drug_diag, "reac": reac_diag}


def print_raw_diagnostics(diagnostics: dict) -> None:
    d = diagnostics["demo"]
    print("\n  DEMO diagnostics:")
    print(f"    rows:                  {d['n_rows']:,}")
    print(f"    unique primaryid:      {d['n_unique_primaryid']:,}")
    print(f"    unique caseid:         {d['n_unique_caseid']:,}")
    print(f"    missing primaryid:     {d['n_missing_primaryid']:,}")
    print(f"    missing caseid:        {d['n_missing_caseid']:,}")
    print(f"    missing caseversion:   {d['n_missing_caseversion']:,}")
    print(f"    malformed caseversion: {d['n_malformed_caseversion']:,}")
    print(f"    caseversion min/max:   {d['caseversion_min']} / {d['caseversion_max']}")

    g = diagnostics["drug"]
    print("\n  DRUG diagnostics:")
    print(f"    rows:                  {g['n_rows']:,}")
    print(f"    unique primaryid:      {g['n_unique_primaryid']:,}")
    print(f"    missing primaryid:     {g['n_missing_primaryid']:,}")
    print("    top role_cod values:")
    for k, v in g["top_role_cod"].items():
        print(f"      {k:<20} {v:,}")

    r = diagnostics["reac"]
    print("\n  REAC diagnostics:")
    print(f"    rows:                  {r['n_rows']:,}")
    print(f"    unique primaryid:      {r['n_unique_primaryid']:,}")
    print(f"    missing primaryid:     {r['n_missing_primaryid']:,}")
    print("    top PT values:")
    for k, v in r["top_pt"].items():
        print(f"      {k:<30} {v:,}")


def save_raw_tables(
    out_dir: str,
    demo_df: pd.DataFrame,
    drug_df: pd.DataFrame,
    reac_df: pd.DataFrame,
) -> None:
    demo_path = os.path.join(out_dir, "raw_demo.parquet")
    demo_df.to_parquet(demo_path, index=False)
    print(f"  Saved: {demo_path}  ({len(demo_df):,} rows)")

    drug_path = os.path.join(out_dir, "raw_drug.parquet")
    drug_df.to_parquet(drug_path, index=False)
    print(f"  Saved: {drug_path}  ({len(drug_df):,} rows)")

    reac_path = os.path.join(out_dir, "raw_reac.parquet")
    reac_df.to_parquet(reac_path, index=False)
    print(f"  Saved: {reac_path}  ({len(reac_df):,} rows)")


def save_raw_diagnostics(out_dir: str, diagnostics: dict) -> None:
    path = os.path.join(out_dir, "raw_table_diagnostics.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    print(f"  Saved: {path}")


def resolve_case_versions(demo_df: pd.DataFrame) -> pd.DataFrame:
    """
    TODO (Milestone 3): For each caseid, keep only the primaryid belonging to
    the highest caseversion (the most recent amendment/resubmission),
    producing the authoritative primaryid -> caseid mapping ("retained
    primaryids") used to filter DRUG and REAC downstream.
    """
    print("  [TODO] resolve_case_versions() not yet implemented (Milestone 3)")
    return pd.DataFrame(columns=["primaryid", "caseid"])


# ---------------------------------------------------------------------------
# Milestone 3 (TODO) -- filter to retained PRIMARYIDs; dedupe reactions/drugs
# ---------------------------------------------------------------------------

def filter_to_retained_primaryids(df: pd.DataFrame, retained_primaryids) -> pd.DataFrame:
    """TODO (Milestone 3): filter a DRUG or REAC dataframe down to retained primaryids only."""
    print("  [TODO] filter_to_retained_primaryids() not yet implemented (Milestone 3)")
    return df


def build_case_reactions(reac_df: pd.DataFrame, case_map_df: pd.DataFrame) -> pd.DataFrame:
    """
    TODO (Milestone 3): Join filtered REAC with the primaryid->caseid map and
    drop duplicate (caseid, pt) rows so each case reports a given reaction once.
    Output columns: caseid, pt.  Saved later as case_reactions.parquet.
    """
    print("  [TODO] build_case_reactions() not yet implemented (Milestone 3)")
    return pd.DataFrame(columns=["caseid", "pt"])


def build_case_drugs(drug_df: pd.DataFrame, case_map_df: pd.DataFrame) -> pd.DataFrame:
    """
    TODO (Milestone 3): Join filtered DRUG with the primaryid->caseid map,
    normalize drug names (normalize_drug_string), and drop duplicate
    (caseid, normalized_drug) rows, merging role_cod where needed.
    Output columns: caseid, drug, role_cod.  Saved later as case_drugs.parquet.
    """
    print("  [TODO] build_case_drugs() not yet implemented (Milestone 3)")
    return pd.DataFrame(columns=["caseid", "drug", "role_cod"])


# ---------------------------------------------------------------------------
# Milestone 4 (TODO) -- canonical drug-pair generation per case
# ---------------------------------------------------------------------------

def build_case_drug_pairs(case_drugs_df: pd.DataFrame) -> pd.DataFrame:
    """
    TODO (Milestone 4): Per case, enumerate canonical (drug_i, drug_j) pairs
    with drug_i < drug_j (no self-pairs, no direction duplicates).
    Output columns: caseid, drug_i, drug_j.  Saved later as case_drug_pairs.parquet.
    """
    print("  [TODO] build_case_drug_pairs() not yet implemented (Milestone 4)")
    return pd.DataFrame(columns=["caseid", "drug_i", "drug_j"])


# ---------------------------------------------------------------------------
# Milestone 5 (TODO) -- drug_pair_reaction_counts aggregation (no PRR yet)
# ---------------------------------------------------------------------------

def build_drug_pair_reaction_counts(
    case_drug_pairs_df: pd.DataFrame,
    case_reactions_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    TODO (Milestone 5): Join case_drug_pairs with case_reactions on caseid and
    count distinct cases per (drug_i, drug_j, pt). This is a count table only;
    PRR / disproportionality statistics are explicitly out of scope here.
    Output columns: drug_i, drug_j, pt, n_cases.
    Saved later as drug_pair_reaction_counts.parquet.
    """
    print("  [TODO] build_drug_pair_reaction_counts() not yet implemented (Milestone 5)")
    return pd.DataFrame(columns=["drug_i", "drug_j", "pt", "n_cases"])


# ---------------------------------------------------------------------------
# Milestone 6 (TODO) -- validation checks
# ---------------------------------------------------------------------------

def run_validation_checks(
    demo_df: pd.DataFrame,
    case_map_df: pd.DataFrame,
    case_reactions_df: pd.DataFrame,
    case_drugs_df: pd.DataFrame,
) -> None:
    """
    TODO (Milestone 6): Sanity checks, e.g.:
      - distinct caseid count in case_map_df matches distinct caseid in DEMO
      - no primaryid appears under two different caseids
      - no case has zero reactions or zero drugs after filtering
      - case_drug_pairs has no (drug_i == drug_j) rows and no reversed duplicates
    """
    print("  [TODO] run_validation_checks() not yet implemented (Milestone 6)")


# ---------------------------------------------------------------------------
# Save artifacts (TODO -- wired for Milestone 2+, no-op for now)
# ---------------------------------------------------------------------------

def save_case_level_artifacts(
    out_dir: str,
    case_reactions_df: pd.DataFrame,
    case_drugs_df: pd.DataFrame,
    case_drug_pairs_df: pd.DataFrame,
    drug_pair_reaction_counts_df: pd.DataFrame,
) -> None:
    """
    TODO: persist case_reactions.parquet, case_drugs.parquet,
    case_drug_pairs.parquet, and drug_pair_reaction_counts.parquet to out_dir.
    """
    print(f"  [TODO] save_case_level_artifacts() not yet implemented -- target dir: {out_dir}")


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Corrected case-level FAERS preprocessing pipeline (NEW). "
            "Resolves CASEID/CASEVERSION so each real-world case is counted "
            "once, ahead of drug-only R-GCN construction. Milestone 2: loads "
            "raw DEMO/DRUG/REAC and reports raw-table diagnostics only -- "
            "see file docstring for milestone status."
        )
    )
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "--faers-zip", default=_cfg.FAERS_ZIP_PATH,
        help="Path to the quarterly FAERS ASCII zip. Default: config.FAERS_ZIP_PATH.",
    )
    input_group.add_argument(
        "--faers-dir", default=None,
        help=(
            "Path to an already-extracted FAERS ASCII directory, used instead "
            "of --faers-zip. Mutually exclusive with --faers-zip."
        ),
    )
    parser.add_argument(
        "--quarter", default=None,
        help=(
            "Quarter label, e.g. 2025Q4. Used for output labeling and to "
            "derive default inner file names (DEMO/DRUG/REAC{qq}.txt). "
            "If omitted, inferred from the zip filename when possible."
        ),
    )
    parser.add_argument(
        "--out-dir", default=DEFAULT_OUT_DIR,
        help=f"Output directory for case-level artifacts. Default: {DEFAULT_OUT_DIR}",
    )
    parser.add_argument(
        "--mapping-file", default=DEFAULT_MAPPING_FILE,
        help=(
            "Path to the FAERS drug name -> ChEMBL ID mapping parquet "
            f"produced by map_faers_to_chembl.py. Default: {DEFAULT_MAPPING_FILE}"
        ),
    )
    parser.add_argument(
        "--include-fuzzy", action="store_true",
        help=(
            "Include fuzzy-matched drug names from the mapping file. Reserved "
            "for the drug-normalization step in a later milestone; unused in "
            "Milestone 1."
        ),
    )
    parser.add_argument(
        "--max-rows", type=int, default=None,
        help="Limit DEMO/DRUG/REAC tables to the first N rows each. Smoke-testing only.",
    )
    args = parser.parse_args()

    quarter = args.quarter or _infer_quarter_from_zip(args.faers_zip)

    print("=" * 60)
    print("FAERS CASE-LEVEL Preprocessing Pipeline (NEW)")
    print("=" * 60)
    print("This is the corrected, case-level pipeline. It resolves")
    print("CASEID/CASEVERSION instead of treating every PRIMARYID as an")
    print("independent case. It does not modify process_faers.py.")
    print()
    print(f"  FAERS input:     {args.faers_dir or args.faers_zip}")
    print(f"  Input mode:      {'directory' if args.faers_dir else 'zip'}")
    print(f"  Quarter:         {quarter}")
    print(f"  Mapping file:    {args.mapping_file}")
    print(f"  Include fuzzy:   {args.include_fuzzy}")
    print(f"  Output dir:      {args.out_dir}")
    if args.max_rows:
        print(f"  Max rows:        {args.max_rows:,}  (TEST MODE)")

    os.makedirs(args.out_dir, exist_ok=True)

    print("\n--- Milestone 2: load raw DEMO/DRUG/REAC ---")
    demo_df = load_demo_table(args.faers_zip, args.faers_dir, _cfg.FAERS_DEMO_INNER, max_rows=args.max_rows)
    drug_df = load_drug_table(args.faers_zip, args.faers_dir, _cfg.FAERS_DRUG_INNER, max_rows=args.max_rows)
    reac_df = load_reac_table(args.faers_zip, args.faers_dir, _cfg.FAERS_REAC_INNER, max_rows=args.max_rows)

    print("\n--- Milestone 2: raw-table diagnostics ---")
    diagnostics = compute_raw_diagnostics(demo_df, drug_df, reac_df)
    print_raw_diagnostics(diagnostics)

    print("\n--- Milestone 2: saving raw tables + diagnostics ---")
    save_raw_tables(args.out_dir, demo_df, drug_df, reac_df)
    save_raw_diagnostics(args.out_dir, diagnostics)

    print("\n--- Milestone 3 (pending): resolve CASEID/CASEVERSION; filter to retained PRIMARYIDs; dedupe reactions/drugs ---")
    case_map_df = resolve_case_versions(demo_df)
    case_reactions_df = build_case_reactions(reac_df, case_map_df)
    case_drugs_df = build_case_drugs(drug_df, case_map_df)

    print("\n--- Milestone 4 (pending): canonical drug-pair generation ---")
    case_drug_pairs_df = build_case_drug_pairs(case_drugs_df)

    print("\n--- Milestone 5 (pending): drug_pair_reaction_counts aggregation ---")
    drug_pair_reaction_counts_df = build_drug_pair_reaction_counts(case_drug_pairs_df, case_reactions_df)

    print("\n--- Milestone 6 (pending): validation checks ---")
    run_validation_checks(demo_df, case_map_df, case_reactions_df, case_drugs_df)

    print("\n--- Save (pending): case-level artifacts ---")
    save_case_level_artifacts(
        args.out_dir, case_reactions_df, case_drugs_df,
        case_drug_pairs_df, drug_pair_reaction_counts_df,
    )

    print("\n" + "=" * 60)
    print("Milestone 2 complete: raw DEMO/DRUG/REAC loaded and diagnosed.")
    print("No CASEID/CASEVERSION resolution, drug mapping, deduplication, or")
    print("drug-pair generation was performed yet.")
    print(f"Raw artifacts in: {args.out_dir}")
    print("Next: implement Milestone 3 (resolve CASEID/CASEVERSION; dedupe).")
    print("=" * 60)


if __name__ == "__main__":
    main()
