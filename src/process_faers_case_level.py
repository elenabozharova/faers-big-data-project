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
  [x] Milestone 2 -- load raw DEMO/DRUG/REAC; raw-table diagnostics
  [x] Milestone 3 -- resolve CASEID/CASEVERSION; filter to retained primaryids
  [x] Milestone 4 -- dedupe reactions per case (case_reactions.parquet)          (this file)
  [ ] Milestone 4b -- dedupe drugs per case (case_drugs.parquet)
  [ ] Milestone 5 -- canonical (drug_i, drug_j) pair generation per case
  [ ] Milestone 6 -- drug_pair_reaction_counts aggregation (no PRR)
  [ ] Milestone 7 -- global validation checks / PRR

Milestone 2 scope
------------------
Loads DEMO, DRUG, and REAC as-is and reports raw-table diagnostics only:
  - no CASEID/CASEVERSION resolution
  - no drug-name mapping/normalization
  - no reaction deduplication
  - no drug-pair generation
`caseversion` is coerced to a nullable numeric column; rows that fail
coercion are counted and flagged, not dropped, so the raw row count is
never silently changed by this step.

Milestone 3 scope
------------------
Resolves one retained PRIMARYID per CASEID (highest CASEVERSION, tie-broken
by latest fda_dt then largest PRIMARYID) and filters raw DRUG/REAC down to
rows whose PRIMARYID was retained. Still explicitly out of scope:
  - reaction/drug deduplication within a case (Milestone 4)
  - drug-name mapping/normalization (Milestone 4b+)
  - drug-pair generation (Milestone 5)
  - PRR / disproportionality statistics (Milestone 7)

Milestone 4 scope (this file)
------------------------------
Builds case_reactions.parquet: one row per (caseid, reaction_pt). Reuses
process_faers._SE_EXCLUDE (administrative/non-clinical PT strings) read-only,
purely for PT filtering -- no old row/primaryid-counting logic is reused.
Still explicitly out of scope:
  - drug deduplication within a case (Milestone 4b)
  - drug-name mapping/normalization (Milestone 4b+)
  - drug-pair generation (Milestone 5)
  - PRR / disproportionality statistics (Milestone 7)

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

# Reused as-is (read-only helpers/constants; not modified, not copied):
from process_faers import load_faers_table  # noqa: E402  -- used in Milestone 2
from normalize_drug_names import normalize_drug_string  # noqa: E402  -- used in Milestone 4b+

try:
    # Administrative/non-clinical MedDRA PT strings (plain frozenset of
    # data, no row/primaryid-counting logic). Used only to filter reaction
    # PTs in build_case_reactions(); process_faers.py itself is untouched.
    from process_faers import _SE_EXCLUDE  # noqa: E402
except ImportError:
    _SE_EXCLUDE = None  # exclusion list unavailable -- Milestone 4 will skip this filter and say so

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
    Resolve one retained primaryid per caseid from CASEID/CASEVERSION.

    Priority order (each column breaks ties left over from the previous one):
      1. highest caseversion_numeric.
         Missing/malformed caseversion is ranked -1 (lowest priority) via
         fillna(-1), so a malformed-caseversion row is only ever retained
         when it is the *only* row available for that caseid.
      2. latest fda_dt, if the column exists. FAERS dates are usually
         YYYYMMDD; partial dates and unparseable values are coerced to NaT
         and ranked lowest (same "least priority, never crash" rule).
      3. largest primaryid (numeric).

    Rows with a missing caseid cannot be assigned to any case and are
    excluded entirely -- counted here (and again in the diagnostics) rather
    than silently dropped.
    """
    working = demo_df.copy()

    missing_caseid_mask = _is_missing(working["caseid"])
    n_missing_caseid = int(missing_caseid_mask.sum())
    if n_missing_caseid:
        print(f"  [WARN] DEMO: {n_missing_caseid:,} rows have missing caseid; excluded from case resolution")
    working = working[~missing_caseid_mask].copy()

    working["_caseversion_rank"] = working["caseversion_numeric"].fillna(-1)

    if "fda_dt" in working.columns:
        fda_dt_parsed = pd.to_datetime(working["fda_dt"], format="%Y%m%d", errors="coerce")
        working["_fda_dt_rank"] = fda_dt_parsed.fillna(pd.Timestamp.min)
    else:
        working["_fda_dt_rank"] = 0  # no fda_dt available: tie-break falls through to primaryid

    working["_primaryid_rank"] = pd.to_numeric(working["primaryid"], errors="coerce").fillna(-1)

    working = working.sort_values(
        by=["caseid", "_caseversion_rank", "_fda_dt_rank", "_primaryid_rank"],
        ascending=[True, False, False, False],
        kind="mergesort",  # stable -> deterministic winner on exact ties
    )

    n_versions = working.groupby("caseid").size().rename("n_versions_for_case")

    retained = working.drop_duplicates(subset="caseid", keep="first").copy()
    retained = retained.merge(n_versions, on="caseid", how="left")

    out_cols = {
        "caseid": retained["caseid"],
        "retained_primaryid": retained["primaryid"],
        "caseversion": retained["caseversion_numeric"],
        "n_versions_for_case": retained["n_versions_for_case"],
    }
    if "fda_dt" in retained.columns:
        out_cols["fda_dt"] = retained["fda_dt"]

    retained_cases_df = pd.DataFrame(out_cols).reset_index(drop=True)

    print(
        f"  Resolved {len(retained_cases_df):,} retained cases from "
        f"{len(working):,} usable DEMO rows ({n_missing_caseid:,} rows excluded for missing caseid)"
    )
    return retained_cases_df


# ---------------------------------------------------------------------------
# Milestone 3 -- filter to retained PRIMARYIDs
# ---------------------------------------------------------------------------

def filter_to_retained_primaryids(
    df: pd.DataFrame,
    retained_cases_df: pd.DataFrame,
    table_name: str,
) -> pd.DataFrame:
    """
    Filter a raw DRUG or REAC dataframe down to rows whose primaryid is a
    retained_primaryid, attaching the *resolved* caseid via merge.

    DRUG/REAC already carry their own raw caseid column, but it is dropped
    before merging so the resolved caseid from retained_cases_df (the
    authoritative one) doesn't collide with it as caseid_x/caseid_y.
    """
    retained_map = (
        retained_cases_df[["retained_primaryid", "caseid"]]
        .drop_duplicates()
        .rename(columns={"retained_primaryid": "primaryid", "caseid": "caseid_resolved"})
    )
    before = len(df)
    filtered = df.merge(retained_map, on="primaryid", how="inner")

    # DRUG/REAC carry their own raw caseid. It should always agree with the
    # DEMO-derived (authoritative) one -- log if it doesn't, since a
    # mismatch would indicate a data-quality problem worth investigating,
    # not something to silently paper over.
    if "caseid" in filtered.columns:
        raw_caseid = filtered["caseid"].astype(str).str.strip()
        resolved_caseid = filtered["caseid_resolved"].astype(str).str.strip()
        n_mismatched = int((raw_caseid != resolved_caseid).sum())
        if n_mismatched:
            print(
                f"  [WARN] {table_name}: {n_mismatched:,} rows have a raw caseid that "
                f"disagrees with the DEMO-derived caseid for the same primaryid "
                f"(DEMO-derived caseid is authoritative and is what gets kept)"
            )
        filtered = filtered.drop(columns=["caseid"])

    filtered = filtered.rename(columns={"caseid_resolved": "caseid"})

    dropped = before - len(filtered)
    print(f"  {table_name}: {len(filtered):,} / {before:,} rows retained ({dropped:,} dropped; primaryid not retained)")
    return filtered


def compute_case_resolution_diagnostics(
    demo_df: pd.DataFrame,
    retained_cases_df: pd.DataFrame,
    drug_df: pd.DataFrame,
    resolved_drug_df: pd.DataFrame,
    reac_df: pd.DataFrame,
    resolved_reac_df: pd.DataFrame,
) -> dict:
    """Compute the Milestone 3 diagnostics: case resolution + DRUG/REAC filtering counts."""
    retained_primaryids = set(retained_cases_df["retained_primaryid"])

    caseid_not_missing = ~_is_missing(demo_df["caseid"])
    caseversion_bad = demo_df["caseversion_malformed"] | demo_df["caseversion_numeric"].isna()

    return {
        "n_demo_rows": int(len(demo_df)),
        "n_unique_raw_caseid": int(demo_df.loc[caseid_not_missing, "caseid"].nunique()),
        "n_unique_raw_primaryid": int(demo_df["primaryid"].nunique()),
        "n_retained_cases": int(len(retained_cases_df)),
        "n_retained_primaryids": int(retained_cases_df["retained_primaryid"].nunique()),
        "n_cases_with_multiple_versions": int((retained_cases_df["n_versions_for_case"] > 1).sum()),
        "max_versions_for_a_case": (
            int(retained_cases_df["n_versions_for_case"].max()) if len(retained_cases_df) else None
        ),
        "n_demo_rows_dropped_missing_caseid": int((~caseid_not_missing).sum()),
        "n_demo_rows_missing_or_malformed_caseversion": int(caseversion_bad.sum()),
        "n_drug_rows_raw": int(len(drug_df)),
        "n_drug_rows_resolved": int(len(resolved_drug_df)),
        "n_reac_rows_raw": int(len(reac_df)),
        "n_reac_rows_resolved": int(len(resolved_reac_df)),
        "n_drug_primaryids_not_retained": int(len(set(drug_df["primaryid"].unique()) - retained_primaryids)),
        "n_reac_primaryids_not_retained": int(len(set(reac_df["primaryid"].unique()) - retained_primaryids)),
    }


def print_case_resolution_diagnostics(diagnostics: dict) -> None:
    d = diagnostics
    print("\n  Case resolution diagnostics:")
    print(f"    raw DEMO rows:                          {d['n_demo_rows']:,}")
    print(f"    unique raw caseid:                      {d['n_unique_raw_caseid']:,}")
    print(f"    unique raw primaryid:                   {d['n_unique_raw_primaryid']:,}")
    print(f"    retained cases:                         {d['n_retained_cases']:,}")
    print(f"    retained primaryids:                    {d['n_retained_primaryids']:,}")
    print(f"    cases with >1 version:                  {d['n_cases_with_multiple_versions']:,}")
    print(f"    max versions for a single case:          {d['max_versions_for_a_case']}")
    print(f"    DEMO rows dropped (missing caseid):      {d['n_demo_rows_dropped_missing_caseid']:,}")
    print(f"    DEMO rows missing/malformed caseversion: {d['n_demo_rows_missing_or_malformed_caseversion']:,}")
    print(f"    raw DRUG rows:                           {d['n_drug_rows_raw']:,}")
    print(f"    resolved DRUG rows:                      {d['n_drug_rows_resolved']:,}")
    print(f"    raw REAC rows:                           {d['n_reac_rows_raw']:,}")
    print(f"    resolved REAC rows:                      {d['n_reac_rows_resolved']:,}")
    print(f"    DRUG primaryids not retained:            {d['n_drug_primaryids_not_retained']:,}")
    print(f"    REAC primaryids not retained:            {d['n_reac_primaryids_not_retained']:,}")


def run_case_resolution_validation(
    demo_df: pd.DataFrame,
    retained_cases_df: pd.DataFrame,
    resolved_drug_df: pd.DataFrame,
    resolved_reac_df: pd.DataFrame,
) -> None:
    """
    Validate the Milestone 3 case resolution and PRIMARYID filtering.

    Raises AssertionError listing every failed check rather than warning --
    a resolution bug here (e.g. a duplicated caseid, or a drug/reaction row
    that leaked in from a non-retained primaryid) would silently corrupt
    every later milestone, so these checks fail loudly instead.
    """
    failures = []

    if retained_cases_df["caseid"].duplicated().any():
        failures.append("retained_cases.caseid is not unique")

    if retained_cases_df["retained_primaryid"].isna().any():
        failures.append("retained_cases.retained_primaryid has null values")

    retained_primaryids = set(retained_cases_df["retained_primaryid"])
    if not set(resolved_drug_df["primaryid"]).issubset(retained_primaryids):
        failures.append("resolved_drug_rows contains a primaryid outside the retained set")
    if not set(resolved_reac_df["primaryid"]).issubset(retained_primaryids):
        failures.append("resolved_reac_rows contains a primaryid outside the retained set")

    if resolved_drug_df["caseid"].isna().any():
        failures.append("resolved_drug_rows has rows with null caseid")
    if resolved_reac_df["caseid"].isna().any():
        failures.append("resolved_reac_rows has rows with null caseid")

    # Retained caseversion should equal the max *valid* (non-malformed,
    # non-missing) caseversion for that case -- unless every version for
    # that case was missing/malformed, in which case there is no valid
    # value to compare against and the check is skipped for that case.
    valid_demo = demo_df[~demo_df["caseversion_malformed"] & demo_df["caseversion_numeric"].notna()]
    max_valid_by_case = valid_demo.groupby("caseid")["caseversion_numeric"].max()
    check_df = retained_cases_df.merge(
        max_valid_by_case.rename("max_valid_caseversion"), on="caseid", how="left"
    )
    has_valid_version = check_df["max_valid_caseversion"].notna()
    mismatched = has_valid_version & (check_df["caseversion"] != check_df["max_valid_caseversion"])
    if mismatched.any():
        failures.append(
            f"{int(mismatched.sum()):,} retained cases have caseversion != max valid caseversion for that case"
        )

    if failures:
        raise AssertionError("Milestone 3 validation failed:\n  - " + "\n  - ".join(failures))

    print("  All Milestone 3 validation checks passed:")
    print("    - retained_cases.caseid is unique")
    print("    - retained_cases.retained_primaryid has no nulls")
    print("    - resolved_drug_rows.primaryid subset of retained primaryids")
    print("    - resolved_reac_rows.primaryid subset of retained primaryids")
    print("    - resolved_drug_rows.caseid has no nulls")
    print("    - resolved_reac_rows.caseid has no nulls")
    print("    - retained caseversion == max valid caseversion per case (where one exists)")


def save_case_resolution_artifacts(
    out_dir: str,
    retained_cases_df: pd.DataFrame,
    resolved_drug_df: pd.DataFrame,
    resolved_reac_df: pd.DataFrame,
    diagnostics: dict,
) -> None:
    retained_path = os.path.join(out_dir, "retained_cases.parquet")
    retained_cases_df.to_parquet(retained_path, index=False)
    print(f"  Saved: {retained_path}  ({len(retained_cases_df):,} rows)")

    drug_path = os.path.join(out_dir, "resolved_drug_rows.parquet")
    resolved_drug_df.to_parquet(drug_path, index=False)
    print(f"  Saved: {drug_path}  ({len(resolved_drug_df):,} rows)")

    reac_path = os.path.join(out_dir, "resolved_reac_rows.parquet")
    resolved_reac_df.to_parquet(reac_path, index=False)
    print(f"  Saved: {reac_path}  ({len(resolved_reac_df):,} rows)")

    diagnostics_path = os.path.join(out_dir, "case_resolution_diagnostics.json")
    with open(diagnostics_path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    print(f"  Saved: {diagnostics_path}")


# ---------------------------------------------------------------------------
# Milestone 4 -- dedupe reactions per case (case_reactions.parquet)
# ---------------------------------------------------------------------------

def build_case_reactions(resolved_reac_df: pd.DataFrame) -> tuple:
    """
    Build the case-level adverse-event table from resolved_reac_rows (which
    already carries the authoritative, DEMO-derived caseid).

    Steps:
      1. Normalize pt -> reaction_pt: lowercase, strip, collapse repeated
         internal whitespace; empty strings become missing.
      2. Drop rows with missing caseid or missing/blank reaction_pt.
      3. Optionally drop administrative/non-clinical PTs using
         process_faers._SE_EXCLUDE (reused read-only, PT filtering only --
         no row/primaryid-counting logic from the old pipeline is reused).
      4. Deduplicate to one row per (caseid, reaction_pt) -- a case
         contributes a given MedDRA PT at most once. source_primaryid keeps
         whichever primaryid sorts first within the case for provenance
         only; it has no bearing on which reaction_pt values survive.

    Returns (case_reactions_df, diagnostics_dict).
    """
    _require_columns(resolved_reac_df, ["caseid", "primaryid", "pt"], "resolved_reac_rows")

    n_input_rows = len(resolved_reac_df)
    working = resolved_reac_df.copy()

    # Step 1: normalize pt -> reaction_pt
    normalized = (
        working["pt"].fillna("").astype(str).str.strip().str.lower()
        .str.replace(r"\s+", " ", regex=True)
    )
    working["reaction_pt"] = normalized.replace("", pd.NA)

    # Step 2: drop missing caseid / missing reaction_pt
    missing_caseid_mask = _is_missing(working["caseid"])
    n_missing_caseid = int(missing_caseid_mask.sum())

    missing_pt_mask = working["reaction_pt"].isna()
    n_missing_pt = int(missing_pt_mask.sum())

    working = working[~missing_caseid_mask & ~missing_pt_mask].copy()

    # Step 3: optional administrative/non-clinical PT exclusion
    exclusion_reused = _SE_EXCLUDE is not None
    n_excluded = 0
    if exclusion_reused:
        excluded_mask = working["reaction_pt"].isin(_SE_EXCLUDE)
        n_excluded = int(excluded_mask.sum())
        print(
            f"  Reused process_faers._SE_EXCLUDE ({len(_SE_EXCLUDE)} PTs, PT filtering only): "
            f"removed {n_excluded:,} administrative/non-clinical rows"
        )
        working = working[~excluded_mask].copy()
    else:
        print(
            "  [TODO] process_faers._SE_EXCLUDE unavailable -- no administrative/non-clinical "
            "PT exclusion list applied (not inventing a new one)"
        )

    n_before_dedup = len(working)

    # Step 4: one row per (caseid, reaction_pt)
    working = working.sort_values(["caseid", "reaction_pt", "primaryid"], kind="mergesort")
    deduped = working.drop_duplicates(subset=["caseid", "reaction_pt"], keep="first")

    case_reactions_df = (
        deduped.rename(columns={"primaryid": "source_primaryid", "pt": "pt_original"})
        [["caseid", "reaction_pt", "source_primaryid", "pt_original"]]
        .reset_index(drop=True)
    )

    n_after_dedup = len(case_reactions_df)

    diagnostics = {
        "n_input_rows": n_input_rows,
        "n_missing_caseid": n_missing_caseid,
        "n_missing_or_blank_pt": n_missing_pt,
        "exclusion_list_reused": exclusion_reused,
        "n_removed_by_exclusion_list": n_excluded,
        "n_rows_before_dedup": n_before_dedup,
        "n_rows_after_dedup": n_after_dedup,
        "n_duplicate_rows_removed": n_before_dedup - n_after_dedup,
        "n_unique_cases_with_reactions": int(case_reactions_df["caseid"].nunique()),
        "n_unique_reaction_pt": int(case_reactions_df["reaction_pt"].nunique()),
        "top_20_reaction_pt": {
            str(k): int(v) for k, v in case_reactions_df["reaction_pt"].value_counts().head(20).items()
        },
    }

    print(
        f"  case_reactions: {n_after_dedup:,} rows "
        f"({diagnostics['n_unique_cases_with_reactions']:,} cases, "
        f"{diagnostics['n_unique_reaction_pt']:,} unique reaction_pt)"
    )

    return case_reactions_df, diagnostics


def print_reaction_dedup_diagnostics(diagnostics: dict) -> None:
    d = diagnostics
    print("\n  Reaction dedup diagnostics:")
    print(f"    input resolved REAC rows:       {d['n_input_rows']:,}")
    print(f"    rows with missing caseid:       {d['n_missing_caseid']:,}")
    print(f"    rows missing/blank PT:          {d['n_missing_or_blank_pt']:,}")
    print(f"    exclusion list reused:          {d['exclusion_list_reused']}")
    print(f"    rows removed by exclusion list: {d['n_removed_by_exclusion_list']:,}")
    print(f"    rows before dedup:              {d['n_rows_before_dedup']:,}")
    print(f"    rows after dedup:               {d['n_rows_after_dedup']:,}")
    print(f"    duplicate rows removed:         {d['n_duplicate_rows_removed']:,}")
    print(f"    unique cases with reactions:    {d['n_unique_cases_with_reactions']:,}")
    print(f"    unique reaction_pt:             {d['n_unique_reaction_pt']:,}")
    print("    top 20 reaction_pt:")
    for k, v in d["top_20_reaction_pt"].items():
        print(f"      {k:<40} {v:,}")


def run_case_reactions_validation(resolved_reac_df: pd.DataFrame, case_reactions_df: pd.DataFrame) -> None:
    """
    Validate case_reactions. Raises AssertionError listing every failed
    check -- a bug here would silently corrupt drug-pair/reaction counting
    in later milestones, so it fails loudly instead of warning.
    """
    failures = []

    required_cols = ["caseid", "reaction_pt", "source_primaryid"]
    missing_cols = [c for c in required_cols if c not in case_reactions_df.columns]
    if missing_cols:
        failures.append(f"case_reactions is missing required column(s): {missing_cols}")
    else:
        if case_reactions_df["caseid"].isna().any():
            failures.append("case_reactions has rows with null caseid")
        if case_reactions_df["reaction_pt"].isna().any():
            failures.append("case_reactions has rows with null reaction_pt")
        if case_reactions_df.duplicated(subset=["caseid", "reaction_pt"]).any():
            failures.append("case_reactions has duplicate (caseid, reaction_pt) rows")
        valid_caseids = set(resolved_reac_df["caseid"])
        if not set(case_reactions_df["caseid"]).issubset(valid_caseids):
            failures.append("case_reactions contains a caseid not present in resolved_reac_rows")

    if failures:
        raise AssertionError("Milestone 4 (case_reactions) validation failed:\n  - " + "\n  - ".join(failures))

    print("  All Milestone 4 (case_reactions) validation checks passed:")
    print("    - case_reactions has required columns")
    print("    - no missing caseid")
    print("    - no missing reaction_pt")
    print("    - no duplicate (caseid, reaction_pt)")
    print("    - every case_reactions.caseid present in resolved_reac_rows")


def save_case_reactions_artifacts(out_dir: str, case_reactions_df: pd.DataFrame, diagnostics: dict) -> None:
    path = os.path.join(out_dir, "case_reactions.parquet")
    case_reactions_df.to_parquet(path, index=False)
    print(f"  Saved: {path}  ({len(case_reactions_df):,} rows)")

    diagnostics_path = os.path.join(out_dir, "reaction_dedup_diagnostics.json")
    with open(diagnostics_path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    print(f"  Saved: {diagnostics_path}")


# ---------------------------------------------------------------------------
# Milestone 4b (TODO) -- dedupe drugs per case
# ---------------------------------------------------------------------------

def build_case_drugs(drug_df: pd.DataFrame, case_map_df: pd.DataFrame) -> pd.DataFrame:
    """
    TODO (Milestone 4b): Given resolved_drug_rows (already filtered to
    retained primaryids and carrying caseid), normalize drug names
    (normalize_drug_string) and drop duplicate (caseid, normalized_drug)
    rows, merging role_cod where needed.
    Output columns: caseid, drug, role_cod.  Saved later as case_drugs.parquet.
    """
    print("  [TODO] build_case_drugs() not yet implemented (Milestone 4b)")
    return pd.DataFrame(columns=["caseid", "drug", "role_cod"])


# ---------------------------------------------------------------------------
# Milestone 5 (TODO) -- canonical drug-pair generation per case
# ---------------------------------------------------------------------------

def build_case_drug_pairs(case_drugs_df: pd.DataFrame) -> pd.DataFrame:
    """
    TODO (Milestone 5): Per case, enumerate canonical (drug_i, drug_j) pairs
    with drug_i < drug_j (no self-pairs, no direction duplicates).
    Output columns: caseid, drug_i, drug_j.  Saved later as case_drug_pairs.parquet.
    """
    print("  [TODO] build_case_drug_pairs() not yet implemented (Milestone 5)")
    return pd.DataFrame(columns=["caseid", "drug_i", "drug_j"])


# ---------------------------------------------------------------------------
# Milestone 6 (TODO) -- drug_pair_reaction_counts aggregation (no PRR yet)
# ---------------------------------------------------------------------------

def build_drug_pair_reaction_counts(
    case_drug_pairs_df: pd.DataFrame,
    case_reactions_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    TODO (Milestone 6): Join case_drug_pairs with case_reactions on caseid and
    count distinct cases per (drug_i, drug_j, pt). This is a count table only;
    PRR / disproportionality statistics are explicitly out of scope here.
    Output columns: drug_i, drug_j, pt, n_cases.
    Saved later as drug_pair_reaction_counts.parquet.
    """
    print("  [TODO] build_drug_pair_reaction_counts() not yet implemented (Milestone 6)")
    return pd.DataFrame(columns=["drug_i", "drug_j", "pt", "n_cases"])


# ---------------------------------------------------------------------------
# Milestone 7 (TODO) -- global validation checks / PRR
# ---------------------------------------------------------------------------

def run_validation_checks(
    demo_df: pd.DataFrame,
    case_map_df: pd.DataFrame,
    case_reactions_df: pd.DataFrame,
    case_drugs_df: pd.DataFrame,
) -> None:
    """
    TODO (Milestone 7): Sanity checks, e.g.:
      - distinct caseid count in case_map_df matches distinct caseid in DEMO
      - no primaryid appears under two different caseids
      - no case has zero reactions or zero drugs after filtering
      - case_drug_pairs has no (drug_i == drug_j) rows and no reversed duplicates
    """
    print("  [TODO] run_validation_checks() not yet implemented (Milestone 7)")


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
            "once, ahead of drug-only R-GCN construction. Milestone 4: "
            "builds case_reactions.parquet (one row per caseid, reaction_pt) "
            "-- see file docstring for milestone status."
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

    print("\n--- Milestone 3: resolve CASEID/CASEVERSION ---")
    retained_cases_df = resolve_case_versions(demo_df)

    print("\n--- Milestone 3: filter DRUG/REAC to retained PRIMARYIDs ---")
    resolved_drug_df = filter_to_retained_primaryids(drug_df, retained_cases_df, "DRUG")
    resolved_reac_df = filter_to_retained_primaryids(reac_df, retained_cases_df, "REAC")

    print("\n--- Milestone 3: case-resolution diagnostics ---")
    case_resolution_diagnostics = compute_case_resolution_diagnostics(
        demo_df, retained_cases_df, drug_df, resolved_drug_df, reac_df, resolved_reac_df,
    )
    print_case_resolution_diagnostics(case_resolution_diagnostics)

    print("\n--- Milestone 3: validation checks ---")
    run_case_resolution_validation(demo_df, retained_cases_df, resolved_drug_df, resolved_reac_df)

    print("\n--- Milestone 3: saving case-resolution artifacts ---")
    save_case_resolution_artifacts(
        args.out_dir, retained_cases_df, resolved_drug_df, resolved_reac_df, case_resolution_diagnostics,
    )

    print("\n--- Milestone 4: build case_reactions ---")
    case_reactions_df, reaction_dedup_diagnostics = build_case_reactions(resolved_reac_df)
    print_reaction_dedup_diagnostics(reaction_dedup_diagnostics)

    print("\n--- Milestone 4: case_reactions validation checks ---")
    run_case_reactions_validation(resolved_reac_df, case_reactions_df)

    print("\n--- Milestone 4: saving case_reactions artifacts ---")
    save_case_reactions_artifacts(args.out_dir, case_reactions_df, reaction_dedup_diagnostics)

    print("\n--- Milestone 4b (pending): dedupe drugs per case ---")
    case_drugs_df = build_case_drugs(resolved_drug_df, retained_cases_df)

    print("\n--- Milestone 5 (pending): canonical drug-pair generation ---")
    case_drug_pairs_df = build_case_drug_pairs(case_drugs_df)

    print("\n--- Milestone 6 (pending): drug_pair_reaction_counts aggregation ---")
    drug_pair_reaction_counts_df = build_drug_pair_reaction_counts(case_drug_pairs_df, case_reactions_df)

    print("\n--- Milestone 7 (pending): global validation checks ---")
    run_validation_checks(demo_df, retained_cases_df, case_reactions_df, case_drugs_df)

    print("\n--- Save (pending): case-level artifacts ---")
    save_case_level_artifacts(
        args.out_dir, case_reactions_df, case_drugs_df,
        case_drug_pairs_df, drug_pair_reaction_counts_df,
    )

    print("\n" + "=" * 60)
    print("Milestone 4 complete: case_reactions.parquet built (one row per")
    print("caseid, reaction_pt).")
    print("No drug deduplication, drug mapping, or drug-pair generation")
    print("was performed yet.")
    print(f"Artifacts in: {args.out_dir}")
    print("Next: implement Milestone 4b (dedupe drugs per case).")
    print("=" * 60)


if __name__ == "__main__":
    main()
