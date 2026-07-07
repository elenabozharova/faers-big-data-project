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
  [x] Milestone 4 -- dedupe reactions per case (case_reactions.parquet)
  [x] Milestone 4b -- map + dedupe drugs per case (case_drugs.parquet)
  [x] Milestone 4c -- selectable drug ID mode (decagon_idx / chembl_id)        (this file)
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

Milestone 4 scope
------------------
Builds case_reactions.parquet: one row per (caseid, reaction_pt). Reuses
process_faers._SE_EXCLUDE (administrative/non-clinical PT strings) read-only,
purely for PT filtering -- no old row/primaryid-counting logic is reused.

Milestone 4b scope
-------------------
Builds case_drugs_long.parquet (row-level mapping provenance) and
case_drugs.parquet (one row per caseid, mapped_drug_id). Reuses, read-only:
  - the Tier 0 raw-name selection rule (prod_ai vs drugname) and Tier 1
    normalize_drug_string(), both from process_faers.py's
    map_drug_rows_to_decagon()
  - the FAERS->ChEMBL mapping file and the Decagon drug_to_idx.json universe
Does NOT reuse build_name_to_decagon_lookup()/map_drug_rows_to_decagon()
themselves, since those collapse "why a drug didn't map" into aggregate
counters -- this milestone needs a per-row mapping_status, so the row-level
join against the mapping table is reimplemented here (same underlying
rules: no ChEMBL id / ambiguous / fuzzy / outside Decagon universe).

Milestone 4c scope (this file)
-------------------------------
Adds --drug-id-mode, so `mapped_drug_id` can be either:
  - decagon_idx (default, backward-compatible with Milestone 4b): a row is
    only "mapped" if its ChEMBL id also falls inside the old Decagon
    drug_to_idx.json universe; mapped_drug_id is that integer index.
  - chembl_id (recommended for the drug-only R-GCN direction the thesis is
    now taking): mapped_drug_id is the ChEMBL id string itself, with no
    Decagon-universe restriction. A row is "mapped" as soon as it has a
    valid, unambiguous ChEMBL id that passes the fuzzy-inclusion rule --
    `unmapped_outside_drug_universe` essentially disappears in this mode.
Both modes share the same Tier 0/Tier 1 name resolution and the same
FAERS->ChEMBL mapping file; only the "is this row mapped, and to what ID"
step differs. `drug_id_mode` is stamped onto every output row so a run's
provenance is unambiguous even if two runs' artifacts are later combined.
Still explicitly out of scope:
  - drug-pair generation (Milestone 5)
  - PRR / disproportionality statistics (Milestone 7)

This file does NOT:
  - implement full preprocessing
  - generate drug pairs
  - compute PRR (Proportional Reporting Ratio)
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

# Milestone 4c: how mapped_drug_id is derived. 'decagon_idx' is the default
# for backward compatibility with Milestone 4b; 'chembl_id' is recommended
# for the drug-only R-GCN direction (no Decagon universe restriction).
DRUG_ID_MODES = ("decagon_idx", "chembl_id")
DEFAULT_DRUG_ID_MODE = "decagon_idx"
_CHEMBL_ID_RE = re.compile(r"^CHEMBL\d+$")


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
# Milestone 4b -- map + dedupe drugs per case (case_drugs.parquet)
# ---------------------------------------------------------------------------

# Drug universe used for mapped_drug_id. Mirrors process_faers.py's default
# (--mode all, the broad ChEMBL activity-based universe); this pipeline has
# no --mode flag yet.
_DRUG_UNIVERSE_MODE = "all"


def _choose_raw_drug_name(drug_df: pd.DataFrame) -> pd.Series:
    """
    Tier 0 raw-name selection -- identical rule to
    process_faers.map_drug_rows_to_decagon(): prefer prod_ai when present,
    non-empty (>2 chars), and not purely numeric; otherwise fall back to
    drugname.
    """
    prod_ai = (
        drug_df["prod_ai"].fillna("").astype(str).str.strip()
        if "prod_ai" in drug_df.columns
        else pd.Series([""] * len(drug_df), index=drug_df.index, dtype=str)
    )
    drugname = (
        drug_df["drugname"].fillna("").astype(str).str.strip()
        if "drugname" in drug_df.columns
        else pd.Series([""] * len(drug_df), index=drug_df.index, dtype=str)
    )
    is_meaningful = (prod_ai.str.len() > 2) & (~prod_ai.str.replace(".", "", regex=False).str.isnumeric())
    return prod_ai.where(is_meaningful, drugname)


def load_drug_mapping_file(mapping_file: str) -> pd.DataFrame:
    """
    Load the FAERS drug name -> ChEMBL mapping parquet from an explicit file
    path (this pipeline's --mapping-file is a full path, unlike
    process_faers.load_drug_mapping() which takes a containing directory).
    """
    if not os.path.exists(mapping_file):
        raise FileNotFoundError(
            f"Drug mapping not found at:\n  {mapping_file}\n"
            "Run map_faers_to_chembl.py first to produce this file."
        )
    df = pd.read_parquet(mapping_file)
    print(f"  Loaded drug mapping: {len(df):,} rows from {mapping_file}")
    n_dupe_normalized = int(df["normalized"].duplicated().sum())
    if n_dupe_normalized:
        print(f"  [WARN] drug mapping has {n_dupe_normalized:,} duplicate 'normalized' rows; keeping first")
        df = df.drop_duplicates(subset="normalized", keep="first")
    return df


def load_drug_universe(mode: str = _DRUG_UNIVERSE_MODE) -> dict:
    """
    Load the Decagon drug_to_idx.json for the given universe mode.

    mode='strict' -> data/processed/decagon/drug_to_idx.json   (mechanism-based, smaller)
    mode='all'    -> data/processed/decagon_all/drug_to_idx.json (broad activity-based; default here)
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


def build_case_drugs(
    resolved_drug_df: pd.DataFrame,
    mapping_file: str,
    drug_id_mode: str = DEFAULT_DRUG_ID_MODE,
    include_fuzzy: bool = False,
) -> tuple:
    """
    Build case-level drug exposure tables from resolved_drug_rows (which
    already carries the authoritative, DEMO-derived caseid).

    Per row: choose chosen_raw_name (Tier 0), normalize it (Tier 1) into
    normalized_drug_name, then classify mapping_status as one of:
      mapped / unmapped_blank_name / unmapped_no_match /
      unmapped_ambiguous / excluded_fuzzy / unmapped_outside_drug_universe
    (This is more granular than process_faers.py's build_name_to_decagon_lookup,
    which lumps "no chembl_id" and "ambiguous" into a single counter.)

    drug_id_mode (Milestone 4c) controls what mapped_drug_id is and whether
    the Decagon universe restricts "mapped" at all:
      - 'decagon_idx': mapped_drug_id is the Decagon drug_to_idx.json index;
        a row outside that universe is 'unmapped_outside_drug_universe'.
      - 'chembl_id': mapped_drug_id is the ChEMBL id itself; there is no
        Decagon universe restriction, so 'unmapped_outside_drug_universe'
        will not occur (every row with a valid, unambiguous, non-fuzzy-
        excluded ChEMBL id is mapped).

    Returns (case_drugs_long_df, case_drugs_df, diagnostics_dict, unmapped_drugs_top_df).
    """
    if drug_id_mode not in DRUG_ID_MODES:
        raise ValueError(f"drug_id_mode must be one of {DRUG_ID_MODES}, got {drug_id_mode!r}")

    _require_columns(resolved_drug_df, ["caseid", "primaryid"], "resolved_drug_rows")

    n_input_rows = len(resolved_drug_df)
    working = resolved_drug_df.reset_index(drop=True).copy()
    for col in ("drugname", "prod_ai", "role_cod"):
        if col not in working.columns:
            working[col] = ""

    mapping_df = load_drug_mapping_file(mapping_file)

    decagon_drug_idx = None
    decagon_chembl_ids = set()
    if drug_id_mode == "decagon_idx":
        decagon_drug_idx = load_drug_universe(_DRUG_UNIVERSE_MODE)
        decagon_chembl_ids = set(decagon_drug_idx.keys())
    else:
        print("  drug_id_mode='chembl_id': skipping Decagon drug_to_idx.json universe restriction")

    working["chosen_raw_name"] = _choose_raw_drug_name(working)

    print("  Normalizing drug names (Tier 1)...")
    blank_name_mask = working["chosen_raw_name"].str.strip() == ""
    working["normalized_drug_name"] = pd.NA
    working.loc[~blank_name_mask, "normalized_drug_name"] = (
        working.loc[~blank_name_mask, "chosen_raw_name"].apply(normalize_drug_string)
    )

    mapping_lookup = mapping_df[["normalized", "chembl_id", "is_ambiguous", "match_method"]].rename(
        columns={"normalized": "normalized_drug_name"}
    )
    working = working.merge(mapping_lookup, on="normalized_drug_name", how="left", indicator="_merge_status")

    blank_name_mask = working["chosen_raw_name"].str.strip() == ""  # recompute post-merge for safe alignment
    has_match = working["_merge_status"] == "both"
    has_chembl_id = has_match & working["chembl_id"].notna()
    is_ambiguous = has_chembl_id & working["is_ambiguous"].fillna(False).astype(bool)
    is_fuzzy = has_chembl_id & (working["match_method"] == "fuzzy")
    if drug_id_mode == "decagon_idx":
        in_universe = has_chembl_id & working["chembl_id"].isin(decagon_chembl_ids)
    else:
        in_universe = has_chembl_id  # chembl_id mode: no Decagon universe restriction

    mapping_status = pd.Series("mapped", index=working.index, dtype=object)
    mapping_status[blank_name_mask] = "unmapped_blank_name"
    mapping_status[~blank_name_mask & ~has_chembl_id] = "unmapped_no_match"
    mapping_status[~blank_name_mask & has_chembl_id & is_ambiguous] = "unmapped_ambiguous"
    excluded_fuzzy_mask = ~blank_name_mask & has_chembl_id & ~is_ambiguous & is_fuzzy & (not include_fuzzy)
    mapping_status[excluded_fuzzy_mask] = "excluded_fuzzy"
    still_candidate = ~blank_name_mask & has_chembl_id & ~is_ambiguous & (include_fuzzy | ~is_fuzzy)
    mapping_status[still_candidate & ~in_universe] = "unmapped_outside_drug_universe"

    working["mapping_status"] = mapping_status
    mapped_mask = working["mapping_status"] == "mapped"

    working["mapped_drug_id"] = pd.NA
    if drug_id_mode == "decagon_idx":
        working.loc[mapped_mask, "mapped_drug_id"] = working.loc[mapped_mask, "chembl_id"].map(decagon_drug_idx)
    else:
        working.loc[mapped_mask, "mapped_drug_id"] = working.loc[mapped_mask, "chembl_id"]
    working["mapping_method"] = working["match_method"]
    working["drug_id_mode"] = drug_id_mode

    case_drugs_long_df = working[
        [
            "caseid", "primaryid", "role_cod", "drugname", "prod_ai",
            "chosen_raw_name", "normalized_drug_name", "mapped_drug_id",
            "mapping_status", "mapping_method", "chembl_id", "drug_id_mode",
        ]
    ].reset_index(drop=True)

    n_mapped = int(mapped_mask.sum())
    n_unmapped = n_input_rows - n_mapped
    print(
        f"  case_drugs_long: {len(case_drugs_long_df):,} rows -- "
        f"{n_mapped:,} mapped ({100 * n_mapped / n_input_rows:.1f}%), {n_unmapped:,} unmapped"
    )

    # ---- case_drugs: one row per (caseid, mapped_drug_id), mapped rows only ----
    mapped_rows = case_drugs_long_df[case_drugs_long_df["mapping_status"] == "mapped"].copy()
    grouped = mapped_rows.groupby(["caseid", "mapped_drug_id"], dropna=False)
    case_drugs_df = grouped.agg(
        representative_name=("normalized_drug_name", "first"),
        roles_seen=("role_cod", lambda s: sorted(set(s.dropna().astype(str)) - {""})),
        mapping_methods_seen=("mapping_method", lambda s: sorted(set(s.dropna().astype(str)))),
        n_raw_rows_for_case_drug=("primaryid", "size"),
    ).reset_index()
    case_drugs_df["drug_id_mode"] = drug_id_mode

    n_dedup_removed = len(mapped_rows) - len(case_drugs_df)
    print(
        f"  case_drugs: {len(case_drugs_df):,} rows "
        f"({case_drugs_df['caseid'].nunique():,} cases, "
        f"{case_drugs_df['mapped_drug_id'].nunique():,} unique mapped drugs; "
        f"{n_dedup_removed:,} duplicate (caseid, mapped_drug_id) rows removed)"
    )

    # ---- unmapped_drugs_top ----
    unmapped_rows = case_drugs_long_df[case_drugs_long_df["mapping_status"] != "mapped"]
    unmapped_drugs_top_df = (
        unmapped_rows.groupby(["normalized_drug_name", "chosen_raw_name", "mapping_status"], dropna=False)
        .size()
        .reset_index(name="n_rows")
        .sort_values("n_rows", ascending=False)
        .reset_index(drop=True)
    )

    # ---- cases with >=2 distinct raw (normalized) drugs but <2 mapped drugs ----
    raw_drug_count_per_case = (
        working.loc[~blank_name_mask].groupby("caseid")["normalized_drug_name"].nunique()
    )
    mapped_drug_count_per_case = case_drugs_df.groupby("caseid")["mapped_drug_id"].nunique()
    combined = pd.DataFrame({"raw": raw_drug_count_per_case}).join(
        pd.DataFrame({"mapped": mapped_drug_count_per_case}), how="left"
    )
    combined["mapped"] = combined["mapped"].fillna(0)
    n_cases_lost_multi_drug = int(((combined["raw"] >= 2) & (combined["mapped"] < 2)).sum())

    mapped_drug_count_dist = mapped_drug_count_per_case.value_counts().sort_index()
    n_cases_with_ge2_mapped = int((mapped_drug_count_per_case >= 2).sum())

    diagnostics = {
        "drug_id_mode": drug_id_mode,
        "n_input_rows": n_input_rows,
        "n_blank_chosen_name": int(blank_name_mask.sum()),
        "n_mapped_rows": n_mapped,
        "n_unmapped_rows": n_unmapped,
        "mapping_rate": (n_mapped / n_input_rows) if n_input_rows else None,
        "n_unique_raw_cases": int(working["caseid"].nunique()),
        "n_unique_cases_with_mapped_drug": int(case_drugs_df["caseid"].nunique()),
        "n_unique_chembl_ids_mapped": int(working.loc[mapped_mask, "chembl_id"].nunique()),
        "n_unique_mapped_drugs": int(case_drugs_df["mapped_drug_id"].nunique()),
        "include_fuzzy": include_fuzzy,
        "n_removed_fuzzy_excluded": int((mapping_status == "excluded_fuzzy").sum()),
        # Only relevant in drug_id_mode='decagon_idx'; structurally 0 in 'chembl_id' mode
        # since that mode applies no Decagon-universe restriction.
        "n_removed_outside_drug_universe": int((mapping_status == "unmapped_outside_drug_universe").sum()),
        "n_removed_ambiguous": int((mapping_status == "unmapped_ambiguous").sum()),
        "n_removed_no_match": int((mapping_status == "unmapped_no_match").sum()),
        "n_duplicate_case_drug_rows_removed": int(n_dedup_removed),
        "n_cases_with_ge2_mapped_drugs": n_cases_with_ge2_mapped,
        "n_cases_with_raw_ge2_mapped_lt2": n_cases_lost_multi_drug,
        "mapped_drug_count_per_case_distribution": {
            str(int(k)): int(v) for k, v in mapped_drug_count_dist.items()
        },
    }

    return case_drugs_long_df, case_drugs_df, diagnostics, unmapped_drugs_top_df


def print_drug_mapping_diagnostics(diagnostics: dict) -> None:
    d = diagnostics
    print("\n  Drug mapping diagnostics:")
    print(f"    drug_id_mode:                          {d['drug_id_mode']}")
    print(f"    input resolved DRUG rows:             {d['n_input_rows']:,}")
    print(f"    rows with blank chosen name:           {d['n_blank_chosen_name']:,}")
    print(f"    mapped rows:                           {d['n_mapped_rows']:,}")
    print(f"    unmapped rows:                          {d['n_unmapped_rows']:,}")
    rate = d["mapping_rate"]
    print(f"    mapping rate:                          {rate:.1%}" if rate is not None else "    mapping rate:                          n/a")
    print(f"    unique raw cases:                      {d['n_unique_raw_cases']:,}")
    print(f"    unique cases with >=1 mapped drug:     {d['n_unique_cases_with_mapped_drug']:,}")
    print(f"    unique ChEMBL ids mapped:              {d['n_unique_chembl_ids_mapped']:,}")
    print(f"    unique mapped drug ids:                {d['n_unique_mapped_drugs']:,}")
    print(f"    include_fuzzy:                         {d['include_fuzzy']}")
    print(f"    removed (fuzzy excluded):              {d['n_removed_fuzzy_excluded']:,}")
    universe_note = "" if d["drug_id_mode"] == "decagon_idx" else "  (n/a in chembl_id mode)"
    print(f"    removed (outside drug universe):       {d['n_removed_outside_drug_universe']:,}{universe_note}")
    print(f"    removed (ambiguous):                   {d['n_removed_ambiguous']:,}")
    print(f"    removed (no match):                    {d['n_removed_no_match']:,}")
    print(f"    duplicate (caseid, mapped_drug_id):    {d['n_duplicate_case_drug_rows_removed']:,}")
    print(f"    cases with >=2 mapped drugs:           {d['n_cases_with_ge2_mapped_drugs']:,}")
    print(f"    cases with raw>=2 but mapped<2:         {d['n_cases_with_raw_ge2_mapped_lt2']:,}")
    print("    mapped-drug-count-per-case distribution:")
    for k, v in d["mapped_drug_count_per_case_distribution"].items():
        print(f"      {k:>3} drug(s): {v:,} cases")


def run_case_drugs_validation(
    resolved_drug_df: pd.DataFrame,
    case_drugs_long_df: pd.DataFrame,
    case_drugs_df: pd.DataFrame,
    drug_id_mode: str,
) -> None:
    """
    Validate case_drugs_long/case_drugs. Raises AssertionError listing every
    failed check -- a mapping/dedup bug here would silently corrupt
    drug-pair generation in later milestones, so it fails loudly.
    """
    failures = []

    required_long_cols = [
        "caseid", "primaryid", "role_cod", "drugname", "prod_ai",
        "chosen_raw_name", "normalized_drug_name", "mapped_drug_id",
        "mapping_status", "mapping_method", "chembl_id", "drug_id_mode",
    ]
    missing_long_cols = [c for c in required_long_cols if c not in case_drugs_long_df.columns]
    if missing_long_cols:
        failures.append(f"case_drugs_long is missing required column(s): {missing_long_cols}")
    else:
        long_modes = set(case_drugs_long_df["drug_id_mode"].dropna().unique())
        if long_modes != {drug_id_mode}:
            failures.append(
                f"case_drugs_long.drug_id_mode has value(s) {sorted(long_modes)}, expected only {drug_id_mode!r}"
            )

    if "drug_id_mode" not in case_drugs_df.columns:
        failures.append("case_drugs is missing required column: drug_id_mode")
    else:
        short_modes = set(case_drugs_df["drug_id_mode"].dropna().unique())
        if short_modes and short_modes != {drug_id_mode}:
            failures.append(
                f"case_drugs.drug_id_mode has value(s) {sorted(short_modes)}, expected only {drug_id_mode!r}"
            )

    if case_drugs_df["caseid"].isna().any():
        failures.append("case_drugs has rows with null caseid")
    if case_drugs_df["mapped_drug_id"].isna().any():
        failures.append("case_drugs has rows with null mapped_drug_id")
    if case_drugs_df.duplicated(subset=["caseid", "mapped_drug_id"]).any():
        failures.append("case_drugs has duplicate (caseid, mapped_drug_id) rows")

    valid_caseids = set(resolved_drug_df["caseid"])
    if not set(case_drugs_df["caseid"]).issubset(valid_caseids):
        failures.append("case_drugs contains a caseid not present in resolved_drug_rows")

    mapped_long = case_drugs_long_df[case_drugs_long_df["mapping_status"] == "mapped"]
    mapped_pairs_long = set(map(tuple, mapped_long[["caseid", "mapped_drug_id"]].itertuples(index=False, name=None)))
    case_drugs_pairs = set(map(tuple, case_drugs_df[["caseid", "mapped_drug_id"]].itertuples(index=False, name=None)))
    if not mapped_pairs_long.issubset(case_drugs_pairs):
        failures.append(
            "some mapped (caseid, mapped_drug_id) pairs in case_drugs_long are missing from case_drugs after dedup"
        )

    mapped_ids = case_drugs_df["mapped_drug_id"].dropna()
    if drug_id_mode == "chembl_id":
        bad_ids = mapped_ids[~mapped_ids.astype(str).str.match(_CHEMBL_ID_RE.pattern)]
        if len(bad_ids):
            failures.append(
                f"{len(bad_ids):,} mapped_drug_id values in chembl_id mode don't look like ChEMBL ids "
                f"(expected 'CHEMBL<digits>'), e.g. {sorted(bad_ids.astype(str).unique())[:5]}"
            )
    else:  # decagon_idx
        def _is_decagon_index(v) -> bool:
            if isinstance(v, bool):
                return False
            if isinstance(v, (int,)):
                return True
            try:
                return float(v).is_integer()
            except (TypeError, ValueError):
                return False

        bad_ids = mapped_ids[~mapped_ids.apply(_is_decagon_index)]
        if len(bad_ids):
            failures.append(
                f"{len(bad_ids):,} mapped_drug_id values in decagon_idx mode don't look like Decagon "
                f"drug_to_idx.json integer indices, e.g. {sorted(bad_ids.astype(str).unique())[:5]}"
            )

    if failures:
        raise AssertionError("Milestone 4c (case_drugs) validation failed:\n  - " + "\n  - ".join(failures))

    print("  All Milestone 4c (case_drugs) validation checks passed:")
    print("    - case_drugs_long has required provenance columns")
    print("    - case_drugs_long.drug_id_mode has exactly one value, matching the CLI mode")
    print("    - case_drugs.drug_id_mode present and matches the CLI mode")
    print("    - case_drugs has no missing caseid")
    print("    - case_drugs has no missing mapped_drug_id")
    print("    - case_drugs has no duplicate (caseid, mapped_drug_id)")
    print("    - every case_drugs.caseid present in resolved_drug_rows")
    print("    - every mapped (caseid, mapped_drug_id) in case_drugs_long survives into case_drugs")
    print(f"    - mapped_drug_id values look like the expected ID type for drug_id_mode={drug_id_mode!r}")


# ---------------------------------------------------------------------------
# Milestone 4c -- decagon_idx vs chembl_id comparison helper
# ---------------------------------------------------------------------------

# (metric key in drug_mapping_diagnostics.json, human-readable label)
_COMPARISON_METRICS = [
    ("drug_id_mode", "drug_id_mode"),
    ("n_mapped_rows", "mapped row count"),
    ("mapping_rate", "mapping rate"),
    ("n_unique_mapped_drugs", "unique mapped drugs"),
    ("n_unique_chembl_ids_mapped", "unique ChEMBL ids mapped"),
    ("n_unique_cases_with_mapped_drug", "cases with >=1 mapped drug"),
    ("n_cases_with_ge2_mapped_drugs", "cases with >=2 mapped drugs"),
    ("n_cases_with_raw_ge2_mapped_lt2", "cases with raw>=2 but mapped<2"),
    ("n_removed_outside_drug_universe", "unmapped_outside_drug_universe count"),
    ("n_removed_ambiguous", "unmapped_ambiguous count"),
    ("n_removed_no_match", "unmapped_no_match count"),
    ("n_removed_fuzzy_excluded", "excluded_fuzzy count"),
]


def compare_drug_id_mode_diagnostics(diag_a: dict, diag_b: dict, label_a: str = "A", label_b: str = "B") -> None:
    """
    Print a side-by-side comparison of two drug_mapping_diagnostics.json
    payloads (Milestone 4c requirement 11), typically one produced with
    --drug-id-mode decagon_idx and one with --drug-id-mode chembl_id on the
    same --max-rows sample.

    This does not run the pipeline itself. Produce the two diagnostics
    files first, e.g.:
      python src/process_faers_case_level.py --max-rows 5000 --drug-id-mode decagon_idx --out-dir data/processed/faers_case_level_decagon
      python src/process_faers_case_level.py --max-rows 5000 --drug-id-mode chembl_id  --out-dir data/processed/faers_case_level_chembl
    then pass the two resulting drug_mapping_diagnostics.json paths via
    --compare-diagnostics (see main()).
    """
    print(f"\n  drug_id_mode comparison: {label_a} vs {label_b}")
    header = f"    {'metric':<38} {label_a:>18} {label_b:>18}"
    print(header)
    print("    " + "-" * (len(header) - 4))
    for key, label in _COMPARISON_METRICS:
        va, vb = diag_a.get(key, "n/a"), diag_b.get(key, "n/a")
        if key == "mapping_rate":
            va = f"{va:.1%}" if isinstance(va, (int, float)) else va
            vb = f"{vb:.1%}" if isinstance(vb, (int, float)) else vb
        print(f"    {label:<38} {str(va):>18} {str(vb):>18}")


def save_case_drugs_artifacts(
    out_dir: str,
    case_drugs_long_df: pd.DataFrame,
    case_drugs_df: pd.DataFrame,
    diagnostics: dict,
    unmapped_drugs_top_df: pd.DataFrame,
) -> None:
    long_path = os.path.join(out_dir, "case_drugs_long.parquet")
    case_drugs_long_df.to_parquet(long_path, index=False)
    print(f"  Saved: {long_path}  ({len(case_drugs_long_df):,} rows)")

    case_drugs_path = os.path.join(out_dir, "case_drugs.parquet")
    case_drugs_df.to_parquet(case_drugs_path, index=False)
    print(f"  Saved: {case_drugs_path}  ({len(case_drugs_df):,} rows)")

    diagnostics_path = os.path.join(out_dir, "drug_mapping_diagnostics.json")
    with open(diagnostics_path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    print(f"  Saved: {diagnostics_path}")

    unmapped_path = os.path.join(out_dir, "unmapped_drugs_top.parquet")
    unmapped_drugs_top_df.to_parquet(unmapped_path, index=False)
    print(f"  Saved: {unmapped_path}  ({len(unmapped_drugs_top_df):,} rows)")


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
            "once, ahead of drug-only R-GCN construction. Milestone 4c: "
            "builds case_drugs.parquet (one row per caseid, mapped_drug_id) "
            "with a selectable --drug-id-mode -- see file docstring for "
            "milestone status."
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
            "Include fuzzy-matched drug names from the mapping file when "
            "building case_drugs (Milestone 4b). By default these are "
            "excluded (mapping_status='excluded_fuzzy')."
        ),
    )
    parser.add_argument(
        "--drug-id-mode", choices=DRUG_ID_MODES, default=DEFAULT_DRUG_ID_MODE,
        help=(
            "How mapped_drug_id is derived (Milestone 4c). 'decagon_idx' "
            "(default, backward-compatible with Milestone 4b) restricts "
            "mapped drugs to the old Decagon drug_to_idx.json universe and "
            "uses its integer index as mapped_drug_id. 'chembl_id' "
            "(RECOMMENDED for the new drug-only R-GCN pipeline) uses the "
            "ChEMBL id string itself as mapped_drug_id, with no Decagon "
            "universe restriction."
        ),
    )
    parser.add_argument(
        "--max-rows", type=int, default=None,
        help="Limit DEMO/DRUG/REAC tables to the first N rows each. Smoke-testing only.",
    )
    parser.add_argument(
        "--compare-diagnostics", nargs=2, metavar=("DIAG_A_JSON", "DIAG_B_JSON"),
        default=None,
        help=(
            "Skip the pipeline entirely and just print a side-by-side "
            "comparison of two previously-saved drug_mapping_diagnostics.json "
            "files (Milestone 4c requirement 11), e.g. one from a "
            "--drug-id-mode decagon_idx run and one from a --drug-id-mode "
            "chembl_id run on the same --max-rows sample."
        ),
    )
    args = parser.parse_args()

    if args.compare_diagnostics:
        path_a, path_b = args.compare_diagnostics
        with open(path_a, encoding="utf-8") as f:
            diag_a = json.load(f)
        with open(path_b, encoding="utf-8") as f:
            diag_b = json.load(f)
        compare_drug_id_mode_diagnostics(
            diag_a, diag_b,
            label_a=diag_a.get("drug_id_mode", path_a),
            label_b=diag_b.get("drug_id_mode", path_b),
        )
        return

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
    print(f"  Drug ID mode:    {args.drug_id_mode}" + ("  (recommended for R-GCN)" if args.drug_id_mode == "chembl_id" else ""))
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

    print(f"\n--- Milestone 4c: build case_drugs (drug_id_mode={args.drug_id_mode}) ---")
    case_drugs_long_df, case_drugs_df, drug_mapping_diagnostics, unmapped_drugs_top_df = build_case_drugs(
        resolved_drug_df, args.mapping_file, drug_id_mode=args.drug_id_mode, include_fuzzy=args.include_fuzzy,
    )
    print_drug_mapping_diagnostics(drug_mapping_diagnostics)

    print("\n--- Milestone 4c: case_drugs validation checks ---")
    run_case_drugs_validation(resolved_drug_df, case_drugs_long_df, case_drugs_df, args.drug_id_mode)

    print("\n--- Milestone 4c: saving case_drugs artifacts ---")
    save_case_drugs_artifacts(
        args.out_dir, case_drugs_long_df, case_drugs_df, drug_mapping_diagnostics, unmapped_drugs_top_df,
    )

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
    print("Milestone 4c complete: case_drugs.parquet built (one row per")
    print(f"caseid, mapped_drug_id) using drug_id_mode={args.drug_id_mode!r}.")
    print("No drug-pair generation or PRR was performed yet.")
    print(f"Artifacts in: {args.out_dir}")
    print("Next: implement Milestone 5 (canonical drug-pair generation).")
    print("=" * 60)


if __name__ == "__main__":
    main()
