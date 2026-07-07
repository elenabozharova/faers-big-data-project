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
  [x] Milestone 4c -- selectable drug ID mode (decagon_idx / chembl_id)
  [x] Milestone 5 -- canonical (drug_a, drug_b) pair generation per case
  [x] Milestone 5b -- expected pair-reaction join-size audit (no join)
  [x] Milestone 5c -- role-aware pair audit (PS/SS/C/I, no filtering)
  [x] Milestone 5d -- configurable pair-role filter (case_drug_pairs_selected)
  [x] Milestone 6 -- case_pair_reactions + drug_pair_reaction_counts (no PRR)  (this file)
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

Milestone 5 scope
------------------
Builds case_drug_pairs.parquet: per case, every unordered (drug_a, drug_b)
combination of that case's unique mapped_drug_id values from case_drugs.parquet
(drug_a < drug_b by string order, no self-pairs, no reversed duplicates).
Cases with fewer than 2 mapped drugs emit no pairs. High-polypharmacy cases
(many mapped drugs -> combinatorially many pairs) are audited via
drug_pair_generation_diagnostics.json and high_polypharmacy_cases.parquet,
not capped or filtered -- that decision is deferred, not made here.

Milestone 5b scope
--------------------
Builds expected_pair_reaction_join_audit.parquet: per case, the *expected*
row count of the eventual case_pair_reactions join (n_pairs_for_case *
n_reactions_for_case), WITHOUT performing that join. Exists to answer "how
big will Milestone 6's output be, and is it dominated by a handful of
high-polypharmacy cases" before paying the cost of materializing it.

Milestone 5c scope
-------------------
Builds case_drug_pair_role_audit.parquet: one row per case_drug_pairs row,
with role-aware metadata (PS/SS/C/I per drug in the pair, C-C flag,
has-suspect-or-interacting flag) joined in from case_drugs.roles_seen, plus
drug_pair_role_audit_diagnostics.json, a sensitivity-analysis-only report
of how many pairs/cases would remain under an "exclude C-C" or "require
>=1 PS/SS/I" rule. Deliberately a separate table rather than new columns
on case_drug_pairs.parquet, so Milestone 5's output and validation are
untouched. No filtering is applied to case_drug_pairs.parquet itself --
concomitant (C) drugs can still be clinically relevant, so this is
reported as a what-if, not treated as automatically correct.

Milestone 5d scope
-------------------
Adds --pair-role-filter {all_pairs, non_cc} and builds
case_drug_pairs_selected.parquet: the actual set of pairs that proceed to
the pair-reaction join under the chosen filter, computed from
case_drug_pairs.parquet + case_drug_pair_role_audit without ever mutating
either. 'all_pairs' (default) keeps everything; 'non_cc' (recommended for
the thesis comparison) drops C-C pairs.

Milestone 6 scope (this file)
-------------------------------
Builds case_pair_reactions.parquet (caseid's selected pair x that case's
reactions, i.e. CO-REPORTING only -- the pair and the reaction were both
reported in the same FAERS case; this is NOT a causal or mechanistic
claim) and drug_pair_reaction_counts.parquet (n_cases = distinct caseid
count per (drug_a, drug_b, reaction_pt), never a raw row count), plus
min-support-filtered variants (n_cases >= 3/5/10). The join can be tens of
millions of rows (measured ~41M for --max-rows 300000 --pair-role-filter
non_cc), so it is computed and written in bounded-size chunks
(build_case_pair_reactions_chunked) rather than as one in-memory pandas
merge -- see that function's docstring for why this is still exact (not
an approximation) despite never materializing the full table at once.
Still explicitly out of scope:
  - PRR / disproportionality statistics (Milestone 7)
  - R-GCN export, train/validation/test splitting
  - collapsing ChEMBL salt/form ids to a parent active ingredient

This file does NOT:
  - implement full preprocessing
  - claim or compute causality between a drug pair and a reaction
  - compute PRR (Proportional Reporting Ratio)
  - touch R-GCN, ChEMBL, STRING, or Decagon adjacency / training code

Usage
-----
  python src/process_faers_case_level.py
  python src/process_faers_case_level.py --quarter 2025Q4 --max-rows 50000
  python src/process_faers_case_level.py --faers-dir data/raw/faers/extracted_2025Q4
"""

import argparse
import itertools
import json
import os
import re
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

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
# Milestone 5 -- canonical drug-pair generation per case (case_drug_pairs.parquet)
# ---------------------------------------------------------------------------

# Audit-only thresholds for "high polypharmacy" cases (Milestone 5 requirement
# 9): cases at/above these mapped-drug counts are reported in diagnostics and
# saved to high_polypharmacy_cases.parquet, but never capped or filtered out.
_HIGH_POLYPHARMACY_THRESHOLD = 10
_HIGH_POLYPHARMACY_THRESHOLD_STRICT = 15


def build_case_drug_pairs(case_drugs_df: pd.DataFrame, drug_id_mode: str) -> tuple:
    """
    Build canonical (drug_a, drug_b) pairs per case from case_drugs_df (one
    row per (caseid, mapped_drug_id), mapped rows only -- Milestone 4b/4c).

    Per case: take the case's unique mapped_drug_id values, sort them by
    string value for a deterministic canonical order (works uniformly for
    ChEMBL id strings in chembl_id mode and for Decagon integer indices in
    decagon_idx mode -- the sort only needs to be stable, not numeric), then
    emit every unordered combination with drug_a < drug_b. A case with fewer
    than 2 mapped drugs emits no pairs. No joining with reactions and no
    capping/filtering of high-polypharmacy cases happens here.

    Returns (case_drug_pairs_df, diagnostics_dict, high_polypharmacy_df).
    """
    _require_columns(case_drugs_df, ["caseid", "mapped_drug_id"], "case_drugs")

    n_input_rows = len(case_drugs_df)
    n_unique_cases = int(case_drugs_df["caseid"].nunique())

    drugs_per_case = (
        case_drugs_df[["caseid", "mapped_drug_id"]]
        .drop_duplicates()
        .assign(mapped_drug_id=lambda d: d["mapped_drug_id"].astype(str))
        .groupby("caseid")["mapped_drug_id"]
        .apply(lambda s: sorted(s))
    )
    n_mapped_per_case = drugs_per_case.apply(len).rename("n_mapped_drugs_in_case")

    print(f"  Generating pairs for {n_unique_cases:,} cases ({int((n_mapped_per_case >= 2).sum()):,} with >=2 mapped drugs)...")

    pair_rows = [
        (caseid, drug_a, drug_b, len(drugs))
        for caseid, drugs in drugs_per_case.items()
        for drug_a, drug_b in itertools.combinations(drugs, 2)
    ]
    case_drug_pairs_df = pd.DataFrame(
        pair_rows, columns=["caseid", "drug_a", "drug_b", "n_mapped_drugs_in_case"]
    )
    case_drug_pairs_df["drug_id_mode"] = drug_id_mode

    n_pairs = len(case_drug_pairs_df)
    n_unique_pairs = (
        int(case_drug_pairs_df[["drug_a", "drug_b"]].drop_duplicates().shape[0]) if n_pairs else 0
    )
    pairs_per_case = (
        case_drug_pairs_df.groupby("caseid").size() if n_pairs else pd.Series(dtype=int)
    )

    print(
        f"  case_drug_pairs: {n_pairs:,} pair rows "
        f"({n_unique_pairs:,} unique (drug_a, drug_b) pairs across all cases)"
    )

    n_cases_exactly_1 = int((n_mapped_per_case == 1).sum())
    n_cases_ge2 = int((n_mapped_per_case >= 2).sum())
    n_high_10 = int((n_mapped_per_case >= _HIGH_POLYPHARMACY_THRESHOLD).sum())
    n_high_15 = int((n_mapped_per_case >= _HIGH_POLYPHARMACY_THRESHOLD_STRICT).sum())

    dist = n_mapped_per_case.value_counts().sort_index()

    top_20 = n_mapped_per_case.sort_values(ascending=False).head(20).reset_index()
    top_20_cases = [
        {"caseid": str(row["caseid"]), "n_mapped_drugs_in_case": int(row["n_mapped_drugs_in_case"])}
        for _, row in top_20.iterrows()
    ]

    diagnostics = {
        "drug_id_mode": drug_id_mode,
        "n_input_case_drugs_rows": n_input_rows,
        "n_unique_cases_in_case_drugs": n_unique_cases,
        "n_cases_with_exactly_1_mapped_drug": n_cases_exactly_1,
        "n_cases_with_ge2_mapped_drugs": n_cases_ge2,
        "n_generated_pair_rows": n_pairs,
        "n_unique_drug_pairs": n_unique_pairs,
        "max_mapped_drugs_in_one_case": int(n_mapped_per_case.max()) if len(n_mapped_per_case) else 0,
        "max_generated_pairs_in_one_case": int(pairs_per_case.max()) if len(pairs_per_case) else 0,
        "n_mapped_drugs_in_case_distribution": {str(int(k)): int(v) for k, v in dist.items()},
        f"n_high_polypharmacy_cases_ge{_HIGH_POLYPHARMACY_THRESHOLD}": n_high_10,
        f"n_high_polypharmacy_cases_ge{_HIGH_POLYPHARMACY_THRESHOLD_STRICT}": n_high_15,
        "top_20_cases_by_mapped_drug_count": top_20_cases,
    }

    # ---- high_polypharmacy_cases: audit only, never capped/filtered above ----
    high_mask = n_mapped_per_case >= _HIGH_POLYPHARMACY_THRESHOLD
    high_polypharmacy_df = n_mapped_per_case[high_mask].reset_index()
    if len(high_polypharmacy_df):
        high_polypharmacy_df["n_pairs_generated"] = (
            high_polypharmacy_df["n_mapped_drugs_in_case"]
            * (high_polypharmacy_df["n_mapped_drugs_in_case"] - 1) // 2
        )
    else:
        high_polypharmacy_df["n_pairs_generated"] = pd.Series(dtype=int)
    high_polypharmacy_df["drug_id_mode"] = drug_id_mode
    high_polypharmacy_df = high_polypharmacy_df.sort_values(
        "n_mapped_drugs_in_case", ascending=False
    ).reset_index(drop=True)

    return case_drug_pairs_df, diagnostics, high_polypharmacy_df


def print_drug_pair_diagnostics(diagnostics: dict) -> None:
    d = diagnostics
    print("\n  Drug pair generation diagnostics:")
    print(f"    drug_id_mode:                          {d['drug_id_mode']}")
    print(f"    input case_drugs rows:                 {d['n_input_case_drugs_rows']:,}")
    print(f"    unique cases in case_drugs:             {d['n_unique_cases_in_case_drugs']:,}")
    print(f"    cases with exactly 1 mapped drug:       {d['n_cases_with_exactly_1_mapped_drug']:,}")
    print(f"    cases with >=2 mapped drugs:            {d['n_cases_with_ge2_mapped_drugs']:,}")
    print(f"    generated pair rows:                    {d['n_generated_pair_rows']:,}")
    print(f"    unique (drug_a, drug_b) pairs:          {d['n_unique_drug_pairs']:,}")
    print(f"    max mapped drugs in one case:           {d['max_mapped_drugs_in_one_case']:,}")
    print(f"    max generated pairs in one case:        {d['max_generated_pairs_in_one_case']:,}")
    print(f"    high-polypharmacy cases (>={_HIGH_POLYPHARMACY_THRESHOLD}):          {d[f'n_high_polypharmacy_cases_ge{_HIGH_POLYPHARMACY_THRESHOLD}']:,}")
    print(f"    high-polypharmacy cases (>={_HIGH_POLYPHARMACY_THRESHOLD_STRICT}):          {d[f'n_high_polypharmacy_cases_ge{_HIGH_POLYPHARMACY_THRESHOLD_STRICT}']:,}")
    print("    n_mapped_drugs_in_case distribution (top 20 bucket values shown):")
    dist_items = sorted(d["n_mapped_drugs_in_case_distribution"].items(), key=lambda kv: int(kv[0]))
    for k, v in dist_items[:20]:
        print(f"      {k:>4} drug(s): {v:,} cases")
    if len(dist_items) > 20:
        print(f"      ... ({len(dist_items) - 20} more distinct counts)")
    print("    top 20 cases by mapped drug count:")
    for row in d["top_20_cases_by_mapped_drug_count"]:
        print(f"      caseid={row['caseid']:<12} n_mapped_drugs_in_case={row['n_mapped_drugs_in_case']:,}")


def run_case_drug_pairs_validation(
    case_drugs_df: pd.DataFrame,
    case_drug_pairs_df: pd.DataFrame,
    drug_id_mode: str,
) -> None:
    """
    Validate case_drug_pairs. Raises AssertionError listing every failed
    check -- a pair-generation bug here would silently corrupt every
    downstream pair-reaction count and PRR calculation, so it fails loudly.
    """
    failures = []

    required_cols = ["caseid", "drug_a", "drug_b", "n_mapped_drugs_in_case", "drug_id_mode"]
    missing_cols = [c for c in required_cols if c not in case_drug_pairs_df.columns]
    if missing_cols:
        failures.append(f"case_drug_pairs is missing required column(s): {missing_cols}")
        raise AssertionError("Milestone 5 (case_drug_pairs) validation failed:\n  - " + "\n  - ".join(failures))

    if case_drug_pairs_df["caseid"].isna().any():
        failures.append("case_drug_pairs has rows with null caseid")
    if case_drug_pairs_df["drug_a"].isna().any():
        failures.append("case_drug_pairs has rows with null drug_a")
    if case_drug_pairs_df["drug_b"].isna().any():
        failures.append("case_drug_pairs has rows with null drug_b")

    if (case_drug_pairs_df["drug_a"] == case_drug_pairs_df["drug_b"]).any():
        failures.append("case_drug_pairs has self-pairs (drug_a == drug_b)")

    if (case_drug_pairs_df["drug_a"] > case_drug_pairs_df["drug_b"]).any():
        failures.append("case_drug_pairs has rows where drug_a > drug_b (canonical ordering violated)")

    if case_drug_pairs_df.duplicated(subset=["caseid", "drug_a", "drug_b"]).any():
        failures.append("case_drug_pairs has duplicate (caseid, drug_a, drug_b) rows")

    # Generated pair count per case must equal n * (n - 1) / 2, n = unique
    # mapped drugs for that case in case_drugs (not just cases that produced
    # pairs -- a case with 1 mapped drug must produce exactly 0).
    n_unique_drugs_per_case = (
        case_drugs_df[["caseid", "mapped_drug_id"]].drop_duplicates().groupby("caseid").size()
    )
    expected_pairs_per_case = (n_unique_drugs_per_case * (n_unique_drugs_per_case - 1) // 2).rename(
        "expected_pairs"
    )
    actual_pairs_per_case = case_drug_pairs_df.groupby("caseid").size().rename("actual_pairs")
    check_df = pd.DataFrame({"expected_pairs": expected_pairs_per_case}).join(
        pd.DataFrame({"actual_pairs": actual_pairs_per_case}), how="left"
    )
    check_df["actual_pairs"] = check_df["actual_pairs"].fillna(0).astype(int)
    mismatched = check_df["expected_pairs"] != check_df["actual_pairs"]
    if mismatched.any():
        failures.append(
            f"{int(mismatched.sum()):,} cases have generated pair count != n*(n-1)/2 "
            f"for their unique mapped-drug count in case_drugs"
        )

    valid_caseids = set(case_drugs_df["caseid"])
    if not set(case_drug_pairs_df["caseid"]).issubset(valid_caseids):
        failures.append("case_drug_pairs contains a caseid not present in case_drugs")

    modes = set(case_drug_pairs_df["drug_id_mode"].dropna().unique())
    if modes and modes != {drug_id_mode}:
        failures.append(
            f"case_drug_pairs.drug_id_mode has value(s) {sorted(modes)}, expected only {drug_id_mode!r}"
        )
    elif not modes and len(case_drug_pairs_df):
        failures.append("case_drug_pairs.drug_id_mode is entirely null")

    if failures:
        raise AssertionError("Milestone 5 (case_drug_pairs) validation failed:\n  - " + "\n  - ".join(failures))

    print("  All Milestone 5 (case_drug_pairs) validation checks passed:")
    print("    - case_drug_pairs has required columns")
    print("    - no missing caseid / drug_a / drug_b")
    print("    - no self-pairs (drug_a == drug_b)")
    print("    - no rows with drug_a > drug_b")
    print("    - no duplicate (caseid, drug_a, drug_b)")
    print("    - generated pair count == n*(n-1)/2 for every case")
    print("    - every case_drug_pairs.caseid present in case_drugs")
    print("    - drug_id_mode present and consistent with the current run")


def save_case_drug_pairs_artifacts(
    out_dir: str,
    case_drug_pairs_df: pd.DataFrame,
    diagnostics: dict,
    high_polypharmacy_df: pd.DataFrame,
) -> None:
    pairs_path = os.path.join(out_dir, "case_drug_pairs.parquet")
    case_drug_pairs_df.to_parquet(pairs_path, index=False)
    print(f"  Saved: {pairs_path}  ({len(case_drug_pairs_df):,} rows)")

    diagnostics_path = os.path.join(out_dir, "drug_pair_generation_diagnostics.json")
    with open(diagnostics_path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    print(f"  Saved: {diagnostics_path}")

    high_path = os.path.join(out_dir, "high_polypharmacy_cases.parquet")
    high_polypharmacy_df.to_parquet(high_path, index=False)
    print(f"  Saved: {high_path}  ({len(high_polypharmacy_df):,} rows, >= {_HIGH_POLYPHARMACY_THRESHOLD} mapped drugs; audit only, not filtered)")


# ---------------------------------------------------------------------------
# Milestone 5b -- expected pair-reaction join-size audit (audit only, no join)
# ---------------------------------------------------------------------------
#
# TODO (deferred, not addressed by this milestone):
#   - high-polypharmacy cases may dominate the eventual case_pair_reactions
#     table (this audit exists to quantify that before we pay the cost of
#     materializing the join)
#   - a future sensitivity analysis may compare capping thresholds such as
#     max 50 or max 100 mapped drugs per case and re-measure join size
#   - a future role-aware analysis may compare "all pairs" against pairs
#     filtered/weighted by role_cod (PS/SS/I/C), since this milestone's
#     n_pairs_for_case counts every mapped-drug pair regardless of role
#   - a future mapping refinement may collapse ChEMBL salt/form ids (e.g.
#     "methotrexate" vs "methotrexate sodium") to a parent active
#     ingredient, which would shrink n_mapped_drugs_in_case (and therefore
#     n_pairs_for_case) for cases like caseid 14638057 (see manual audit)

# Spot-check values for the manually-audited high-polypharmacy case from the
# --max-rows 300000 --drug-id-mode chembl_id run (200 mapped drugs -> 19,900
# pairs; 196 reactions). Used by run_expected_pair_reaction_join_audit_validation()
# as a regression check, not a general invariant -- skipped if this caseid
# isn't present in the current run (e.g. a smaller --max-rows sample).
_KNOWN_AUDIT_CASEID = "14638057"
_KNOWN_AUDIT_CASEID_N_PAIRS = 19900
_KNOWN_AUDIT_CASEID_N_REACTIONS = 196


def build_expected_pair_reaction_join_audit(
    case_drug_pairs_df: pd.DataFrame,
    case_reactions_df: pd.DataFrame,
) -> tuple:
    """
    Audit-only: estimate the size of the eventual case_pair_reactions join
    (Milestone 6) before materializing it, and quantify whether a small
    number of high-polypharmacy cases dominate that size. Does NOT perform
    the actual pair/reaction join -- no drug or reaction row is combined
    here, only per-case counts and their product.

    Per caseid (union of caseids appearing in case_drug_pairs_df and/or
    case_reactions_df):
      n_pairs_for_case          -- row count for that caseid in case_drug_pairs_df (0 if none)
      n_reactions_for_case      -- row count for that caseid in case_reactions_df (0 if none)
      n_mapped_drugs_in_case    -- carried from case_drug_pairs_df (constant per case there);
                                    0 for a case with no pairs, since a case with 0 or 1 mapped
                                    drug never appears in case_drug_pairs_df and this function
                                    (by design) does not take case_drugs_df as input to
                                    distinguish those two. This does not affect
                                    expected_pair_reaction_rows, which is correctly 0 for such
                                    cases regardless (n_pairs_for_case is already 0).
      expected_pair_reaction_rows = n_pairs_for_case * n_reactions_for_case

    Returns (audit_df, diagnostics_dict).
    """
    _require_columns(case_drug_pairs_df, ["caseid", "n_mapped_drugs_in_case"], "case_drug_pairs")
    _require_columns(case_reactions_df, ["caseid"], "case_reactions")

    if len(case_drug_pairs_df):
        pairs_per_case = case_drug_pairs_df.groupby("caseid").agg(
            n_pairs_for_case=("caseid", "size"),
            n_mapped_drugs_in_case=("n_mapped_drugs_in_case", "first"),
        )
    else:
        pairs_per_case = pd.DataFrame(
            columns=["n_pairs_for_case", "n_mapped_drugs_in_case"]
        ).rename_axis("caseid")

    if len(case_reactions_df):
        reactions_per_case = (
            case_reactions_df.groupby("caseid").size().rename("n_reactions_for_case").to_frame()
        )
    else:
        reactions_per_case = pd.DataFrame(columns=["n_reactions_for_case"]).rename_axis("caseid")

    audit_df = pairs_per_case.join(reactions_per_case, how="outer")
    audit_df["n_pairs_for_case"] = audit_df["n_pairs_for_case"].fillna(0).astype(int)
    audit_df["n_reactions_for_case"] = audit_df["n_reactions_for_case"].fillna(0).astype(int)
    audit_df["n_mapped_drugs_in_case"] = audit_df["n_mapped_drugs_in_case"].fillna(0).astype(int)
    audit_df["expected_pair_reaction_rows"] = (
        audit_df["n_pairs_for_case"].astype("int64") * audit_df["n_reactions_for_case"].astype("int64")
    )
    audit_df = audit_df.reset_index().rename(columns={"index": "caseid"})[
        ["caseid", "n_mapped_drugs_in_case", "n_pairs_for_case", "n_reactions_for_case", "expected_pair_reaction_rows"]
    ]

    n_cases_with_pairs = int((audit_df["n_pairs_for_case"] > 0).sum())
    n_cases_with_reactions = int((audit_df["n_reactions_for_case"] > 0).sum())
    n_cases_with_both = int(
        ((audit_df["n_pairs_for_case"] > 0) & (audit_df["n_reactions_for_case"] > 0)).sum()
    )

    total_expected = int(audit_df["expected_pair_reaction_rows"].sum())
    print(
        f"  Expected pair-reaction join size: {total_expected:,} rows across "
        f"{len(audit_df):,} cases ({n_cases_with_both:,} with both pairs and reactions)"
    )

    sorted_by_expected = audit_df.sort_values("expected_pair_reaction_rows", ascending=False)
    top_20 = sorted_by_expected.head(20)
    top_20_cases = [
        {
            "caseid": str(row["caseid"]),
            "n_mapped_drugs_in_case": int(row["n_mapped_drugs_in_case"]),
            "n_pairs_for_case": int(row["n_pairs_for_case"]),
            "n_reactions_for_case": int(row["n_reactions_for_case"]),
            "expected_pair_reaction_rows": int(row["expected_pair_reaction_rows"]),
        }
        for _, row in top_20.iterrows()
    ]

    def _pct_from_mask(mask: pd.Series) -> float:
        return (float(audit_df.loc[mask, "expected_pair_reaction_rows"].sum()) / total_expected * 100.0) if total_expected else 0.0

    def _pct_top_n(n: int) -> float:
        return (float(sorted_by_expected["expected_pair_reaction_rows"].head(n).sum()) / total_expected * 100.0) if total_expected else 0.0

    diagnostics = {
        "total_expected_pair_reaction_rows": total_expected,
        "n_cases_in_audit": int(len(audit_df)),
        "n_cases_with_pairs": n_cases_with_pairs,
        "n_cases_with_reactions": n_cases_with_reactions,
        "n_cases_with_both_pairs_and_reactions": n_cases_with_both,
        "max_mapped_drugs_in_one_case": int(audit_df["n_mapped_drugs_in_case"].max()) if len(audit_df) else 0,
        "max_pairs_in_one_case": int(audit_df["n_pairs_for_case"].max()) if len(audit_df) else 0,
        "max_reactions_in_one_case": int(audit_df["n_reactions_for_case"].max()) if len(audit_df) else 0,
        "max_expected_rows_from_one_case": int(audit_df["expected_pair_reaction_rows"].max()) if len(audit_df) else 0,
        "pct_expected_rows_from_cases_with_mapped_drugs_ge10": _pct_from_mask(audit_df["n_mapped_drugs_in_case"] >= 10),
        "pct_expected_rows_from_cases_with_mapped_drugs_ge20": _pct_from_mask(audit_df["n_mapped_drugs_in_case"] >= 20),
        "pct_expected_rows_from_cases_with_mapped_drugs_ge50": _pct_from_mask(audit_df["n_mapped_drugs_in_case"] >= 50),
        "pct_expected_rows_from_cases_with_mapped_drugs_ge100": _pct_from_mask(audit_df["n_mapped_drugs_in_case"] >= 100),
        "pct_expected_rows_from_top_1_case": _pct_top_n(1),
        "pct_expected_rows_from_top_10_cases": _pct_top_n(10),
        "pct_expected_rows_from_top_100_cases": _pct_top_n(100),
        "top_20_cases_by_expected_pair_reaction_rows": top_20_cases,
    }

    return audit_df, diagnostics


def print_expected_pair_reaction_join_diagnostics(diagnostics: dict) -> None:
    d = diagnostics
    print("\n  Expected pair-reaction join-size audit:")
    print(f"    total expected pair-reaction rows:     {d['total_expected_pair_reaction_rows']:,}")
    print(f"    cases in audit:                        {d['n_cases_in_audit']:,}")
    print(f"    cases with pairs:                      {d['n_cases_with_pairs']:,}")
    print(f"    cases with reactions:                  {d['n_cases_with_reactions']:,}")
    print(f"    cases with both pairs and reactions:   {d['n_cases_with_both_pairs_and_reactions']:,}")
    print(f"    max mapped drugs in one case:           {d['max_mapped_drugs_in_one_case']:,}")
    print(f"    max pairs in one case:                  {d['max_pairs_in_one_case']:,}")
    print(f"    max reactions in one case:               {d['max_reactions_in_one_case']:,}")
    print(f"    max expected rows from one case:        {d['max_expected_rows_from_one_case']:,}")
    print("    % of expected rows from cases with n_mapped_drugs_in_case >=:")
    print(f"      >=10:  {d['pct_expected_rows_from_cases_with_mapped_drugs_ge10']:.2f}%")
    print(f"      >=20:  {d['pct_expected_rows_from_cases_with_mapped_drugs_ge20']:.2f}%")
    print(f"      >=50:  {d['pct_expected_rows_from_cases_with_mapped_drugs_ge50']:.2f}%")
    print(f"      >=100: {d['pct_expected_rows_from_cases_with_mapped_drugs_ge100']:.2f}%")
    print("    % of expected rows from:")
    print(f"      top 1 case:    {d['pct_expected_rows_from_top_1_case']:.2f}%")
    print(f"      top 10 cases:  {d['pct_expected_rows_from_top_10_cases']:.2f}%")
    print(f"      top 100 cases: {d['pct_expected_rows_from_top_100_cases']:.2f}%")
    print("    top 20 cases by expected pair-reaction rows:")
    for row in d["top_20_cases_by_expected_pair_reaction_rows"]:
        print(
            f"      caseid={row['caseid']:<12} "
            f"n_mapped_drugs={row['n_mapped_drugs_in_case']:>4} "
            f"n_pairs={row['n_pairs_for_case']:>7,} "
            f"n_reactions={row['n_reactions_for_case']:>4} "
            f"expected_rows={row['expected_pair_reaction_rows']:>10,}"
        )


def run_expected_pair_reaction_join_audit_validation(
    case_drug_pairs_df: pd.DataFrame,
    case_reactions_df: pd.DataFrame,
    audit_df: pd.DataFrame,
) -> None:
    """
    Validate expected_pair_reaction_join_audit. Raises AssertionError listing
    every failed check -- this audit exists specifically to catch a
    dominating-case problem before Milestone 6 pays the cost of the real
    join, so a silent bug here would defeat its purpose.
    """
    failures = []

    required_cols = [
        "caseid", "n_mapped_drugs_in_case", "n_pairs_for_case",
        "n_reactions_for_case", "expected_pair_reaction_rows",
    ]
    missing_cols = [c for c in required_cols if c not in audit_df.columns]
    if missing_cols:
        failures.append(f"audit_df is missing required column(s): {missing_cols}")
        raise AssertionError(
            "Expected pair-reaction join audit validation failed:\n  - " + "\n  - ".join(failures)
        )

    if audit_df["caseid"].isna().any():
        failures.append("audit_df has rows with null caseid")
    if (audit_df["n_pairs_for_case"] < 0).any():
        failures.append("audit_df has rows with n_pairs_for_case < 0")
    if (audit_df["n_reactions_for_case"] < 0).any():
        failures.append("audit_df has rows with n_reactions_for_case < 0")

    recomputed = audit_df["n_pairs_for_case"].astype("int64") * audit_df["n_reactions_for_case"].astype("int64")
    if (audit_df["expected_pair_reaction_rows"].astype("int64") != recomputed).any():
        failures.append("audit_df.expected_pair_reaction_rows != n_pairs_for_case * n_reactions_for_case for some rows")

    valid_caseids = set(case_drug_pairs_df["caseid"]) | set(case_reactions_df["caseid"])
    if not set(audit_df["caseid"]).issubset(valid_caseids):
        failures.append("audit_df contains a caseid not present in case_drug_pairs or case_reactions")

    known_case_rows = audit_df[audit_df["caseid"] == _KNOWN_AUDIT_CASEID]
    if len(known_case_rows):
        row = known_case_rows.iloc[0]
        if int(row["n_pairs_for_case"]) != _KNOWN_AUDIT_CASEID_N_PAIRS:
            failures.append(
                f"caseid {_KNOWN_AUDIT_CASEID} (known audit case) has n_pairs_for_case="
                f"{int(row['n_pairs_for_case']):,}, expected {_KNOWN_AUDIT_CASEID_N_PAIRS:,}"
            )
        if int(row["n_reactions_for_case"]) != _KNOWN_AUDIT_CASEID_N_REACTIONS:
            failures.append(
                f"caseid {_KNOWN_AUDIT_CASEID} (known audit case) has n_reactions_for_case="
                f"{int(row['n_reactions_for_case']):,}, expected {_KNOWN_AUDIT_CASEID_N_REACTIONS:,}"
            )

    if failures:
        raise AssertionError(
            "Expected pair-reaction join audit validation failed:\n  - " + "\n  - ".join(failures)
        )

    print("  All expected pair-reaction join audit validation checks passed:")
    print("    - audit_df has required columns")
    print("    - no missing caseid")
    print("    - n_pairs_for_case >= 0 and n_reactions_for_case >= 0")
    print("    - expected_pair_reaction_rows == n_pairs_for_case * n_reactions_for_case")
    print("    - every audit caseid present in case_drug_pairs or case_reactions")
    print(f"    - known audit case {_KNOWN_AUDIT_CASEID} matches expected pairs/reactions (if present in this run)")


def save_expected_pair_reaction_join_audit_artifacts(
    out_dir: str,
    audit_df: pd.DataFrame,
    diagnostics: dict,
) -> None:
    audit_path = os.path.join(out_dir, "expected_pair_reaction_join_audit.parquet")
    audit_df.to_parquet(audit_path, index=False)
    print(f"  Saved: {audit_path}  ({len(audit_df):,} rows)")

    diagnostics_path = os.path.join(out_dir, "expected_pair_reaction_join_diagnostics.json")
    with open(diagnostics_path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    print(f"  Saved: {diagnostics_path}")


# ---------------------------------------------------------------------------
# Milestone 5c -- role-aware pair metadata / sensitivity audit (no filtering)
# ---------------------------------------------------------------------------
#
# TODO (deferred, not addressed by this milestone):
#   - a future final-preprocessing decision may compare three variants:
#     all pairs / exclude C-C / require >=1 PS, SS, or I -- this milestone
#     only reports what each variant's size would be, it picks none of them
#   - role-aware filtering (dropping C-C pairs, or requiring a suspect role)
#     is NOT automatically "more correct": concomitant (C) drugs can still
#     be clinically relevant to an adverse event (e.g. a drug-drug
#     interaction where neither drug is the "suspect"), so this is a
#     sensitivity analysis to inform that decision, not a recommendation

_ROLE_PS, _ROLE_SS, _ROLE_I, _ROLE_C = "PS", "SS", "I", "C"


def _role_set(roles) -> set:
    """Normalize a case_drugs.roles_seen cell (list/tuple/ndarray/None) to a plain set."""
    if roles is None or (isinstance(roles, float) and pd.isna(roles)):
        return set()
    try:
        return set(roles)
    except TypeError:
        return set()


def build_case_drug_pair_role_audit(
    case_drug_pairs_df: pd.DataFrame,
    case_drugs_df: pd.DataFrame,
    drug_id_mode: str,
) -> tuple:
    """
    Sensitivity-analysis-only audit: attach FAERS role_cod metadata
    (PS/SS/C/I, from case_drugs.roles_seen) to each row of case_drug_pairs
    and report how many pairs/cases would remain under an "exclude C-C" or
    "require >=1 PS/SS/I" rule.

    Deliberately a SEPARATE table (case_drug_pair_role_audit.parquet)
    rather than new columns bolted onto case_drug_pairs.parquet: keeps
    Milestone 5's case_drug_pairs.parquet and its validation completely
    unchanged, and keeps this audit's role-metadata columns out of the
    core pair table that Milestone 6 will actually join against. Does NOT
    filter case_drug_pairs.parquet -- see the module-level TODO above.

    Returns (role_audit_df, diagnostics_dict).
    """
    _require_columns(case_drug_pairs_df, ["caseid", "drug_a", "drug_b"], "case_drug_pairs")
    _require_columns(case_drugs_df, ["caseid", "mapped_drug_id", "roles_seen"], "case_drugs")

    roles_only = case_drugs_df[["caseid", "mapped_drug_id", "roles_seen"]].copy()
    roles_only["mapped_drug_id"] = roles_only["mapped_drug_id"].astype(str)

    audit_df = case_drug_pairs_df[["caseid", "drug_a", "drug_b"]].copy()
    audit_df["drug_a"] = audit_df["drug_a"].astype(str)
    audit_df["drug_b"] = audit_df["drug_b"].astype(str)

    audit_df = audit_df.merge(
        roles_only.rename(columns={"mapped_drug_id": "drug_a", "roles_seen": "roles_a"}),
        on=["caseid", "drug_a"], how="left",
    ).merge(
        roles_only.rename(columns={"mapped_drug_id": "drug_b", "roles_seen": "roles_b"}),
        on=["caseid", "drug_b"], how="left",
    )

    role_sets_a = audit_df["roles_a"].apply(_role_set)
    role_sets_b = audit_df["roles_b"].apply(_role_set)

    audit_df["drug_a_has_ps"] = role_sets_a.apply(lambda s: _ROLE_PS in s)
    audit_df["drug_b_has_ps"] = role_sets_b.apply(lambda s: _ROLE_PS in s)
    audit_df["drug_a_has_ss"] = role_sets_a.apply(lambda s: _ROLE_SS in s)
    audit_df["drug_b_has_ss"] = role_sets_b.apply(lambda s: _ROLE_SS in s)
    audit_df["drug_a_has_i"] = role_sets_a.apply(lambda s: _ROLE_I in s)
    audit_df["drug_b_has_i"] = role_sets_b.apply(lambda s: _ROLE_I in s)
    audit_df["drug_a_only_c"] = role_sets_a.apply(lambda s: s == {_ROLE_C})
    audit_df["drug_b_only_c"] = role_sets_b.apply(lambda s: s == {_ROLE_C})
    audit_df["is_cc_pair"] = audit_df["drug_a_only_c"] & audit_df["drug_b_only_c"]

    pair_has_ps = audit_df["drug_a_has_ps"] | audit_df["drug_b_has_ps"]
    pair_has_ss = audit_df["drug_a_has_ss"] | audit_df["drug_b_has_ss"]
    pair_has_i = audit_df["drug_a_has_i"] | audit_df["drug_b_has_i"]
    audit_df["has_suspect_or_interacting"] = pair_has_ps | pair_has_ss | pair_has_i
    audit_df["drug_id_mode"] = drug_id_mode

    audit_df = audit_df[[
        "caseid", "drug_a", "drug_b", "roles_a", "roles_b",
        "drug_a_has_ps", "drug_b_has_ps", "drug_a_has_ss", "drug_b_has_ss",
        "drug_a_has_i", "drug_b_has_i", "drug_a_only_c", "drug_b_only_c",
        "is_cc_pair", "has_suspect_or_interacting", "drug_id_mode",
    ]]

    n_total = len(audit_df)
    n_cc = int(audit_df["is_cc_pair"].sum())
    n_with_ps = int(pair_has_ps.sum())
    n_with_ss = int(pair_has_ss.sum())
    n_with_i = int(pair_has_i.sum())
    n_suspect_or_interacting = int(audit_df["has_suspect_or_interacting"].sum())

    def _pct(n: int) -> float:
        return (n / n_total * 100.0) if n_total else 0.0

    print(
        f"  Role-aware pair audit: {n_total:,} pairs -- {n_cc:,} C-C ({_pct(n_cc):.1f}%), "
        f"{n_suspect_or_interacting:,} with >=1 PS/SS/I ({_pct(n_suspect_or_interacting):.1f}%)"
    )

    # ---- per-case aggregation ----
    per_case_src = audit_df[["caseid", "is_cc_pair", "has_suspect_or_interacting"]].copy()
    per_case_src["pair_has_ps"] = pair_has_ps.values
    per_case_src["pair_has_ss"] = pair_has_ss.values
    per_case_src["pair_has_i"] = pair_has_i.values
    per_case = per_case_src.groupby("caseid").agg(
        n_pairs=("caseid", "size"),
        n_cc_pairs=("is_cc_pair", "sum"),
        n_pairs_with_ps=("pair_has_ps", "sum"),
        n_pairs_with_ss=("pair_has_ss", "sum"),
        n_pairs_with_i=("pair_has_i", "sum"),
        n_pairs_suspect_or_interacting=("has_suspect_or_interacting", "sum"),
    ).reset_index()
    per_case["pct_cc_pairs"] = per_case.apply(
        lambda r: (r["n_cc_pairs"] / r["n_pairs"] * 100.0) if r["n_pairs"] else 0.0, axis=1
    )
    per_case["n_pairs_remaining_exclude_cc"] = per_case["n_pairs"] - per_case["n_cc_pairs"]

    def _case_rows(df: pd.DataFrame) -> list:
        return [
            {
                "caseid": str(row["caseid"]),
                "n_pairs": int(row["n_pairs"]),
                "n_cc_pairs": int(row["n_cc_pairs"]),
                "pct_cc_pairs": float(row["pct_cc_pairs"]),
            }
            for _, row in df.iterrows()
        ]

    top_20_by_n_cc = per_case.sort_values("n_cc_pairs", ascending=False).head(20)
    top_20_by_pct_cc = per_case.sort_values(
        ["pct_cc_pairs", "n_pairs"], ascending=[False, False]
    ).head(20)

    known_case_row = per_case[per_case["caseid"] == _KNOWN_AUDIT_CASEID]
    known_case_summary = None
    if len(known_case_row):
        r = known_case_row.iloc[0]
        known_case_summary = {
            "caseid": _KNOWN_AUDIT_CASEID,
            "total_pairs": int(r["n_pairs"]),
            "cc_pairs": int(r["n_cc_pairs"]),
            "pairs_with_ge1_ps_ss_or_i": int(r["n_pairs_suspect_or_interacting"]),
        }

    n_cases_total = int(per_case["caseid"].nunique())
    n_cases_remaining_exclude_cc = int((per_case["n_pairs_remaining_exclude_cc"] > 0).sum())
    n_cases_remaining_require_suspect = int((per_case["n_pairs_suspect_or_interacting"] > 0).sum())

    diagnostics = {
        "drug_id_mode": drug_id_mode,
        "total_pair_rows": n_total,
        "n_cc_pairs": n_cc,
        "pct_cc_pairs": _pct(n_cc),
        "n_pairs_with_ge1_ps": n_with_ps,
        "pct_pairs_with_ge1_ps": _pct(n_with_ps),
        "n_pairs_with_ge1_ss": n_with_ss,
        "pct_pairs_with_ge1_ss": _pct(n_with_ss),
        "n_pairs_with_ge1_i": n_with_i,
        "pct_pairs_with_ge1_i": _pct(n_with_i),
        "n_pairs_with_ge1_ps_ss_or_i": n_suspect_or_interacting,
        "pct_pairs_with_ge1_ps_ss_or_i": _pct(n_suspect_or_interacting),
        "n_cases_in_audit": n_cases_total,
        "per_case_pct_cc_pairs_summary": {
            "mean": float(per_case["pct_cc_pairs"].mean()) if n_cases_total else None,
            "median": float(per_case["pct_cc_pairs"].median()) if n_cases_total else None,
            "min": float(per_case["pct_cc_pairs"].min()) if n_cases_total else None,
            "max": float(per_case["pct_cc_pairs"].max()) if n_cases_total else None,
        },
        "top_20_cases_by_n_cc_pairs": _case_rows(top_20_by_n_cc),
        "top_20_cases_by_pct_cc_pairs": _case_rows(top_20_by_pct_cc),
        f"case_{_KNOWN_AUDIT_CASEID}_summary": known_case_summary,
        # Sensitivity analysis only -- NOT applied to case_drug_pairs.parquet.
        "sensitivity_exclude_cc": {
            "n_pairs_remaining": n_total - n_cc,
            "n_cases_with_ge1_pair_remaining": n_cases_remaining_exclude_cc,
        },
        "sensitivity_require_ps_ss_or_i": {
            "n_pairs_remaining": n_suspect_or_interacting,
            "n_cases_with_ge1_pair_remaining": n_cases_remaining_require_suspect,
        },
    }

    return audit_df, diagnostics


def print_drug_pair_role_audit_diagnostics(diagnostics: dict) -> None:
    d = diagnostics
    print("\n  Role-aware pair audit diagnostics:")
    print(f"    drug_id_mode:                          {d['drug_id_mode']}")
    print(f"    total pair rows:                       {d['total_pair_rows']:,}")
    print(f"    C-C pairs:                              {d['n_cc_pairs']:,} ({d['pct_cc_pairs']:.2f}%)")
    print(f"    pairs with >=1 PS:                      {d['n_pairs_with_ge1_ps']:,} ({d['pct_pairs_with_ge1_ps']:.2f}%)")
    print(f"    pairs with >=1 SS:                      {d['n_pairs_with_ge1_ss']:,} ({d['pct_pairs_with_ge1_ss']:.2f}%)")
    print(f"    pairs with >=1 I:                       {d['n_pairs_with_ge1_i']:,} ({d['pct_pairs_with_ge1_i']:.2f}%)")
    print(f"    pairs with >=1 PS/SS/I:                 {d['n_pairs_with_ge1_ps_ss_or_i']:,} ({d['pct_pairs_with_ge1_ps_ss_or_i']:.2f}%)")
    print(f"    cases in audit:                         {d['n_cases_in_audit']:,}")
    s = d["per_case_pct_cc_pairs_summary"]
    if s["mean"] is not None:
        print(f"    per-case %% C-C pairs: mean={s['mean']:.1f}% median={s['median']:.1f}% min={s['min']:.1f}% max={s['max']:.1f}%")
    print("    top 20 cases by # C-C pairs:")
    for row in d["top_20_cases_by_n_cc_pairs"]:
        print(f"      caseid={row['caseid']:<12} n_cc_pairs={row['n_cc_pairs']:>6,} / n_pairs={row['n_pairs']:>6,} ({row['pct_cc_pairs']:.1f}%)")
    print("    top 20 cases by %% C-C pairs:")
    for row in d["top_20_cases_by_pct_cc_pairs"]:
        print(f"      caseid={row['caseid']:<12} pct_cc={row['pct_cc_pairs']:>6.1f}% (n_pairs={row['n_pairs']:,}, n_cc_pairs={row['n_cc_pairs']:,})")
    known_key = f"case_{_KNOWN_AUDIT_CASEID}_summary"
    if d.get(known_key):
        k = d[known_key]
        print(f"    caseid {_KNOWN_AUDIT_CASEID}: total_pairs={k['total_pairs']:,}  cc_pairs={k['cc_pairs']:,}  pairs_with_ge1_ps_ss_or_i={k['pairs_with_ge1_ps_ss_or_i']:,}")
    print("    sensitivity analysis (NOT applied to case_drug_pairs.parquet):")
    exc = d["sensitivity_exclude_cc"]
    req = d["sensitivity_require_ps_ss_or_i"]
    print(f"      exclude C-C:              {exc['n_pairs_remaining']:,} pairs remain across {exc['n_cases_with_ge1_pair_remaining']:,} cases")
    print(f"      require >=1 PS/SS/I:      {req['n_pairs_remaining']:,} pairs remain across {req['n_cases_with_ge1_pair_remaining']:,} cases")


def run_drug_pair_role_audit_validation(
    case_drug_pairs_df: pd.DataFrame,
    role_audit_df: pd.DataFrame,
    diagnostics: dict,
) -> None:
    """
    Validate case_drug_pair_role_audit. Raises AssertionError listing every
    failed check -- since this is a sensitivity analysis that will inform a
    real filtering decision later, a silently-wrong percentage or pair count
    here would be worse than no audit at all.
    """
    failures = []

    required_cols = [
        "caseid", "drug_a", "drug_b", "roles_a", "roles_b",
        "drug_a_has_ps", "drug_b_has_ps", "drug_a_has_ss", "drug_b_has_ss",
        "drug_a_has_i", "drug_b_has_i", "drug_a_only_c", "drug_b_only_c",
        "is_cc_pair", "has_suspect_or_interacting", "drug_id_mode",
    ]
    missing_cols = [c for c in required_cols if c not in role_audit_df.columns]
    if missing_cols:
        failures.append(f"role_audit_df is missing required column(s): {missing_cols}")
        raise AssertionError("Role-aware pair audit validation failed:\n  - " + "\n  - ".join(failures))

    if role_audit_df["caseid"].isna().any():
        failures.append("role_audit_df has rows with null caseid")
    if role_audit_df["drug_a"].isna().any():
        failures.append("role_audit_df has rows with null drug_a")
    if role_audit_df["drug_b"].isna().any():
        failures.append("role_audit_df has rows with null drug_b")
    if role_audit_df.duplicated(subset=["caseid", "drug_a", "drug_b"]).any():
        failures.append("role_audit_df has duplicate (caseid, drug_a, drug_b) rows")

    if role_audit_df["roles_a"].isna().any() or role_audit_df["roles_b"].isna().any():
        failures.append("role_audit_df has rows where role metadata (roles_a/roles_b) could not be looked up")

    pct_fields = [
        "pct_cc_pairs", "pct_pairs_with_ge1_ps", "pct_pairs_with_ge1_ss",
        "pct_pairs_with_ge1_i", "pct_pairs_with_ge1_ps_ss_or_i",
    ]
    for field in pct_fields:
        value = diagnostics.get(field)
        if value is not None and not (0.0 <= value <= 100.0):
            failures.append(f"diagnostics[{field!r}] = {value} is outside the valid [0, 100] range")

    if len(role_audit_df) != len(case_drug_pairs_df):
        failures.append(
            f"role_audit_df has {len(role_audit_df):,} rows, expected to match "
            f"case_drug_pairs row count {len(case_drug_pairs_df):,}"
        )

    if failures:
        raise AssertionError("Role-aware pair audit validation failed:\n  - " + "\n  - ".join(failures))

    print("  All role-aware pair audit validation checks passed:")
    print("    - role_audit_df has required columns")
    print("    - no missing caseid / drug_a / drug_b")
    print("    - no duplicate (caseid, drug_a, drug_b)")
    print("    - role metadata (roles_a/roles_b) present for every pair")
    print("    - all reported percentages are within [0, 100]")
    print("    - role_audit_df row count matches case_drug_pairs row count")


def save_drug_pair_role_audit_artifacts(
    out_dir: str,
    role_audit_df: pd.DataFrame,
    diagnostics: dict,
) -> None:
    audit_path = os.path.join(out_dir, "case_drug_pair_role_audit.parquet")
    role_audit_df.to_parquet(audit_path, index=False)
    print(f"  Saved: {audit_path}  ({len(role_audit_df):,} rows)")

    diagnostics_path = os.path.join(out_dir, "drug_pair_role_audit_diagnostics.json")
    with open(diagnostics_path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    print(f"  Saved: {diagnostics_path}")


# ---------------------------------------------------------------------------
# Milestone 5d -- configurable pair-role filter (case_drug_pairs_selected.parquet)
# ---------------------------------------------------------------------------

PAIR_ROLE_FILTERS = ("all_pairs", "non_cc")
DEFAULT_PAIR_ROLE_FILTER = "all_pairs"

_SELECTED_PAIR_COLUMNS = [
    "caseid", "drug_a", "drug_b", "n_mapped_drugs_in_case", "drug_id_mode",
    "pair_role_filter", "drug_a_roles_seen", "drug_b_roles_seen",
    "is_cc_pair", "has_ps", "has_ss", "has_i", "has_suspect_or_interacting",
]


def apply_pair_role_filter(
    case_drug_pairs_df: pd.DataFrame,
    case_drug_pair_role_audit_df: pd.DataFrame,
    pair_role_filter: str,
    drug_id_mode: str,
) -> tuple:
    """
    Select which case_drug_pairs rows proceed to the (not-yet-implemented)
    pair-reaction join, based on --pair-role-filter. Does NOT mutate or
    overwrite case_drug_pairs.parquet -- this produces a new, separate
    table (case_drug_pairs_selected.parquet), so the full unfiltered pair
    set and its audit tables remain available.

      'all_pairs' (default, backward-compatible/fully auditable): keep
        every row from case_drug_pairs_df unchanged.
      'non_cc' (recommended for the thesis comparison): drop pairs where
        both drugs are only concomitant (C-C), i.e. keep pairs with >=1
        PS/SS/I drug. Equivalent to case_drug_pair_role_audit_df's
        has_suspect_or_interacting flag (see Milestone 5c, where "exclude
        C-C" and "require >=1 PS/SS/I" were shown to produce identical
        counts, since FAERS only has 4 role codes).

    Role metadata (roles_seen per drug, is_cc_pair, has_ps/ss/i,
    has_suspect_or_interacting) is carried into the selected output from
    case_drug_pair_role_audit_df for downstream traceability -- so it's
    still visible *why* a pair was selected/dropped without re-joining.

    Returns (selected_df, diagnostics_dict).
    """
    if pair_role_filter not in PAIR_ROLE_FILTERS:
        raise ValueError(f"pair_role_filter must be one of {PAIR_ROLE_FILTERS}, got {pair_role_filter!r}")

    _require_columns(case_drug_pairs_df, ["caseid", "drug_a", "drug_b", "n_mapped_drugs_in_case"], "case_drug_pairs")
    _require_columns(
        case_drug_pair_role_audit_df,
        ["caseid", "drug_a", "drug_b", "roles_a", "roles_b", "is_cc_pair",
         "drug_a_has_ps", "drug_b_has_ps", "drug_a_has_ss", "drug_b_has_ss",
         "drug_a_has_i", "drug_b_has_i", "has_suspect_or_interacting"],
        "case_drug_pair_role_audit",
    )

    role_cols = case_drug_pair_role_audit_df[[
        "caseid", "drug_a", "drug_b", "roles_a", "roles_b", "is_cc_pair",
        "drug_a_has_ps", "drug_b_has_ps", "drug_a_has_ss", "drug_b_has_ss",
        "drug_a_has_i", "drug_b_has_i", "has_suspect_or_interacting",
    ]].rename(columns={"roles_a": "drug_a_roles_seen", "roles_b": "drug_b_roles_seen"})

    merged = case_drug_pairs_df.merge(role_cols, on=["caseid", "drug_a", "drug_b"], how="left")
    merged["has_ps"] = merged["drug_a_has_ps"] | merged["drug_b_has_ps"]
    merged["has_ss"] = merged["drug_a_has_ss"] | merged["drug_b_has_ss"]
    merged["has_i"] = merged["drug_a_has_i"] | merged["drug_b_has_i"]
    merged["pair_role_filter"] = pair_role_filter
    merged["drug_id_mode"] = drug_id_mode

    if pair_role_filter == "all_pairs":
        selected_df = merged
    else:  # non_cc
        selected_df = merged[~merged["is_cc_pair"]]

    selected_df = selected_df[_SELECTED_PAIR_COLUMNS].reset_index(drop=True)

    n_input = len(case_drug_pairs_df)
    n_selected = len(selected_df)
    n_removed = n_input - n_selected

    def _pct(n: int) -> float:
        return (n / n_input * 100.0) if n_input else 0.0

    print(
        f"  Pair-role filter '{pair_role_filter}': {n_selected:,} / {n_input:,} pairs selected "
        f"({_pct(n_selected):.1f}%), {n_removed:,} removed ({_pct(n_removed):.1f}%)"
    )

    n_cases_before = int(case_drug_pairs_df["caseid"].nunique())
    n_cases_after = int(selected_df["caseid"].nunique())
    n_unique_pairs_before = int(case_drug_pairs_df[["drug_a", "drug_b"]].drop_duplicates().shape[0])
    n_unique_pairs_after = int(selected_df[["drug_a", "drug_b"]].drop_duplicates().shape[0]) if n_selected else 0

    def _high_polypharmacy_buckets(df: pd.DataFrame) -> dict:
        return {
            "rows_from_cases_ge10": int((df["n_mapped_drugs_in_case"] >= 10).sum()),
            "rows_from_cases_ge20": int((df["n_mapped_drugs_in_case"] >= 20).sum()),
            "rows_from_cases_ge50": int((df["n_mapped_drugs_in_case"] >= 50).sum()),
            "rows_from_cases_ge100": int((df["n_mapped_drugs_in_case"] >= 100).sum()),
        }

    known_case_input = case_drug_pairs_df[case_drug_pairs_df["caseid"] == _KNOWN_AUDIT_CASEID]
    known_case_selected = selected_df[selected_df["caseid"] == _KNOWN_AUDIT_CASEID]
    known_case_summary = None
    if len(known_case_input):
        n_in = len(known_case_input)
        n_sel = len(known_case_selected)
        known_case_summary = {
            "caseid": _KNOWN_AUDIT_CASEID,
            "input_pairs": n_in,
            "selected_pairs": n_sel,
            "removed_pairs": n_in - n_sel,
            "pct_selected": (n_sel / n_in * 100.0) if n_in else 0.0,
            "pct_removed": ((n_in - n_sel) / n_in * 100.0) if n_in else 0.0,
        }

    diagnostics = {
        "pair_role_filter": pair_role_filter,
        "drug_id_mode": drug_id_mode,
        "n_input_pairs": n_input,
        "n_selected_pairs": n_selected,
        "n_removed_pairs": n_removed,
        "pct_selected_pairs": _pct(n_selected),
        "pct_removed_pairs": _pct(n_removed),
        "n_cases_before": n_cases_before,
        "n_cases_after": n_cases_after,
        "n_cases_lost_all_pairs_filtered": n_cases_before - n_cases_after,
        "n_unique_drug_pairs_before": n_unique_pairs_before,
        "n_unique_drug_pairs_after": n_unique_pairs_after,
        "high_polypharmacy_contribution_before": _high_polypharmacy_buckets(case_drug_pairs_df),
        "high_polypharmacy_contribution_after": (
            _high_polypharmacy_buckets(selected_df) if n_selected
            else {"rows_from_cases_ge10": 0, "rows_from_cases_ge20": 0, "rows_from_cases_ge50": 0, "rows_from_cases_ge100": 0}
        ),
        f"case_{_KNOWN_AUDIT_CASEID}_summary": known_case_summary,
    }

    return selected_df, diagnostics


def print_pair_role_filter_diagnostics(diagnostics: dict) -> None:
    d = diagnostics
    print("\n  Pair-role filter diagnostics:")
    print(f"    pair_role_filter:                      {d['pair_role_filter']}")
    print(f"    drug_id_mode:                          {d['drug_id_mode']}")
    print(f"    input pairs:                            {d['n_input_pairs']:,}")
    print(f"    selected pairs:                         {d['n_selected_pairs']:,} ({d['pct_selected_pairs']:.2f}%)")
    print(f"    removed pairs:                          {d['n_removed_pairs']:,} ({d['pct_removed_pairs']:.2f}%)")
    print(f"    cases before:                            {d['n_cases_before']:,}")
    print(f"    cases after:                             {d['n_cases_after']:,}")
    print(f"    cases lost (all pairs filtered out):     {d['n_cases_lost_all_pairs_filtered']:,}")
    print(f"    unique drug pairs before:                {d['n_unique_drug_pairs_before']:,}")
    print(f"    unique drug pairs after:                 {d['n_unique_drug_pairs_after']:,}")
    b, a = d["high_polypharmacy_contribution_before"], d["high_polypharmacy_contribution_after"]
    print("    high-polypharmacy pair-row contribution (before -> after):")
    for key, label in [
        ("rows_from_cases_ge10", ">=10 mapped drugs"),
        ("rows_from_cases_ge20", ">=20 mapped drugs"),
        ("rows_from_cases_ge50", ">=50 mapped drugs"),
        ("rows_from_cases_ge100", ">=100 mapped drugs"),
    ]:
        print(f"      {label:<20} {b[key]:>10,} -> {a[key]:>10,}")
    known_key = f"case_{_KNOWN_AUDIT_CASEID}_summary"
    if d.get(known_key):
        k = d[known_key]
        print(
            f"    caseid {_KNOWN_AUDIT_CASEID}: input={k['input_pairs']:,}  selected={k['selected_pairs']:,} "
            f"({k['pct_selected']:.1f}%)  removed={k['removed_pairs']:,} ({k['pct_removed']:.1f}%)"
        )


def run_pair_role_filter_validation(
    case_drug_pairs_df: pd.DataFrame,
    case_drug_pair_role_audit_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    pair_role_filter: str,
) -> None:
    """
    Validate case_drug_pairs_selected. Raises AssertionError listing every
    failed check -- this table decides what actually reaches the
    pair-reaction join, so a silent bug here (e.g. leaking a C-C pair
    through in non_cc mode) would defeat the whole point of the filter.
    """
    failures = []

    required_cols = [
        "caseid", "drug_a", "drug_b", "n_mapped_drugs_in_case",
        "drug_id_mode", "pair_role_filter",
    ]
    missing_cols = [c for c in required_cols if c not in selected_df.columns]
    if missing_cols:
        failures.append(f"selected_df is missing required column(s): {missing_cols}")
        raise AssertionError("Pair-role filter validation failed:\n  - " + "\n  - ".join(failures))

    if selected_df["caseid"].isna().any():
        failures.append("selected_df has rows with null caseid")
    if selected_df["drug_a"].isna().any():
        failures.append("selected_df has rows with null drug_a")
    if selected_df["drug_b"].isna().any():
        failures.append("selected_df has rows with null drug_b")
    if (selected_df["drug_a"] == selected_df["drug_b"]).any():
        failures.append("selected_df has self-pairs (drug_a == drug_b)")
    if selected_df.duplicated(subset=["caseid", "drug_a", "drug_b"]).any():
        failures.append("selected_df has duplicate (caseid, drug_a, drug_b) rows")

    original_pairs = set(map(tuple, case_drug_pairs_df[["caseid", "drug_a", "drug_b"]].itertuples(index=False, name=None)))
    selected_pairs = set(map(tuple, selected_df[["caseid", "drug_a", "drug_b"]].itertuples(index=False, name=None)))
    if not selected_pairs.issubset(original_pairs):
        failures.append("selected_df contains a (caseid, drug_a, drug_b) not present in case_drug_pairs")

    modes = set(selected_df["pair_role_filter"].dropna().unique())
    if modes != {pair_role_filter}:
        failures.append(
            f"selected_df.pair_role_filter has value(s) {sorted(modes)}, expected only {pair_role_filter!r}"
        )

    if pair_role_filter == "all_pairs":
        if len(selected_df) != len(case_drug_pairs_df):
            failures.append(
                f"all_pairs mode: selected_df has {len(selected_df):,} rows, "
                f"expected to match case_drug_pairs row count {len(case_drug_pairs_df):,}"
            )
    else:  # non_cc
        if "is_cc_pair" in selected_df.columns and selected_df["is_cc_pair"].any():
            failures.append("non_cc mode: selected_df contains a pair classified as C-C (is_cc_pair)")
        n_expected = int(case_drug_pair_role_audit_df["has_suspect_or_interacting"].sum())
        if len(selected_df) != n_expected:
            failures.append(
                f"non_cc mode: selected_df has {len(selected_df):,} rows, expected "
                f"{n_expected:,} (case_drug_pair_role_audit.has_suspect_or_interacting count)"
            )

    if failures:
        raise AssertionError("Pair-role filter validation failed:\n  - " + "\n  - ".join(failures))

    print("  All pair-role filter validation checks passed:")
    print("    - selected_df has required columns")
    print("    - no missing caseid / drug_a / drug_b")
    print("    - no self-pairs")
    print("    - no duplicate (caseid, drug_a, drug_b)")
    print("    - every selected pair exists in case_drug_pairs")
    print("    - pair_role_filter has exactly one value, matching the CLI flag")
    if pair_role_filter == "all_pairs":
        print("    - all_pairs mode: selected row count == case_drug_pairs row count")
    else:
        print("    - non_cc mode: no selected pair is classified as C-C")
        print("    - non_cc mode: selected row count == role-audit has_suspect_or_interacting count")


def save_pair_role_filter_artifacts(
    out_dir: str,
    selected_df: pd.DataFrame,
    diagnostics: dict,
) -> None:
    selected_path = os.path.join(out_dir, "case_drug_pairs_selected.parquet")
    selected_df.to_parquet(selected_path, index=False)
    print(f"  Saved: {selected_path}  ({len(selected_df):,} rows)")

    diagnostics_path = os.path.join(out_dir, "pair_role_filter_diagnostics.json")
    with open(diagnostics_path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    print(f"  Saved: {diagnostics_path}")


# ---------------------------------------------------------------------------
# Milestone 6 -- case_pair_reactions + drug_pair_reaction_counts (co-reporting
# only; NOT causal, NOT PRR)
# ---------------------------------------------------------------------------
#
# TODO (deferred, not addressed by this milestone):
#   - PRR / disproportionality statistics are Milestone 7, not here -- n_cases
#     is a raw co-reporting support count, nothing more
#   - --pair-role-filter all_pairs vs non_cc changes which pairs feed this
#     join (see Milestone 5d); this milestone does not pick one for you
#   - R-GCN export and train/validation/test splitting are both still out of
#     scope

# A case_pair_reactions row means: this (drug_a, drug_b) pair and this
# reaction_pt were both reported within the SAME FAERS case. That is
# co-reporting/co-occurrence evidence only -- it says nothing about which
# drug (if either) caused the reaction, nor about mechanism. Do not read
# "pair X is associated with reaction Y" as a causal statement anywhere
# downstream of this table.

_CASE_PAIR_REACTION_ROLE_COLUMNS = ["is_cc_pair", "has_ps", "has_ss", "has_i", "has_suspect_or_interacting"]
_CASE_PAIR_REACTION_COLUMNS = (
    ["caseid", "drug_a", "drug_b", "reaction_pt", "n_mapped_drugs_in_case", "drug_id_mode", "pair_role_filter"]
    + _CASE_PAIR_REACTION_ROLE_COLUMNS
)
DEFAULT_MAX_ROWS_PER_JOIN_CHUNK = 2_000_000


def _join_pairs_to_reactions(pairs_df: pd.DataFrame, reactions_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pure join step, reused by both the simple (in-memory) and chunked
    builders below: one row per (case's selected pair) x (case's reaction),
    i.e. co-reporting only -- see module comment above.
    """
    reac_cols = reactions_df[["caseid", "reaction_pt"]]
    merged = pairs_df.merge(reac_cols, on="caseid", how="inner")
    for col in _CASE_PAIR_REACTION_ROLE_COLUMNS:
        if col not in merged.columns:
            merged[col] = pd.NA
    return merged[_CASE_PAIR_REACTION_COLUMNS]


def build_case_pair_reactions(
    selected_pairs_df: pd.DataFrame,
    case_reactions_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Join case_drug_pairs_selected with case_reactions on caseid to produce
    case_pair_reactions: co-reporting rows only (this pair and this reaction
    appeared in the same case -- see module comment; NOT a causal claim).

    Simple, single-pass, in-memory reference implementation -- safe for
    small/medium inputs (e.g. the --max-rows 5000 smoke test). The real
    pipeline uses build_case_pair_reactions_chunked() below instead, because
    this join can produce tens of millions of rows (measured: ~41M rows for
    the --max-rows 300000 --pair-role-filter non_cc run) which is not safe
    to hold in memory as a single pandas DataFrame alongside everything else
    already loaded.
    """
    _require_columns(
        selected_pairs_df,
        ["caseid", "drug_a", "drug_b", "n_mapped_drugs_in_case", "drug_id_mode", "pair_role_filter"],
        "case_drug_pairs_selected",
    )
    _require_columns(case_reactions_df, ["caseid", "reaction_pt"], "case_reactions")
    return _join_pairs_to_reactions(selected_pairs_df, case_reactions_df).reset_index(drop=True)


def build_drug_pair_reaction_counts(case_pair_reactions_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate case_pair_reactions to (drug_a, drug_b, reaction_pt) -> n_cases.

    n_cases is the number of DISTINCT caseid values reporting that triple
    (co-reporting support), computed via groupby(...)['caseid'].nunique() --
    deliberately NOT groupby(...).size() / raw row count. In this pipeline
    the two happen to be numerically identical for every triple, because a
    single case can contribute at most one row per (drug_a, drug_b,
    reaction_pt): case_drug_pairs_selected has at most one row per
    (caseid, drug_a, drug_b) (Milestone 5) and case_reactions has at most
    one row per (caseid, reaction_pt) (Milestone 4), so their per-case cross
    product cannot duplicate a triple within one case. nunique() is used
    anyway so the support count is correct *by definition*, not merely by
    this dataset's current upstream guarantees.
    """
    _require_columns(
        case_pair_reactions_df,
        ["caseid", "drug_a", "drug_b", "reaction_pt", "drug_id_mode", "pair_role_filter"],
        "case_pair_reactions",
    )
    counts = (
        case_pair_reactions_df.groupby(["drug_a", "drug_b", "reaction_pt"])["caseid"]
        .nunique()
        .rename("n_cases")
        .reset_index()
    )
    counts["drug_id_mode"] = case_pair_reactions_df["drug_id_mode"].iloc[0] if len(case_pair_reactions_df) else pd.NA
    counts["pair_role_filter"] = (
        case_pair_reactions_df["pair_role_filter"].iloc[0] if len(case_pair_reactions_df) else pd.NA
    )
    return counts[["drug_a", "drug_b", "reaction_pt", "n_cases", "drug_id_mode", "pair_role_filter"]]


def _plan_caseid_chunks(case_stats_df: pd.DataFrame, max_rows_per_chunk: int) -> list:
    """
    Greedily bucket caseids (sorted for determinism) into chunks whose total
    n_pair_reaction_rows_for_case stays <= max_rows_per_chunk, so each
    chunk's join is bounded in size. A single case whose own row count
    already exceeds max_rows_per_chunk gets its own chunk (its rows can't be
    split further without breaking the "one case, one chunk" invariant the
    distinct-caseid-count accumulation below relies on).
    """
    ordered = case_stats_df.sort_values("caseid")
    chunks, current_chunk, current_rows = [], [], 0
    for row in ordered.itertuples():
        case_rows = int(row.n_pair_reaction_rows_for_case)
        if current_chunk and current_rows + case_rows > max_rows_per_chunk:
            chunks.append(current_chunk)
            current_chunk, current_rows = [], 0
        current_chunk.append(row.caseid)
        current_rows += case_rows
    if current_chunk:
        chunks.append(current_chunk)
    return chunks


def build_case_pair_reactions_chunked(
    selected_pairs_df: pd.DataFrame,
    case_reactions_df: pd.DataFrame,
    out_path: str,
    drug_id_mode: str,
    pair_role_filter: str,
    max_rows_per_chunk: int = DEFAULT_MAX_ROWS_PER_JOIN_CHUNK,
) -> tuple:
    """
    Memory-safe version of build_case_pair_reactions() +
    build_drug_pair_reaction_counts(), combined into a single chunked pass:

      1. Compute each case's exact contribution to the join up front
         (n_pairs_for_case * n_reactions_for_case for cases present in both
         inputs) -- cheap, one row per case, no join yet.
      2. Partition caseids into batches bounded by max_rows_per_chunk
         (_plan_caseid_chunks).
      3. For each batch: join just that batch (same _join_pairs_to_reactions
         logic as the simple builder), append its rows directly to out_path
         via a streaming pyarrow ParquetWriter (case_pair_reactions.parquet
         is never held in memory in full), and fold its per-triple row
         counts into a running partial-counts accumulator.

      Because batches partition caseids (each case's rows are produced by
      exactly one batch) and a case contributes at most one row per triple
      (see build_drug_pair_reaction_counts docstring), summing each batch's
      per-triple counts across all batches gives the exact number of
      DISTINCT caseids per triple -- equivalent to a single
      groupby(...).nunique() over the full un-materialized table.

      4. Also validates every batch inline (required columns, no missing
         keys, no self-pairs, no duplicate (caseid, drug_a, drug_b,
         reaction_pt), drug_id_mode/pair_role_filter constant, and in
         non_cc mode no is_cc_pair rows) and accumulates a summary for
         run_case_pair_reactions_validation() -- raises immediately on the
         first batch that fails, same "fail loudly" convention as the rest
         of this file.

    Returns (n_total_rows, counts_df, case_stats_df, row_validation_summary).
    """
    _require_columns(
        selected_pairs_df,
        ["caseid", "drug_a", "drug_b", "n_mapped_drugs_in_case", "drug_id_mode", "pair_role_filter"],
        "case_drug_pairs_selected",
    )
    _require_columns(case_reactions_df, ["caseid", "reaction_pt"], "case_reactions")

    pairs_n = selected_pairs_df.groupby("caseid").size().rename("n_pairs_for_case")
    pairs_mapped = selected_pairs_df.groupby("caseid")["n_mapped_drugs_in_case"].first()
    reac_n = case_reactions_df.groupby("caseid").size().rename("n_reactions_for_case")

    case_stats_df = (
        pd.DataFrame({"n_pairs_for_case": pairs_n})
        .join(pd.DataFrame({"n_mapped_drugs_in_case": pairs_mapped}))
        .join(pd.DataFrame({"n_reactions_for_case": reac_n}), how="inner")  # inner: only cases with both
        .reset_index()
    )
    case_stats_df["n_pair_reaction_rows_for_case"] = (
        case_stats_df["n_pairs_for_case"].astype("int64") * case_stats_df["n_reactions_for_case"].astype("int64")
    )

    chunks = _plan_caseid_chunks(case_stats_df, max_rows_per_chunk)
    print(f"  Joining case_drug_pairs_selected x case_reactions in {len(chunks):,} chunk(s) (cap {max_rows_per_chunk:,} rows/chunk)...")

    schema = pa.schema([
        ("caseid", pa.string()), ("drug_a", pa.string()), ("drug_b", pa.string()),
        ("reaction_pt", pa.string()), ("n_mapped_drugs_in_case", pa.int64()),
        ("drug_id_mode", pa.string()), ("pair_role_filter", pa.string()),
        ("is_cc_pair", pa.bool_()), ("has_ps", pa.bool_()), ("has_ss", pa.bool_()),
        ("has_i", pa.bool_()), ("has_suspect_or_interacting", pa.bool_()),
    ])

    total_rows = 0
    partial_counts = []
    row_validation = {
        "n_missing_caseid": 0, "n_missing_drug_a": 0, "n_missing_drug_b": 0,
        "n_missing_reaction_pt": 0, "n_self_pairs": 0, "n_duplicate_keys": 0,
        "n_cc_pairs_seen_in_non_cc_mode": 0, "drug_id_modes_seen": set(), "pair_role_filters_seen": set(),
    }

    with pq.ParquetWriter(out_path, schema) as writer:
        for chunk_caseids in chunks:
            caseid_set = set(chunk_caseids)
            pairs_subset = selected_pairs_df[selected_pairs_df["caseid"].isin(caseid_set)]
            reac_subset = case_reactions_df[case_reactions_df["caseid"].isin(caseid_set)]
            chunk_df = _join_pairs_to_reactions(pairs_subset, reac_subset)

            row_validation["n_missing_caseid"] += int(chunk_df["caseid"].isna().sum())
            row_validation["n_missing_drug_a"] += int(chunk_df["drug_a"].isna().sum())
            row_validation["n_missing_drug_b"] += int(chunk_df["drug_b"].isna().sum())
            row_validation["n_missing_reaction_pt"] += int(chunk_df["reaction_pt"].isna().sum())
            row_validation["n_self_pairs"] += int((chunk_df["drug_a"] == chunk_df["drug_b"]).sum())
            row_validation["n_duplicate_keys"] += int(
                chunk_df.duplicated(subset=["caseid", "drug_a", "drug_b", "reaction_pt"]).sum()
            )
            if "is_cc_pair" in chunk_df.columns:
                row_validation["n_cc_pairs_seen_in_non_cc_mode"] += (
                    int(chunk_df["is_cc_pair"].fillna(False).astype(bool).sum()) if pair_role_filter == "non_cc" else 0
                )
            row_validation["drug_id_modes_seen"] |= set(chunk_df["drug_id_mode"].dropna().unique())
            row_validation["pair_role_filters_seen"] |= set(chunk_df["pair_role_filter"].dropna().unique())

            chunk_df["is_cc_pair"] = chunk_df["is_cc_pair"].astype("boolean")
            chunk_df["has_ps"] = chunk_df["has_ps"].astype("boolean")
            chunk_df["has_ss"] = chunk_df["has_ss"].astype("boolean")
            chunk_df["has_i"] = chunk_df["has_i"].astype("boolean")
            chunk_df["has_suspect_or_interacting"] = chunk_df["has_suspect_or_interacting"].astype("boolean")

            table = pa.Table.from_pandas(chunk_df, schema=schema, preserve_index=False)
            writer.write_table(table)

            total_rows += len(chunk_df)
            partial_counts.append(
                chunk_df.groupby(["drug_a", "drug_b", "reaction_pt"]).size().rename("n_cases_partial").reset_index()
            )

    if partial_counts:
        all_partial = pd.concat(partial_counts, ignore_index=True)
        counts_df = (
            all_partial.groupby(["drug_a", "drug_b", "reaction_pt"])["n_cases_partial"]
            .sum()
            .rename("n_cases")
            .reset_index()
        )
    else:
        counts_df = pd.DataFrame(columns=["drug_a", "drug_b", "reaction_pt", "n_cases"])
    counts_df["drug_id_mode"] = drug_id_mode
    counts_df["pair_role_filter"] = pair_role_filter
    counts_df = counts_df[["drug_a", "drug_b", "reaction_pt", "n_cases", "drug_id_mode", "pair_role_filter"]]

    print(f"  case_pair_reactions: {total_rows:,} rows written to {out_path}")
    print(f"  drug_pair_reaction_counts: {len(counts_df):,} unique (drug_a, drug_b, reaction_pt) triples")

    return total_rows, counts_df, case_stats_df, row_validation


def compute_pair_reaction_diagnostics(
    selected_pairs_df: pd.DataFrame,
    case_reactions_df: pd.DataFrame,
    total_rows: int,
    counts_df: pd.DataFrame,
    case_stats_df: pd.DataFrame,
    drug_id_mode: str,
    pair_role_filter: str,
) -> dict:
    n_cases_with_pairs = int(selected_pairs_df["caseid"].nunique())
    n_cases_with_reactions = int(case_reactions_df["caseid"].nunique())
    n_cases_with_both = int(len(case_stats_df))

    n_cases = counts_df["n_cases"]
    top_20 = counts_df.sort_values("n_cases", ascending=False).head(20)
    top_20_triples = [
        {
            "drug_a": row["drug_a"], "drug_b": row["drug_b"], "reaction_pt": row["reaction_pt"],
            "n_cases": int(row["n_cases"]),
        }
        for _, row in top_20.iterrows()
    ]

    def _bucket_contribution(min_drugs: int) -> dict:
        rows = int(
            case_stats_df.loc[case_stats_df["n_mapped_drugs_in_case"] >= min_drugs, "n_pair_reaction_rows_for_case"].sum()
        )
        pct = (rows / total_rows * 100.0) if total_rows else 0.0
        return {"rows": rows, "pct": pct}

    sorted_cases = case_stats_df.sort_values("n_pair_reaction_rows_for_case", ascending=False)

    def _top_n_contribution(n: int) -> dict:
        rows = int(sorted_cases["n_pair_reaction_rows_for_case"].head(n).sum())
        pct = (rows / total_rows * 100.0) if total_rows else 0.0
        return {"rows": rows, "pct": pct}

    known_row = case_stats_df[case_stats_df["caseid"] == _KNOWN_AUDIT_CASEID]
    known_summary = None
    if len(known_row):
        r = known_row.iloc[0]
        known_summary = {
            "caseid": _KNOWN_AUDIT_CASEID,
            "selected_pairs": int(r["n_pairs_for_case"]),
            "reactions": int(r["n_reactions_for_case"]),
            "generated_pair_reaction_rows": int(r["n_pair_reaction_rows_for_case"]),
        }

    return {
        "drug_id_mode": drug_id_mode,
        "pair_role_filter": pair_role_filter,
        "n_selected_pair_rows": int(len(selected_pairs_df)),
        "n_case_reaction_rows": int(len(case_reactions_df)),
        "n_cases_with_selected_pairs": n_cases_with_pairs,
        "n_cases_with_reactions": n_cases_with_reactions,
        "n_cases_with_both_pairs_and_reactions": n_cases_with_both,
        "n_case_pair_reaction_rows": int(total_rows),
        "n_unique_drug_pairs": int(counts_df[["drug_a", "drug_b"]].drop_duplicates().shape[0]),
        "n_unique_reaction_pt": int(counts_df["reaction_pt"].nunique()),
        "n_unique_triples": int(len(counts_df)),
        "n_cases_min": int(n_cases.min()) if len(n_cases) else None,
        "n_cases_median": float(n_cases.median()) if len(n_cases) else None,
        "n_cases_max": int(n_cases.max()) if len(n_cases) else None,
        "n_triples_with_n_cases_ge3": int((n_cases >= 3).sum()),
        "n_triples_with_n_cases_ge5": int((n_cases >= 5).sum()),
        "n_triples_with_n_cases_ge10": int((n_cases >= 10).sum()),
        "top_20_triples_by_n_cases": top_20_triples,
        "contribution_from_cases_with_mapped_drugs_ge10": _bucket_contribution(10),
        "contribution_from_cases_with_mapped_drugs_ge20": _bucket_contribution(20),
        "contribution_from_cases_with_mapped_drugs_ge50": _bucket_contribution(50),
        "contribution_from_cases_with_mapped_drugs_ge100": _bucket_contribution(100),
        "contribution_from_top_1_case": _top_n_contribution(1),
        "contribution_from_top_10_cases": _top_n_contribution(10),
        "contribution_from_top_100_cases": _top_n_contribution(100),
        f"case_{_KNOWN_AUDIT_CASEID}_summary": known_summary,
    }


def print_pair_reaction_diagnostics(diagnostics: dict) -> None:
    d = diagnostics
    print("\n  Pair-reaction count diagnostics:")
    print(f"    drug_id_mode:                          {d['drug_id_mode']}")
    print(f"    pair_role_filter:                      {d['pair_role_filter']}")
    print(f"    selected pair rows (input):             {d['n_selected_pair_rows']:,}")
    print(f"    case reaction rows (input):             {d['n_case_reaction_rows']:,}")
    print(f"    cases with selected pairs:              {d['n_cases_with_selected_pairs']:,}")
    print(f"    cases with reactions:                   {d['n_cases_with_reactions']:,}")
    print(f"    cases with both:                        {d['n_cases_with_both_pairs_and_reactions']:,}")
    print(f"    case_pair_reactions rows:               {d['n_case_pair_reaction_rows']:,}")
    print(f"    unique drug pairs:                      {d['n_unique_drug_pairs']:,}")
    print(f"    unique reaction_pt:                     {d['n_unique_reaction_pt']:,}")
    print(f"    unique (drug_a, drug_b, reaction_pt):   {d['n_unique_triples']:,}")
    print(f"    n_cases min/median/max:                 {d['n_cases_min']} / {d['n_cases_median']} / {d['n_cases_max']}")
    print(f"    triples with n_cases >= 3:              {d['n_triples_with_n_cases_ge3']:,}")
    print(f"    triples with n_cases >= 5:              {d['n_triples_with_n_cases_ge5']:,}")
    print(f"    triples with n_cases >= 10:             {d['n_triples_with_n_cases_ge10']:,}")
    print("    high-polypharmacy contribution to case_pair_reaction rows:")
    for key, label in [
        ("contribution_from_cases_with_mapped_drugs_ge10", ">=10 mapped drugs"),
        ("contribution_from_cases_with_mapped_drugs_ge20", ">=20 mapped drugs"),
        ("contribution_from_cases_with_mapped_drugs_ge50", ">=50 mapped drugs"),
        ("contribution_from_cases_with_mapped_drugs_ge100", ">=100 mapped drugs"),
    ]:
        c = d[key]
        print(f"      {label:<20} {c['rows']:>12,} rows ({c['pct']:.2f}%)")
    print("    contribution from top-N cases by generated rows:")
    for key, label in [
        ("contribution_from_top_1_case", "top 1 case"),
        ("contribution_from_top_10_cases", "top 10 cases"),
        ("contribution_from_top_100_cases", "top 100 cases"),
    ]:
        c = d[key]
        print(f"      {label:<15} {c['rows']:>12,} rows ({c['pct']:.2f}%)")
    known_key = f"case_{_KNOWN_AUDIT_CASEID}_summary"
    if d.get(known_key):
        k = d[known_key]
        print(
            f"    caseid {_KNOWN_AUDIT_CASEID}: selected_pairs={k['selected_pairs']:,}  "
            f"reactions={k['reactions']:,}  generated_rows={k['generated_pair_reaction_rows']:,}"
        )
    print("    top 20 triples by n_cases:")
    for row in d["top_20_triples_by_n_cases"]:
        print(f"      {row['drug_a']:<16} {row['drug_b']:<16} {row['reaction_pt']:<35} n_cases={row['n_cases']:,}")


def run_case_pair_reactions_validation(
    row_validation_summary: dict,
    pair_role_filter: str,
    drug_id_mode: str,
) -> None:
    """
    Validate case_pair_reactions from the per-chunk summary accumulated by
    build_case_pair_reactions_chunked() (the full table is never
    materialized, so checks run incrementally per chunk and are summarized
    here). Raises AssertionError listing every failed check.
    """
    failures = []
    s = row_validation_summary

    if s["n_missing_caseid"]:
        failures.append(f"{s['n_missing_caseid']:,} case_pair_reactions rows have null caseid")
    if s["n_missing_drug_a"]:
        failures.append(f"{s['n_missing_drug_a']:,} case_pair_reactions rows have null drug_a")
    if s["n_missing_drug_b"]:
        failures.append(f"{s['n_missing_drug_b']:,} case_pair_reactions rows have null drug_b")
    if s["n_missing_reaction_pt"]:
        failures.append(f"{s['n_missing_reaction_pt']:,} case_pair_reactions rows have null reaction_pt")
    if s["n_self_pairs"]:
        failures.append(f"{s['n_self_pairs']:,} case_pair_reactions rows are self-pairs (drug_a == drug_b)")
    if s["n_duplicate_keys"]:
        failures.append(f"{s['n_duplicate_keys']:,} duplicate (caseid, drug_a, drug_b, reaction_pt) rows found")
    if s["drug_id_modes_seen"] and s["drug_id_modes_seen"] != {drug_id_mode}:
        failures.append(f"case_pair_reactions.drug_id_mode has value(s) {sorted(s['drug_id_modes_seen'])}, expected only {drug_id_mode!r}")
    if s["pair_role_filters_seen"] and s["pair_role_filters_seen"] != {pair_role_filter}:
        failures.append(f"case_pair_reactions.pair_role_filter has value(s) {sorted(s['pair_role_filters_seen'])}, expected only {pair_role_filter!r}")
    if pair_role_filter == "non_cc" and s["n_cc_pairs_seen_in_non_cc_mode"]:
        failures.append(f"non_cc mode: {s['n_cc_pairs_seen_in_non_cc_mode']:,} case_pair_reactions rows are C-C pairs")

    if failures:
        raise AssertionError("case_pair_reactions validation failed:\n  - " + "\n  - ".join(failures))

    print("  All case_pair_reactions validation checks passed (accumulated across all join chunks):")
    print("    - required columns present in every chunk (caseid, drug_a, drug_b, reaction_pt, ...)")
    print("    - no missing caseid / drug_a / drug_b / reaction_pt")
    print("    - no self-pairs")
    print("    - no duplicate (caseid, drug_a, drug_b, reaction_pt)")
    print("    - every caseid present in both case_drug_pairs_selected and case_reactions (by construction)")
    print("    - drug_id_mode consistent across all rows")
    print("    - pair_role_filter consistent across all rows")
    if pair_role_filter == "non_cc":
        print("    - non_cc mode: no row is classified as a C-C pair")


def run_drug_pair_reaction_counts_validation(
    counts_df: pd.DataFrame,
    case_pair_reactions_path: str,
    sample_size: int = 20,
) -> None:
    """
    Validate drug_pair_reaction_counts. Raises AssertionError listing every
    failed check. Includes a deterministic cross-check: re-reads the actual
    saved case_pair_reactions.parquet (via pyarrow predicate pushdown, so
    only matching rows are scanned/loaded) for a fixed sample of triples and
    confirms n_cases == that sample's case_pair_reactions.groupby(...)
    ['caseid'].nunique() -- an independent check on the chunked map-reduce
    aggregation used to build counts_df.
    """
    failures = []

    required_cols = ["drug_a", "drug_b", "reaction_pt", "n_cases", "drug_id_mode", "pair_role_filter"]
    missing_cols = [c for c in required_cols if c not in counts_df.columns]
    if missing_cols:
        failures.append(f"counts_df is missing required column(s): {missing_cols}")
        raise AssertionError("drug_pair_reaction_counts validation failed:\n  - " + "\n  - ".join(failures))

    if counts_df["drug_a"].isna().any():
        failures.append("counts_df has rows with null drug_a")
    if counts_df["drug_b"].isna().any():
        failures.append("counts_df has rows with null drug_b")
    if counts_df["reaction_pt"].isna().any():
        failures.append("counts_df has rows with null reaction_pt")
    if counts_df.duplicated(subset=["drug_a", "drug_b", "reaction_pt"]).any():
        failures.append("counts_df has duplicate (drug_a, drug_b, reaction_pt) rows")
    if (counts_df["n_cases"] <= 0).any():
        failures.append("counts_df has rows with n_cases <= 0 (support count must be positive)")

    # Deterministic sample: top 20 by n_cases (highest-support triples --
    # the ones most likely to matter downstream, and the ones a chunking
    # bug would most plausibly get wrong since they're aggregated from the
    # most chunks/cases).
    sample = counts_df.sort_values(
        ["n_cases", "drug_a", "drug_b", "reaction_pt"], ascending=[False, True, True, True]
    ).head(sample_size)

    if len(sample) and os.path.exists(case_pair_reactions_path):
        filters = [
            [("drug_a", "==", row["drug_a"]), ("drug_b", "==", row["drug_b"]), ("reaction_pt", "==", row["reaction_pt"])]
            for _, row in sample.iterrows()
        ]
        actual = pd.read_parquet(
            case_pair_reactions_path, columns=["caseid", "drug_a", "drug_b", "reaction_pt"], filters=filters,
        )
        actual_nunique = (
            actual.groupby(["drug_a", "drug_b", "reaction_pt"])["caseid"].nunique().rename("n_cases_actual")
        )
        actual_rowcount = (
            actual.groupby(["drug_a", "drug_b", "reaction_pt"]).size().rename("n_rows_actual")
        )
        check = sample.set_index(["drug_a", "drug_b", "reaction_pt"]).join(actual_nunique).join(actual_rowcount)
        mismatched = check[check["n_cases"] != check["n_cases_actual"]]
        if len(mismatched):
            failures.append(
                f"{len(mismatched):,} sampled triples have counts_df.n_cases != re-read "
                f"case_pair_reactions.groupby(...).caseid.nunique() (sample: {mismatched.index[:3].tolist()})"
            )
        # By construction (see build_drug_pair_reaction_counts docstring) raw row
        # count and distinct-caseid count coincide in this dataset; confirmed here
        # rather than assumed, but n_cases itself is always the nunique() value.
        row_vs_nunique_mismatch = check[check["n_rows_actual"] != check["n_cases_actual"]]
        if len(row_vs_nunique_mismatch):
            print(
                f"    [NOTE] {len(row_vs_nunique_mismatch):,} sampled triples have raw row count != "
                f"distinct caseid count -- a case is contributing >1 row to the same triple; "
                f"n_cases (nunique-based) remains correct, this is just informational"
            )

    if failures:
        raise AssertionError("drug_pair_reaction_counts validation failed:\n  - " + "\n  - ".join(failures))

    print("  All drug_pair_reaction_counts validation checks passed:")
    print("    - required columns present")
    print("    - no missing drug_a / drug_b / reaction_pt")
    print("    - no duplicate (drug_a, drug_b, reaction_pt)")
    print("    - n_cases > 0 for every triple")
    print(f"    - {len(sample):,} sampled triples' n_cases verified against a fresh, filtered re-read of "
          f"case_pair_reactions.parquet (groupby(...).caseid.nunique(), not raw row count)")


def save_pair_reaction_artifacts(
    out_dir: str,
    counts_df: pd.DataFrame,
    diagnostics: dict,
) -> None:
    """
    case_pair_reactions.parquet is already written incrementally by
    build_case_pair_reactions_chunked() -- this only saves the (much
    smaller) aggregated counts table, its minimum-support filtered
    variants, and the diagnostics.
    """
    counts_path = os.path.join(out_dir, "drug_pair_reaction_counts.parquet")
    counts_df.to_parquet(counts_path, index=False)
    print(f"  Saved: {counts_path}  ({len(counts_df):,} rows)")

    for min_support in (3, 5, 10):
        filtered = counts_df[counts_df["n_cases"] >= min_support]
        filtered_path = os.path.join(out_dir, f"drug_pair_reaction_counts_min{min_support}.parquet")
        filtered.to_parquet(filtered_path, index=False)
        print(f"  Saved: {filtered_path}  ({len(filtered):,} rows, n_cases >= {min_support})")

    diagnostics_path = os.path.join(out_dir, "pair_reaction_count_diagnostics.json")
    with open(diagnostics_path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    print(f"  Saved: {diagnostics_path}")


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
            "once, ahead of drug-only R-GCN construction. Milestone 6: "
            "builds case_pair_reactions.parquet and "
            "drug_pair_reaction_counts.parquet (co-reporting counts only, "
            "NOT PRR/causal) from --pair-role-filter-selected pairs -- see "
            "file docstring for milestone status."
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
        "--pair-role-filter", choices=PAIR_ROLE_FILTERS, default=DEFAULT_PAIR_ROLE_FILTER,
        help=(
            "Which case_drug_pairs rows are written to "
            "case_drug_pairs_selected.parquet (Milestone 5d); does not "
            "modify case_drug_pairs.parquet itself. 'all_pairs' (default, "
            "backward-compatible/fully auditable) keeps every pair. "
            "'non_cc' (RECOMMENDED for the thesis comparison) drops "
            "concomitant-concomitant (C-C) pairs, keeping only pairs with "
            "at least one PS/SS/I drug."
        ),
    )
    parser.add_argument(
        "--max-rows-per-join-chunk", type=int, default=DEFAULT_MAX_ROWS_PER_JOIN_CHUNK,
        help=(
            "Milestone 6: caps how many case_pair_reactions rows any single "
            "case_drug_pairs_selected x case_reactions join batch may "
            f"produce, so the full join (measured: ~41M rows on the "
            f"--max-rows 300000 --pair-role-filter non_cc run) is never held "
            f"in memory at once. Default: {DEFAULT_MAX_ROWS_PER_JOIN_CHUNK:,}."
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
    print(f"  Pair role filter:{args.pair_role_filter}" + ("  (recommended for thesis comparison)" if args.pair_role_filter == "non_cc" else ""))
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

    print("\n--- Milestone 5: build case_drug_pairs ---")
    case_drug_pairs_df, drug_pair_diagnostics, high_polypharmacy_df = build_case_drug_pairs(
        case_drugs_df, args.drug_id_mode,
    )
    print_drug_pair_diagnostics(drug_pair_diagnostics)

    print("\n--- Milestone 5: case_drug_pairs validation checks ---")
    run_case_drug_pairs_validation(case_drugs_df, case_drug_pairs_df, args.drug_id_mode)

    print("\n--- Milestone 5: saving case_drug_pairs artifacts ---")
    save_case_drug_pairs_artifacts(
        args.out_dir, case_drug_pairs_df, drug_pair_diagnostics, high_polypharmacy_df,
    )

    print("\n--- Milestone 5b: expected pair-reaction join-size audit ---")
    join_audit_df, join_audit_diagnostics = build_expected_pair_reaction_join_audit(
        case_drug_pairs_df, case_reactions_df,
    )
    print_expected_pair_reaction_join_diagnostics(join_audit_diagnostics)

    print("\n--- Milestone 5b: join-size audit validation checks ---")
    run_expected_pair_reaction_join_audit_validation(case_drug_pairs_df, case_reactions_df, join_audit_df)

    print("\n--- Milestone 5b: saving join-size audit artifacts ---")
    save_expected_pair_reaction_join_audit_artifacts(args.out_dir, join_audit_df, join_audit_diagnostics)

    print("\n--- Milestone 5c: role-aware pair audit ---")
    role_audit_df, role_audit_diagnostics = build_case_drug_pair_role_audit(
        case_drug_pairs_df, case_drugs_df, args.drug_id_mode,
    )
    print_drug_pair_role_audit_diagnostics(role_audit_diagnostics)

    print("\n--- Milestone 5c: role-aware pair audit validation checks ---")
    run_drug_pair_role_audit_validation(case_drug_pairs_df, role_audit_df, role_audit_diagnostics)

    print("\n--- Milestone 5c: saving role-aware pair audit artifacts ---")
    save_drug_pair_role_audit_artifacts(args.out_dir, role_audit_df, role_audit_diagnostics)

    print(f"\n--- Milestone 5d: apply pair-role filter (pair_role_filter={args.pair_role_filter}) ---")
    selected_pairs_df, pair_filter_diagnostics = apply_pair_role_filter(
        case_drug_pairs_df, role_audit_df, args.pair_role_filter, args.drug_id_mode,
    )
    print_pair_role_filter_diagnostics(pair_filter_diagnostics)

    print("\n--- Milestone 5d: pair-role filter validation checks ---")
    run_pair_role_filter_validation(case_drug_pairs_df, role_audit_df, selected_pairs_df, args.pair_role_filter)

    print("\n--- Milestone 5d: saving pair-role filter artifacts ---")
    save_pair_role_filter_artifacts(args.out_dir, selected_pairs_df, pair_filter_diagnostics)

    print(f"\n--- Milestone 6: build case_pair_reactions (co-reporting only, chunked) ---")
    case_pair_reactions_path = os.path.join(args.out_dir, "case_pair_reactions.parquet")
    n_pair_reaction_rows, drug_pair_reaction_counts_df, pair_reaction_case_stats_df, row_validation_summary = (
        build_case_pair_reactions_chunked(
            selected_pairs_df, case_reactions_df, case_pair_reactions_path,
            args.drug_id_mode, args.pair_role_filter,
            max_rows_per_chunk=args.max_rows_per_join_chunk,
        )
    )

    print("\n--- Milestone 6: case_pair_reactions validation checks ---")
    run_case_pair_reactions_validation(row_validation_summary, args.pair_role_filter, args.drug_id_mode)

    print("\n--- Milestone 6: pair-reaction diagnostics ---")
    pair_reaction_diagnostics = compute_pair_reaction_diagnostics(
        selected_pairs_df, case_reactions_df, n_pair_reaction_rows, drug_pair_reaction_counts_df,
        pair_reaction_case_stats_df, args.drug_id_mode, args.pair_role_filter,
    )
    print_pair_reaction_diagnostics(pair_reaction_diagnostics)

    print("\n--- Milestone 6: drug_pair_reaction_counts validation checks ---")
    run_drug_pair_reaction_counts_validation(drug_pair_reaction_counts_df, case_pair_reactions_path)

    print("\n--- Milestone 6: saving drug_pair_reaction_counts artifacts ---")
    save_pair_reaction_artifacts(args.out_dir, drug_pair_reaction_counts_df, pair_reaction_diagnostics)

    print("\n--- Milestone 7 (pending): global validation checks ---")
    run_validation_checks(demo_df, retained_cases_df, case_reactions_df, case_drugs_df)

    print("\n" + "=" * 60)
    print("Milestone 6 complete: case_pair_reactions.parquet and drug_pair_reaction_counts.parquet")
    print("(+ min3/min5/min10 support-filtered variants) built from case_drug_pairs_selected.parquet")
    print(f"(drug_id_mode={args.drug_id_mode!r}, pair_role_filter={args.pair_role_filter!r}).")
    print("These are CO-REPORTING counts only -- not PRR, not causal. No R-GCN")
    print("export or train/validation/test split was performed yet.")
    print(f"Artifacts in: {args.out_dir}")
    print("Next: implement Milestone 7 (PRR / disproportionality statistics).")
    print("=" * 60)


if __name__ == "__main__":
    main()
