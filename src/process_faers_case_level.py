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
  [x] Milestone 1 -- CLI skeleton, path wiring, placeholders        (this file)
  [ ] Milestone 2 -- load DEMO/DRUG/REAC; resolve CASEID/CASEVERSION
  [ ] Milestone 3 -- filter DRUG/REAC to retained primaryids; dedupe reactions and drugs per case
  [ ] Milestone 4 -- canonical (drug_i, drug_j) pair generation per case
  [ ] Milestone 5 -- drug_pair_reaction_counts aggregation (PRR comes later, not in Milestone 5 either)
  [ ] Milestone 6 -- validation checks (distinct CASEID counts, no duplicate caseids, etc.)

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
# Milestone 2 (TODO) -- load DEMO / DRUG / REAC and resolve CASEID/CASEVERSION
# ---------------------------------------------------------------------------

def load_demo_table(zip_path: str, inner_path: str, max_rows: int = None) -> pd.DataFrame:
    """
    TODO (Milestone 2): Load the FAERS DEMO table.

    Expected to return at least: primaryid, caseid, caseversion.
    Can reuse process_faers.load_faers_table(zip_path, inner_path, max_rows)
    directly since it is a generic dollar-delimited FAERS table loader.
    """
    print("  [TODO] load_demo_table() not yet implemented (Milestone 2)")
    return pd.DataFrame(columns=["primaryid", "caseid", "caseversion"])


def resolve_case_versions(demo_df: pd.DataFrame) -> pd.DataFrame:
    """
    TODO (Milestone 2): For each caseid, keep only the primaryid belonging to
    the highest caseversion (the most recent amendment/resubmission),
    producing the authoritative primaryid -> caseid mapping ("retained
    primaryids") used to filter DRUG and REAC downstream.
    """
    print("  [TODO] resolve_case_versions() not yet implemented (Milestone 2)")
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
            "once, ahead of drug-only R-GCN construction. Milestone 1: CLI "
            "and placeholders only -- see file docstring for milestone status."
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

    print("\n--- Milestone 2 (pending): load DEMO/DRUG/REAC; resolve CASEID/CASEVERSION ---")
    demo_df = load_demo_table(args.faers_zip, _cfg.FAERS_DEMO_INNER, max_rows=args.max_rows)
    case_map_df = resolve_case_versions(demo_df)

    print("\n--- Milestone 3 (pending): filter to retained PRIMARYIDs; dedupe reactions/drugs ---")
    case_reactions_df = build_case_reactions(pd.DataFrame(), case_map_df)
    case_drugs_df = build_case_drugs(pd.DataFrame(), case_map_df)

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
    print("Milestone 1 complete: pipeline skeleton is wired and runnable.")
    print("No case-level data was actually processed yet.")
    print("Next: implement Milestone 2 (DEMO/DRUG/REAC + CASEID/CASEVERSION).")
    print("=" * 60)


if __name__ == "__main__":
    main()
