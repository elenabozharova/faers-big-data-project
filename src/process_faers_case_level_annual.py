"""
process_faers_case_level_annual.py -- Annual FAERS case-level wrapper (POC).

This is a thin orchestration layer around process_faers_case_level.py. It
does NOT reimplement or modify any of that file's scientific logic --
CASEID/CASEVERSION resolution, case-level reaction dedup, drug dedup,
canonical pair generation, role-aware pair filtering, the chunked
pair-reaction join, and support counting are all called as-is via
`import process_faers_case_level as pfcl`.

Why a year-level wrapper instead of running the existing script four times
and merging its output
---------------------------------------------------------------------------
process_faers_case_level.py resolves one retained PRIMARYID per CASEID from
whatever DEMO/DRUG/REAC rows are loaded in a single run (resolve_case_versions).
FAERS cases are frequently amended in a *later* quarter than the one they were
first reported in (higher CASEVERSION). Running the pipeline separately per
quarter and then summing/merging the four quarterly drug_pair_reaction_counts
tables would resolve the same real-world case as "retained" once per quarter
it appears in, inflating n_cases. This wrapper instead concatenates all four
quarters' raw DEMO/DRUG/REAC for a year BEFORE calling resolve_case_versions,
so CASEID/CASEVERSION resolution and every downstream dedup/count happens
exactly once per year, over the full year's data.

Output
------
For each requested year, writes the same per-run artifacts
process_faers_case_level.py already produces (Milestones 2-6) under
  data/processed/faers_case_level_annual/<year>/
plus:
  - a `year` column added to drug_pair_reaction_counts.parquet and its
    _min3/_min5/_min10 variants
  - annual_summary.json (counts, uniques, duplicate/null checks, runtime,
    output size -- see _build_annual_summary for the exact field list)
retained_cases.parquet (one row per real-world case for that year, with its
CASEID) is preserved per year, un-modified, so a later step can check
whether the same CASEID reappears across different years -- that check
itself is out of scope here.

Explicitly out of scope (unchanged from process_faers_case_level.py):
  - PRR / disproportionality statistics (Milestone 7)
  - any modeling / temporal-persistence prediction layer
  - Spark, Dask, R-GCN, Decagon/ChEMBL/STRING enrichment, Kafka

Usage
-----
  python src/process_faers_case_level_annual.py --years 2025 --max-rows-per-quarter 2000
  python src/process_faers_case_level_annual.py --years 2023 2024 2025
"""

import argparse
import glob
import json
import os
import sys
import time

import pandas as pd

_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(_SRC_DIR)
DATA_RAW = os.path.join(PROJECT_ROOT, "data", "raw")
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")

if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

# Reused as-is (all pipeline logic; nothing in this module is copied or
# reimplemented -- see file docstring).
import process_faers_case_level as pfcl  # noqa: E402

# Assumption: quarterly FAERS zips live directly under data/raw/faers/,
# named faers_ascii_<year><Q|q><1-4>.zip (observed casing is inconsistent,
# e.g. faers_ascii_2023Q3.zip vs faers_ascii_2023q1.zip -- matched
# case-insensitively below). Each zip's inner ASCII files are assumed to
# follow FAERS's standard ASCII/{DEMO,DRUG,REAC}<YY>Q<N>.txt naming
# regardless of the zip filename's own casing (confirmed against the zips
# present in this repo).
FAERS_RAW_DIR = os.path.join(DATA_RAW, "faers")
DEFAULT_OUT_DIR = os.path.join(DATA_PROCESSED, "faers_case_level_annual")
DEFAULT_MAPPING_FILE = pfcl.DEFAULT_MAPPING_FILE

QUARTERS = ("Q1", "Q2", "Q3", "Q4")


# ---------------------------------------------------------------------------
# Locate + load one year's worth of quarterly zips
# ---------------------------------------------------------------------------

def _find_quarter_zip(year: int, quarter: str) -> str:
    """
    Locate the single zip for (year, quarter) under FAERS_RAW_DIR, matching
    the 'Q'/'q' case-insensitively. Returns None if not found (caller
    aggregates missing quarters into one error covering the whole year).
    Raises if more than one zip matches (ambiguous, should never happen).
    """
    quarter_digit = quarter[-1]
    pattern = os.path.join(FAERS_RAW_DIR, f"faers_ascii_{year}[Qq]{quarter_digit}.zip")
    matches = sorted(glob.glob(pattern))
    if len(matches) > 1:
        raise FileNotFoundError(
            f"Ambiguous FAERS zip for {year}{quarter}: multiple files matched {pattern!r}: {matches}"
        )
    return matches[0] if matches else None


def _inner_paths(year: int, quarter: str) -> dict:
    yy = str(year)[2:]
    return {
        "demo": f"ASCII/DEMO{yy}{quarter}.txt",
        "drug": f"ASCII/DRUG{yy}{quarter}.txt",
        "reac": f"ASCII/REAC{yy}{quarter}.txt",
    }


def load_year_raw_tables(year: int, max_rows_per_quarter: int = None) -> tuple:
    """
    Load and concatenate DEMO/DRUG/REAC across all four quarters of `year`.

    Fails loudly (FileNotFoundError) if any of the four quarterly zips is
    missing. A missing DEMO/DRUG/REAC file *inside* a present zip, or a zip
    missing a required column, also fails loudly -- via KeyError from
    zipfile.open() (process_faers.load_faers_table) or the ValueError
    process_faers_case_level._require_columns already raises; neither is
    caught here, so the run stops immediately rather than silently skipping
    a quarter.

    Returns (demo_df, drug_df, reac_df) with a `source_quarter` column
    (e.g. '2023Q3') added for provenance.
    """
    missing = []
    zip_paths = {}
    for quarter in QUARTERS:
        zip_path = _find_quarter_zip(year, quarter)
        if zip_path is None:
            missing.append(quarter)
        else:
            zip_paths[quarter] = zip_path

    if missing:
        raise FileNotFoundError(
            f"Year {year}: missing quarter(s) {missing} under {FAERS_RAW_DIR!r} "
            f"(expected faers_ascii_{year}[Qq]<1-4>.zip for all of {QUARTERS}). "
            f"Found: {sorted(zip_paths.values())}"
        )

    demo_parts, drug_parts, reac_parts = [], [], []
    for quarter in QUARTERS:
        zip_path = zip_paths[quarter]
        inner = _inner_paths(year, quarter)
        print(f"\n  -- {year}{quarter}: {os.path.basename(zip_path)} --")
        demo_q = pfcl.load_demo_table(zip_path, None, inner["demo"], max_rows=max_rows_per_quarter)
        drug_q = pfcl.load_drug_table(zip_path, None, inner["drug"], max_rows=max_rows_per_quarter)
        reac_q = pfcl.load_reac_table(zip_path, None, inner["reac"], max_rows=max_rows_per_quarter)
        for df in (demo_q, drug_q, reac_q):
            df["source_quarter"] = f"{year}{quarter}"
        demo_parts.append(demo_q)
        drug_parts.append(drug_q)
        reac_parts.append(reac_q)

    demo_df = pd.concat(demo_parts, ignore_index=True)
    drug_df = pd.concat(drug_parts, ignore_index=True)
    reac_df = pd.concat(reac_parts, ignore_index=True)
    print(
        f"\n  Year {year} concatenated (Q1-Q4): "
        f"DEMO={len(demo_df):,}  DRUG={len(drug_df):,}  REAC={len(reac_df):,} rows"
    )
    return demo_df, drug_df, reac_df


# ---------------------------------------------------------------------------
# Extra fail-loud checks specific to the annual rollup
# ---------------------------------------------------------------------------

def _validate_annual_counts(counts_df: pd.DataFrame, year: int) -> None:
    """
    Beyond process_faers_case_level.run_drug_pair_reaction_counts_validation
    (already run on this table before the `year` column is added): confirm
    the year-qualified key (year, drug_a, drug_b, reaction_pt) is unique and
    n_cases is a valid positive support count. Raises AssertionError listing
    every failed check.
    """
    failures = []

    n_dup = int(counts_df.duplicated(subset=["year", "drug_a", "drug_b", "reaction_pt"]).sum())
    if n_dup:
        failures.append(f"{n_dup:,} duplicate (year, drug_a, drug_b, reaction_pt) rows")

    if counts_df["n_cases"].isna().any():
        failures.append("n_cases has null values")
    if (counts_df["n_cases"] < 1).any():
        failures.append("n_cases has values < 1")

    if failures:
        raise AssertionError(
            f"Year {year} annual counts validation failed:\n  - " + "\n  - ".join(failures)
        )

    print(f"  Year {year}: annual counts validation passed "
          f"(unique (year, drug_a, drug_b, reaction_pt), n_cases >= 1 for all rows)")


# ---------------------------------------------------------------------------
# Annual summary
# ---------------------------------------------------------------------------

def _dir_size_mb(path: str) -> float:
    total_bytes = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            total_bytes += os.path.getsize(os.path.join(root, name))
    return total_bytes / (1024 * 1024)


def _build_annual_summary(
    year: int,
    year_out_dir: str,
    runtime_seconds: float,
    raw_diagnostics: dict,
    case_resolution_diagnostics: dict,
    reaction_dedup_diagnostics: dict,
    case_drugs_df: pd.DataFrame,
    case_drug_pairs_df: pd.DataFrame,
    selected_pairs_df: pd.DataFrame,
    n_pair_reaction_rows: int,
    counts_df: pd.DataFrame,
) -> dict:
    n_cases = counts_df["n_cases"]
    return {
        "year": year,
        "raw_demo_rows": raw_diagnostics["demo"]["n_rows"],
        "raw_drug_rows": raw_diagnostics["drug"]["n_rows"],
        "raw_reac_rows": raw_diagnostics["reac"]["n_rows"],
        "retained_cases": case_resolution_diagnostics["n_retained_cases"],
        "resolved_drug_rows": case_resolution_diagnostics["n_drug_rows_resolved"],
        "resolved_reac_rows": case_resolution_diagnostics["n_reac_rows_resolved"],
        "case_reactions_rows": reaction_dedup_diagnostics["n_rows_after_dedup"],
        "case_drugs_rows": int(len(case_drugs_df)),
        "case_drug_pairs_rows": int(len(case_drug_pairs_df)),
        "case_drug_pairs_selected_rows": int(len(selected_pairs_df)),
        "case_pair_reactions_rows": int(n_pair_reaction_rows),
        "final_triples": int(len(counts_df)),
        "final_triples_min3": int((n_cases >= 3).sum()),
        "final_triples_min5": int((n_cases >= 5).sum()),
        "final_triples_min10": int((n_cases >= 10).sum()),
        "unique_drugs": int(case_drugs_df["mapped_drug_id"].nunique()),
        "unique_drug_pairs": int(counts_df[["drug_a", "drug_b"]].drop_duplicates().shape[0]),
        "unique_reaction_pts": int(counts_df["reaction_pt"].nunique()),
        "duplicate_triple_count": int(
            counts_df.duplicated(subset=["year", "drug_a", "drug_b", "reaction_pt"]).sum()
        ),
        "null_drug_a_count": int(counts_df["drug_a"].isna().sum()),
        "null_drug_b_count": int(counts_df["drug_b"].isna().sum()),
        "null_reaction_pt_count": int(counts_df["reaction_pt"].isna().sum()),
        "top_20_reactions_by_n_triples": [
            {"reaction_pt": pt, "n_triples": int(n)}
            for pt, n in counts_df["reaction_pt"].value_counts().head(20).items()
        ],
        "runtime_seconds": round(runtime_seconds, 2),
        "output_file_size_mb": round(_dir_size_mb(year_out_dir), 2),
    }


# ---------------------------------------------------------------------------
# Process one year (mirrors process_faers_case_level.main()'s Milestone
# 2-6 sequence, called on year-concatenated data instead of one quarter)
# ---------------------------------------------------------------------------

def process_year(
    year: int,
    out_dir: str,
    mapping_file: str,
    drug_id_mode: str,
    pair_role_filter: str,
    max_rows_per_join_chunk: int,
    max_rows_per_quarter: int = None,
    include_fuzzy: bool = False,
) -> dict:
    t0 = time.time()
    year_out_dir = os.path.join(out_dir, str(year))
    os.makedirs(year_out_dir, exist_ok=True)

    print("=" * 60)
    print(f"FAERS CASE-LEVEL Annual Pipeline (POC) -- year {year}")
    print("=" * 60)
    print(f"  Drug ID mode:     {drug_id_mode}")
    print(f"  Pair role filter: {pair_role_filter}")
    print(f"  Output dir:       {year_out_dir}")
    if max_rows_per_quarter:
        print(f"  Max rows/quarter: {max_rows_per_quarter:,}  (TEST MODE)")

    print(f"\n--- {year}: load + concatenate DEMO/DRUG/REAC across Q1-Q4 ---")
    demo_df, drug_df, reac_df = load_year_raw_tables(year, max_rows_per_quarter)

    print(f"\n--- {year} Milestone 2: raw-table diagnostics ---")
    raw_diagnostics = pfcl.compute_raw_diagnostics(demo_df, drug_df, reac_df)
    pfcl.print_raw_diagnostics(raw_diagnostics)
    pfcl.save_raw_tables(year_out_dir, demo_df, drug_df, reac_df)
    pfcl.save_raw_diagnostics(year_out_dir, raw_diagnostics)

    print(f"\n--- {year} Milestone 3: resolve CASEID/CASEVERSION (year-level) ---")
    retained_cases_df = pfcl.resolve_case_versions(demo_df)
    resolved_drug_df = pfcl.filter_to_retained_primaryids(drug_df, retained_cases_df, "DRUG")
    resolved_reac_df = pfcl.filter_to_retained_primaryids(reac_df, retained_cases_df, "REAC")
    case_resolution_diagnostics = pfcl.compute_case_resolution_diagnostics(
        demo_df, retained_cases_df, drug_df, resolved_drug_df, reac_df, resolved_reac_df,
    )
    pfcl.print_case_resolution_diagnostics(case_resolution_diagnostics)
    pfcl.run_case_resolution_validation(demo_df, retained_cases_df, resolved_drug_df, resolved_reac_df)
    pfcl.save_case_resolution_artifacts(
        year_out_dir, retained_cases_df, resolved_drug_df, resolved_reac_df, case_resolution_diagnostics,
    )

    print(f"\n--- {year} Milestone 4: build case_reactions ---")
    case_reactions_df, reaction_dedup_diagnostics = pfcl.build_case_reactions(resolved_reac_df)
    pfcl.print_reaction_dedup_diagnostics(reaction_dedup_diagnostics)
    pfcl.run_case_reactions_validation(resolved_reac_df, case_reactions_df)
    pfcl.save_case_reactions_artifacts(year_out_dir, case_reactions_df, reaction_dedup_diagnostics)

    print(f"\n--- {year} Milestone 4c: build case_drugs (drug_id_mode={drug_id_mode}) ---")
    case_drugs_long_df, case_drugs_df, drug_mapping_diagnostics, unmapped_drugs_top_df = pfcl.build_case_drugs(
        resolved_drug_df, mapping_file, drug_id_mode=drug_id_mode, include_fuzzy=include_fuzzy,
    )
    pfcl.print_drug_mapping_diagnostics(drug_mapping_diagnostics)
    pfcl.run_case_drugs_validation(resolved_drug_df, case_drugs_long_df, case_drugs_df, drug_id_mode)
    pfcl.save_case_drugs_artifacts(
        year_out_dir, case_drugs_long_df, case_drugs_df, drug_mapping_diagnostics, unmapped_drugs_top_df,
    )

    print(f"\n--- {year} Milestone 5: build case_drug_pairs ---")
    case_drug_pairs_df, drug_pair_diagnostics, high_polypharmacy_df = pfcl.build_case_drug_pairs(
        case_drugs_df, drug_id_mode,
    )
    pfcl.print_drug_pair_diagnostics(drug_pair_diagnostics)
    pfcl.run_case_drug_pairs_validation(case_drugs_df, case_drug_pairs_df, drug_id_mode)
    pfcl.save_case_drug_pairs_artifacts(
        year_out_dir, case_drug_pairs_df, drug_pair_diagnostics, high_polypharmacy_df,
    )

    print(f"\n--- {year} Milestone 5c: role-aware pair audit ---")
    role_audit_df, role_audit_diagnostics = pfcl.build_case_drug_pair_role_audit(
        case_drug_pairs_df, case_drugs_df, drug_id_mode,
    )
    pfcl.print_drug_pair_role_audit_diagnostics(role_audit_diagnostics)
    pfcl.run_drug_pair_role_audit_validation(case_drug_pairs_df, role_audit_df, role_audit_diagnostics)
    pfcl.save_drug_pair_role_audit_artifacts(year_out_dir, role_audit_df, role_audit_diagnostics)

    print(f"\n--- {year} Milestone 5d: apply pair-role filter (pair_role_filter={pair_role_filter}) ---")
    selected_pairs_df, pair_filter_diagnostics = pfcl.apply_pair_role_filter(
        case_drug_pairs_df, role_audit_df, pair_role_filter, drug_id_mode,
    )
    pfcl.print_pair_role_filter_diagnostics(pair_filter_diagnostics)
    pfcl.run_pair_role_filter_validation(case_drug_pairs_df, role_audit_df, selected_pairs_df, pair_role_filter)
    pfcl.save_pair_role_filter_artifacts(year_out_dir, selected_pairs_df, pair_filter_diagnostics)

    print(f"\n--- {year} Milestone 6: build case_pair_reactions (co-reporting only, chunked) ---")
    case_pair_reactions_path = os.path.join(year_out_dir, "case_pair_reactions.parquet")
    n_pair_reaction_rows, drug_pair_reaction_counts_df, pair_reaction_case_stats_df, row_validation_summary = (
        pfcl.build_case_pair_reactions_chunked(
            selected_pairs_df, case_reactions_df, case_pair_reactions_path,
            drug_id_mode, pair_role_filter, max_rows_per_chunk=max_rows_per_join_chunk,
        )
    )
    pfcl.run_case_pair_reactions_validation(row_validation_summary, pair_role_filter, drug_id_mode)

    pair_reaction_diagnostics = pfcl.compute_pair_reaction_diagnostics(
        selected_pairs_df, case_reactions_df, n_pair_reaction_rows, drug_pair_reaction_counts_df,
        pair_reaction_case_stats_df, drug_id_mode, pair_role_filter,
    )
    pfcl.print_pair_reaction_diagnostics(pair_reaction_diagnostics)
    pfcl.run_drug_pair_reaction_counts_validation(drug_pair_reaction_counts_df, case_pair_reactions_path)

    # Add the year column AFTER the base script's own validation (which
    # checks the unqualified (drug_a, drug_b, reaction_pt) key) and BEFORE
    # saving, so it lands in drug_pair_reaction_counts.parquet and every
    # min3/min5/min10 variant written by save_pair_reaction_artifacts.
    drug_pair_reaction_counts_df = drug_pair_reaction_counts_df.copy()
    drug_pair_reaction_counts_df.insert(0, "year", year)

    pfcl.save_pair_reaction_artifacts(year_out_dir, drug_pair_reaction_counts_df, pair_reaction_diagnostics)
    pfcl.run_validation_checks(demo_df, retained_cases_df, case_reactions_df, case_drugs_df)

    _validate_annual_counts(drug_pair_reaction_counts_df, year)

    runtime_seconds = time.time() - t0
    summary = _build_annual_summary(
        year, year_out_dir, runtime_seconds,
        raw_diagnostics, case_resolution_diagnostics, reaction_dedup_diagnostics,
        case_drugs_df, case_drug_pairs_df, selected_pairs_df,
        n_pair_reaction_rows, drug_pair_reaction_counts_df,
    )
    summary_path = os.path.join(year_out_dir, "annual_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Saved: {summary_path}")
    print(f"  Year {year} complete in {runtime_seconds:.1f}s. Artifacts in: {year_out_dir}")

    return summary


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Annual FAERS case-level wrapper (POC). Concatenates each "
            "requested year's four quarterly DEMO/DRUG/REAC tables and runs "
            "process_faers_case_level.py's existing pipeline functions "
            "(unmodified) once per year, so CASEID/CASEVERSION resolution "
            "isn't double-counted across quarter boundaries. See file "
            "docstring for details and scope."
        )
    )
    parser.add_argument(
        "--years", type=int, nargs="+", required=True,
        help="One or more years to process, e.g. --years 2023 2024 2025. "
             "Each year requires all four quarterly zips to be present.",
    )
    parser.add_argument(
        "--out-dir", default=DEFAULT_OUT_DIR,
        help=f"Output directory; a subdirectory per year is created under it. Default: {DEFAULT_OUT_DIR}",
    )
    parser.add_argument(
        "--mapping-file", default=DEFAULT_MAPPING_FILE,
        help=f"Path to the FAERS drug name -> ChEMBL ID mapping parquet. Default: {DEFAULT_MAPPING_FILE}",
    )
    parser.add_argument(
        "--include-fuzzy", action="store_true",
        help="Include fuzzy-matched drug names when building case_drugs. Default: excluded.",
    )
    parser.add_argument(
        "--drug-id-mode", choices=pfcl.DRUG_ID_MODES, default=pfcl.DEFAULT_DRUG_ID_MODE,
        help="How mapped_drug_id is derived -- see process_faers_case_level.py --drug-id-mode.",
    )
    parser.add_argument(
        "--pair-role-filter", choices=pfcl.PAIR_ROLE_FILTERS, default=pfcl.DEFAULT_PAIR_ROLE_FILTER,
        help="Which case_drug_pairs rows proceed to the pair-reaction join -- see "
             "process_faers_case_level.py --pair-role-filter.",
    )
    parser.add_argument(
        "--max-rows-per-join-chunk", type=int, default=pfcl.DEFAULT_MAX_ROWS_PER_JOIN_CHUNK,
        help="Caps case_pair_reactions join batch size -- see process_faers_case_level.py of the same name.",
    )
    parser.add_argument(
        "--max-rows-per-quarter", type=int, default=None,
        help="Smoke-test only: cap DEMO/DRUG/REAC to the first N rows PER QUARTER "
             "before concatenation (so a year's total row cap is ~4x this value).",
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    all_summaries = []
    for year in args.years:
        summary = process_year(
            year=year,
            out_dir=args.out_dir,
            mapping_file=args.mapping_file,
            drug_id_mode=args.drug_id_mode,
            pair_role_filter=args.pair_role_filter,
            max_rows_per_join_chunk=args.max_rows_per_join_chunk,
            max_rows_per_quarter=args.max_rows_per_quarter,
            include_fuzzy=args.include_fuzzy,
        )
        all_summaries.append(summary)

    combined_path = os.path.join(args.out_dir, "annual_summary_all_years.json")
    with open(combined_path, "w", encoding="utf-8") as f:
        json.dump(all_summaries, f, indent=2)
    print(f"\nSaved combined summary: {combined_path}")


if __name__ == "__main__":
    main()
