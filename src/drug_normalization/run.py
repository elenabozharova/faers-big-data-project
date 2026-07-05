"""
CLI entry point — run Tiers 0–2 on a FAERS quarterly ZIP and save results.

Usage
-----
  # From the project root:
  python src/drug_normalization/run.py

  # Override input / output:
  python src/drug_normalization/run.py \\
      --faers-zip  data/raw/faers/faers_ascii_2025Q4.zip \\
      --drug-file  ASCII/DRUG25Q4.txt \\
      --output     data/processed/drug_mapping_rxnorm.parquet

  # Quick smoke test on first 10 000 rows:
  python src/drug_normalization/run.py --max-rows 10000

  # Run test set validation after the cascade:
  python src/drug_normalization/run.py --max-rows 10000 --run-tests

Field selection (Tier 0, Step 0)
---------------------------------
Mirrors the logic in map_faers_to_chembl.py:
  Use prod_ai when present, non-empty (>2 chars), and not purely numeric;
  otherwise fall back to drugname.

Output
------
  --output (parquet)     resolved + unresolved rows, full schema
  Unresolved strings are included with rxcui=None so nothing is silently dropped.
"""
import argparse
import os
import sys
import zipfile

import pandas as pd

# Make src/ importable when running as a script
_SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import config as _prj_cfg
from drug_normalization.pipeline import run_tiers_0_2
from drug_normalization.config import CONFIG


# ---------------------------------------------------------------------------
# Field selection
# ---------------------------------------------------------------------------

def select_drug_field(drug_df: pd.DataFrame) -> pd.Series:
    """
    Return a Series of raw drug strings using prod_ai where meaningful,
    otherwise drugname.  Mirrors the Tier-0 logic in map_faers_to_chembl.py.
    """
    prod_ai  = (
        drug_df["prod_ai"].fillna("").astype(str).str.strip()
        if "prod_ai" in drug_df.columns
        else pd.Series([""] * len(drug_df))
    )
    drugname = (
        drug_df["drugname"].fillna("").astype(str).str.strip()
        if "drugname" in drug_df.columns
        else pd.Series([""] * len(drug_df))
    )
    meaningful = (prod_ai.str.len() > 2) & (
        ~prod_ai.str.replace(".", "", regex=False).str.isnumeric()
    )
    return prod_ai.where(meaningful, drugname)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Map FAERS drug strings to RxCUI (Tiers 0–2)."
    )
    parser.add_argument(
        "--faers-zip",
        default=_prj_cfg.FAERS_ZIP_PATH,
        help="Path to FAERS quarterly ZIP (default: from config).",
    )
    parser.add_argument(
        "--drug-file",
        default=_prj_cfg.FAERS_DRUG_INNER,
        help="Inner path to DRUG table inside the ZIP.",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(
            os.path.dirname(_SRC), "data", "processed", "drug_mapping_rxnorm.parquet"
        ),
        help="Output parquet path.",
    )
    parser.add_argument(
        "--max-rows", type=int, default=None,
        help="Limit DRUG rows for a quick test run.",
    )
    parser.add_argument(
        "--run-tests", action="store_true",
        help="Run the built-in test set after the cascade and report precision.",
    )
    parser.add_argument(
        "--tier2-min-score", type=int, default=CONFIG["tier2_min_score"],
        help=f"Minimum RxNav approximate score (default: {CONFIG['tier2_min_score']}).",
    )
    args = parser.parse_args()

    # Allow per-run threshold override without editing config
    cfg = dict(CONFIG)
    cfg["tier2_min_score"] = args.tier2_min_score

    # ── Load FAERS DRUG table ─────────────────────────────────────────────
    print(f"Loading FAERS: {args.faers_zip} / {args.drug_file}")
    with zipfile.ZipFile(args.faers_zip) as z:
        with z.open(args.drug_file) as f:
            drug_df = pd.read_csv(f, sep="$", dtype=str, low_memory=False,
                                  nrows=args.max_rows)
    drug_df.columns = drug_df.columns.str.lower().str.strip()
    print(f"  Loaded {len(drug_df):,} DRUG rows")

    raw_series = select_drug_field(drug_df)

    # ── Run pipeline ──────────────────────────────────────────────────────
    resolved_df, unresolved_df = run_tiers_0_2(raw_series, cfg=cfg)

    # ── Save ─────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    combined = pd.concat([resolved_df, unresolved_df], ignore_index=True)
    combined.to_parquet(args.output, index=False)
    print(f"\nSaved {len(combined):,} rows -> {args.output}")

    # ── Optional test set ─────────────────────────────────────────────────
    if args.run_tests:
        from drug_normalization.test_cases import run_test_report
        run_test_report(resolved_df, unresolved_df)


if __name__ == "__main__":
    main()
