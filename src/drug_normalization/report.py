"""
Coverage reporting — printed after each tier completes.
"""
from typing import Optional


def tier_report(
    tier: int,
    n_total: int,
    n_this_tier: int,
    n_cumulative: int,
    mean_score: Optional[float] = None,
    elapsed_s: Optional[float] = None,
) -> None:
    """Print a standardised coverage line after a tier finishes."""
    pct_this = 100 * n_this_tier / n_total if n_total else 0.0
    pct_cum  = 100 * n_cumulative / n_total if n_total else 0.0
    remaining = n_total - n_cumulative

    bar = "=" * 60
    print(f"\n{bar}")
    print(f"  Tier {tier} complete")
    print(f"  Resolved this tier:  {n_this_tier:7,}  ({pct_this:5.1f}%)")
    print(f"  Cumulative:          {n_cumulative:7,} / {n_total:,}  ({pct_cum:5.1f}%)")
    print(f"  Remaining:           {remaining:7,}  ({100 - pct_cum:5.1f}%)")
    if mean_score is not None:
        print(f"  Mean score (tier):   {mean_score:.4f}")
    if elapsed_s is not None:
        print(f"  Elapsed:             {elapsed_s:.1f} s")
    print(bar)


def final_report(resolved_df, unresolved_df) -> None:
    """Print overall pipeline summary."""
    import pandas as pd

    total = len(resolved_df) + len(unresolved_df)
    n_res = len(resolved_df)

    print("\n" + "=" * 60)
    print("  PIPELINE SUMMARY")
    print("=" * 60)
    print(f"  Total unique cleaned strings:  {total:,}")
    print(f"  Resolved:                      {n_res:,}  ({100*n_res/total:.1f}%)")
    print(f"  Unresolved (needs manual/T3+): {len(unresolved_df):,}")

    if n_res:
        print(f"\n  Tier breakdown:")
        tier_counts = resolved_df["match_tier"].value_counts().sort_index()
        for tier, cnt in tier_counts.items():
            pct = 100 * cnt / n_res
            print(f"    Tier {tier}: {cnt:6,}  ({pct:.1f}% of resolved)")

        needs_review = resolved_df["needs_review"].sum()
        print(f"\n  Flagged needs_review: {needs_review:,}  ({100*needs_review/n_res:.1f}%)")

    print("=" * 60)
