"""
Pipeline orchestrator — Tiers 0–2.

Each tier has the same contract:
    tier_N(unresolved_df, cache, cfg) -> (newly_resolved_df, still_unresolved_df)

The top-level run_tiers_0_2() assembles them into a cascade and returns:
    resolved_df    — all strings that received a RxCUI
    unresolved_df  — strings still without a match (passed to Tier 3+)

Output schema (both DataFrames share these columns):
    raw_string      original FAERS value
    cleaned_string  post-Tier-0 form
    rxcui           RxNorm concept ID (None for unresolved)
    chembl_id       filled by UniChem in a later step (always None here)
    match_tier      int (1 or 2) for resolved; None for unresolved
    match_score     float in [0, 1] (0.0 for unresolved)
    needs_review    bool
    notes           free text ("combo", "salt_stripped", etc.)
"""
import time
from typing import Tuple

import pandas as pd

from .cache import ResponseCache
from .config import CONFIG
from .preprocess import preprocess
from .report import tier_report, final_report
from .rxnorm import rxcui_exact, rxcui_approximate, rxcui_to_ingredient

_SCHEMA = [
    "raw_string", "cleaned_string", "rxcui", "chembl_id",
    "match_tier", "match_score", "needs_review", "notes",
]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _empty_resolved() -> pd.DataFrame:
    return pd.DataFrame(columns=_SCHEMA)


def _make_row(preprocessed: dict, rxcui: str, score: float, tier: int, note: str = "") -> dict:
    return {
        "raw_string":     preprocessed["raw_string"],
        "cleaned_string": preprocessed["cleaned_string"],
        "rxcui":          rxcui,
        "chembl_id":      None,
        "match_tier":     tier,
        "match_score":    round(score, 4),
        "needs_review":   score < CONFIG["needs_review_threshold"],
        "notes":          note,
    }


def _unresolved_frame(work_df: pd.DataFrame) -> pd.DataFrame:
    """Attach the output schema columns to still-unresolved rows."""
    out = work_df[["raw_string", "cleaned_string"]].copy()
    out["rxcui"]        = None
    out["chembl_id"]    = None
    out["match_tier"]   = None
    out["match_score"]  = 0.0
    out["needs_review"] = True
    out["notes"]        = out.get("notes", "")
    return out[_SCHEMA]


# ---------------------------------------------------------------------------
# Tier 1 — RxNorm exact
# ---------------------------------------------------------------------------

def _tier1(
    unresolved: pd.DataFrame,
    cache: ResponseCache,
    cfg: dict,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    delay    = cfg["rxnav_delay_s"]
    every    = cfg["report_every"]
    n        = len(unresolved)
    resolved = []

    for i, row in enumerate(unresolved.itertuples(index=False), 1):
        if i % every == 0:
            print(f"  Tier 1 {i:,}/{n:,} ...", end="\r")

        p = row._asdict()
        note = "combo" if p["is_combo"] else ""

        # Try cleaned form first
        rxcui = rxcui_exact(p["cleaned_string"], cache, delay)
        if rxcui is None and isinstance(p["salt_stripped"], str):
            rxcui = rxcui_exact(p["salt_stripped"], cache, delay)
            if rxcui:
                note = (note + " salt_stripped").strip()

        if rxcui:
            rxcui = rxcui_to_ingredient(rxcui, cache, delay)
            resolved.append(_make_row(p, rxcui, 1.0, tier=1, note=note))

    print()
    resolved_df    = pd.DataFrame(resolved) if resolved else _empty_resolved()
    resolved_clean = set(resolved_df["cleaned_string"]) if resolved else set()
    remaining      = unresolved[~unresolved["cleaned_string"].isin(resolved_clean)].copy()
    return resolved_df, remaining


# ---------------------------------------------------------------------------
# Tier 2 — RxNorm approximate
# ---------------------------------------------------------------------------

def _tier2(
    unresolved: pd.DataFrame,
    cache: ResponseCache,
    cfg: dict,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    delay      = cfg["rxnav_delay_s"]
    every      = cfg["report_every"]
    max_ent    = cfg["tier2_max_entries"]
    min_score  = cfg["tier2_min_score"]
    n          = len(unresolved)
    resolved   = []

    for i, row in enumerate(unresolved.itertuples(index=False), 1):
        if i % every == 0:
            print(f"  Tier 2 {i:,}/{n:,} ...", end="\r")

        p    = row._asdict()
        note = "combo" if p["is_combo"] else ""

        rxcui, score = rxcui_approximate(
            p["cleaned_string"], cache, max_ent, min_score, delay
        )
        if rxcui is None and isinstance(p["salt_stripped"], str):
            rxcui, score = rxcui_approximate(
                p["salt_stripped"], cache, max_ent, min_score, delay
            )
            if rxcui:
                note = (note + " salt_stripped").strip()

        if rxcui:
            rxcui = rxcui_to_ingredient(rxcui, cache, delay)
            resolved.append(_make_row(p, rxcui, score, tier=2, note=note))

    print()
    resolved_df    = pd.DataFrame(resolved) if resolved else _empty_resolved()
    resolved_clean = set(resolved_df["cleaned_string"]) if resolved else set()
    remaining      = unresolved[~unresolved["cleaned_string"].isin(resolved_clean)].copy()
    return resolved_df, remaining


# ---------------------------------------------------------------------------
# Public entrypoint
# ---------------------------------------------------------------------------

def run_tiers_0_2(
    raw_series: pd.Series,
    cfg: dict = None,
    cache_db: str = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run Tiers 0–2 on a Series of raw FAERS drug strings.

    Parameters
    ----------
    raw_series  : Series of raw strings (already field-selected, e.g. prod_ai or drugname)
    cfg         : override CONFIG dict (useful for testing)
    cache_db    : override cache database path

    Returns
    -------
    resolved_df    : DataFrame with rxcui assigned, schema = _SCHEMA
    unresolved_df  : remaining strings, rxcui=None — passed to Tier 3+
    """
    cfg      = cfg or CONFIG
    db_path  = cache_db or cfg["cache_db"]
    cache    = ResponseCache(db_path)

    # ── Tier 0: preprocess all unique raw strings ────────────────────────────
    unique_raw = [s for s in raw_series.dropna().unique() if str(s).strip()]
    print(f"\n[Tier 0] Preprocessing {len(unique_raw):,} unique raw strings ...")
    t0 = time.perf_counter()

    preprocessed = [preprocess(r) for r in unique_raw]

    # Deduplicate on cleaned form (many raw strings collapse to the same form)
    work_df = (
        pd.DataFrame(preprocessed)
        .drop_duplicates(subset=["cleaned_string"])
        .query("cleaned_string != ''")
        .reset_index(drop=True)
    )
    print(f"  Unique cleaned strings: {len(work_df):,}")
    tier_report(0, n_total=len(work_df), n_this_tier=0, n_cumulative=0,
                elapsed_s=time.perf_counter() - t0)

    n_total      = len(work_df)
    all_resolved = []

    # ── Tier 1 ───────────────────────────────────────────────────────────────
    print(f"\n[Tier 1] RxNorm exact  ({len(work_df):,} strings) ...")
    t1 = time.perf_counter()
    res1, work_df = _tier1(work_df, cache, cfg)
    all_resolved.append(res1)
    mean_score1 = res1["match_score"].mean() if len(res1) else None
    tier_report(1, n_total, len(res1), len(res1),
                mean_score=mean_score1, elapsed_s=time.perf_counter() - t1)

    # ── Tier 2 ───────────────────────────────────────────────────────────────
    print(f"\n[Tier 2] RxNorm approx ({len(work_df):,} strings) ...")
    t2 = time.perf_counter()
    res2, work_df = _tier2(work_df, cache, cfg)
    all_resolved.append(res2)
    n_cum = len(res1) + len(res2)
    mean_score2 = res2["match_score"].mean() if len(res2) else None
    tier_report(2, n_total, len(res2), n_cum,
                mean_score=mean_score2, elapsed_s=time.perf_counter() - t2)

    # ── Assemble outputs ─────────────────────────────────────────────────────
    resolved_df   = pd.concat([df for df in all_resolved if len(df)], ignore_index=True)
    unresolved_df = _unresolved_frame(work_df)

    cs = cache.stats()
    print(f"\n  Cache: {cs['entries']:,} entries  "
          f"hits={cs['cache_hits']:,}  misses={cs['cache_misses']:,}  "
          f"hit_rate={cs['hit_rate']:.1%}")
    cache.close()

    final_report(resolved_df, unresolved_df)
    return resolved_df, unresolved_df
