"""
All tunable thresholds and paths in one place.
Edit here; do not scatter magic numbers through the pipeline code.
"""
import os

_SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ROOT = os.path.dirname(_SRC)

CONFIG: dict = {
    # ── Tier 2: RxNorm approximate match ─────────────────────────────────────
    # RxNav returns an integer score; higher = more similar.  Typical range
    # 0–100 (can exceed 100 for very short exact-ish matches).  75 is a
    # conservative starting point: captures close misspellings while avoiding
    # false positives on short strings.
    "tier2_min_score":   75,
    "tier2_max_entries":  5,    # candidates requested from RxNav

    # ── Tier 3: string similarity (placeholders — tuned after Tier 2 review) ─
    "tier3_tfidf_min":   0.85,  # sparse cosine (character n-gram TF-IDF)
    "tier3_fuzzy_min":   90,    # rapidfuzz token_set_ratio [0–100]
    "tier3_phonetic":    True,  # enable Double Metaphone fallback

    # ── General ───────────────────────────────────────────────────────────────
    # Matches below this score are flagged needs_review=True regardless of tier.
    "needs_review_threshold": 0.80,

    # ── RxNav rate limiting ───────────────────────────────────────────────────
    # RxNav asks for ≤ 20 req/s.  0.06 s gives comfortable headroom.
    "rxnav_delay_s": 0.06,

    # ── Progress reporting ────────────────────────────────────────────────────
    "report_every": 500,   # print a progress line every N strings

    # ── Cache ─────────────────────────────────────────────────────────────────
    "cache_db": os.path.join(_ROOT, "data", "cache", "rxnorm_cache.sqlite"),
}
