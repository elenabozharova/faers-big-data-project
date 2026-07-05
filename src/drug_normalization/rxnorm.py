"""
Tiers 1–2 — RxNorm lookups via RxNav REST API.

Tier 1  exact name match     /rxcui.json?name={string}&search=1
Tier 2  approximate match    /approximateTerm.json?term={string}&maxEntries=N

Every request goes through ResponseCache so re-runs never call the network
for a URL that was already fetched.  Sleep is inserted only on cache misses.

RxNav score note
----------------
/approximateTerm returns an internal integer score that reflects edit-distance
similarity.  Values typically range 0–100, but short exact-ish matches can
score above 100.  We clip at 100 and store as a float in [0, 1].
"""
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Optional, Tuple

from .cache import ResponseCache

_BASE = "https://rxnav.nlm.nih.gov/REST"


# ---------------------------------------------------------------------------
# Internal HTTP helper
# ---------------------------------------------------------------------------

def _fetch(url: str, cache: ResponseCache, delay: float) -> Optional[dict]:
    """
    Return the parsed JSON for url, hitting the cache first.
    Returns None on network or parse error (caller treats as no match).
    """
    cached = cache.get(url)
    if cached is not None:
        return cached

    time.sleep(delay)
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        cache.set(url, data)
        return data
    except (urllib.error.URLError, json.JSONDecodeError, OSError):
        return None


# ---------------------------------------------------------------------------
# Tier 1 — exact match
# ---------------------------------------------------------------------------

def rxcui_exact(
    name: str,
    cache: ResponseCache,
    delay: float = 0.06,
) -> Optional[str]:
    """
    Exact RxNorm concept lookup via /rxcui.json.

    Returns the first RxCUI string on an unambiguous match, else None.
    Note: search=1 enables approximate search within the endpoint; we still
    treat any returned RxCUI as a Tier-1 hit because the caller only sends
    cleaned strings here.
    """
    if not name or not isinstance(name, str):
        return None
    url = f"{_BASE}/rxcui.json?name={urllib.parse.quote(name, safe='')}&search=1"
    data = _fetch(url, cache, delay)
    if data is None:
        return None
    try:
        ids = data["idGroup"]["rxnormId"]
        if ids and ids[0]:
            return str(ids[0])
    except (KeyError, TypeError, IndexError):
        pass
    return None


# ---------------------------------------------------------------------------
# Ingredient normalization — walk up to tty=IN after any tier resolves
# ---------------------------------------------------------------------------

def rxcui_to_ingredient(
    rxcui: str,
    cache: ResponseCache,
    delay: float = 0.06,
) -> str:
    """
    Return the ingredient-level RxCUI for a given RxCUI.

    Queries /rxcui/{rxcui}/related.json?tty=IN.  If no IN-level concept is
    found (e.g. already an ingredient, or no mapping), returns rxcui unchanged.
    Results are cached so re-runs are free.
    """
    if not rxcui or not isinstance(rxcui, str):
        return rxcui
    url = f"{_BASE}/rxcui/{urllib.parse.quote(rxcui, safe='')}/related.json?tty=IN"
    data = _fetch(url, cache, delay)
    if data is None:
        return rxcui
    try:
        groups = data["relatedGroup"]["conceptGroup"]
        for group in groups:
            if group.get("tty") == "IN":
                concepts = group.get("conceptProperties") or []
                if concepts:
                    return str(concepts[0]["rxcui"])
    except (KeyError, TypeError, IndexError):
        pass
    return rxcui


# ---------------------------------------------------------------------------
# Tier 2 — approximate match
# ---------------------------------------------------------------------------

def rxcui_approximate(
    name: str,
    cache: ResponseCache,
    max_entries: int = 5,
    min_score: int = 75,
    delay: float = 0.06,
) -> Tuple[Optional[str], float]:
    """
    Approximate RxNorm match via /approximateTerm.json.

    Returns (rxcui, normalised_score) for the top candidate above min_score,
    or (None, 0.0) when nothing clears the threshold.
    """
    if not name or not isinstance(name, str):
        return None, 0.0

    url = (
        f"{_BASE}/approximateTerm.json"
        f"?term={urllib.parse.quote(name, safe='')}"
        f"&maxEntries={max_entries}"
    )
    data = _fetch(url, cache, delay)
    if data is None:
        return None, 0.0

    try:
        candidates = data["approximateGroup"]["candidate"]
        if not candidates:
            return None, 0.0
        best = candidates[0]
        raw_score = int(best.get("score", 0))
        if raw_score >= min_score:
            normalised = round(min(raw_score, 100) / 100.0, 4)
            return str(best["rxcui"]), normalised
    except (KeyError, TypeError, ValueError, IndexError):
        pass

    return None, 0.0
