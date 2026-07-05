"""
Tier 0 — Preprocessing.

Takes a raw FAERS drug string and returns a structured dict with:
  - cleaned_string  : lowercased, dosage/form tokens removed
  - salt_stripped   : salt/ester suffix removed (None if no change)
  - is_combo        : True when the string encodes a combination product
  - components      : individual component strings for combos

No API calls.  Every transformation is deterministic and side-effect free.

Design note: dosage stripping runs BEFORE combo detection so that
"ibuprofen 200 mg / paracetamol 500 mg" correctly splits into two
components rather than treating "200 mg / paracetamol" as a combo boundary.
"""
import re
import unicodedata
from typing import List, Optional


# ---------------------------------------------------------------------------
# Compiled patterns
# ---------------------------------------------------------------------------

# Dosage amounts: "500 mg", "10mcg", "100 iu", "0.5%", etc.
_DOSAGE_RE = re.compile(
    r"\b\d+(?:\.\d+)?\s*"
    r"(?:mg|mcg|ug|g|ml|l|iu|mmol|meq|mEq|units?|%)"
    r"\b",
    re.IGNORECASE,
)

# Dosage-form and route tokens that add no identity information
_FORM_RE = re.compile(
    r"\b(?:"
    r"tablet|tablets|tab|tabs|"
    r"capsule|capsules|cap|caps|"
    r"injection|inj|injectable|"
    r"solution|sol|"
    r"suspension|syrup|"
    r"cream|ointment|lotion|gel|"
    r"patch|transdermal|"
    r"spray|inhaler|inhalation|"
    r"er|xr|sr|cr|ir|la|xl|cd|dr|"    # extended / sustained release
    r"oral|intravenous|iv|im|sc|"
    r"hcl"                             # salt abbreviation also stripped here
    r")\b",
    re.IGNORECASE,
)

# Separators that mark a combination product
# Handles: "/", "\", ";", "+", "and", "with"
# Anchored with \b to avoid splitting on partial words.
_COMBO_SEP_RE = re.compile(
    r"[/;+\\]|\band\b|\bwith\b",
    re.IGNORECASE,
)

# Salt / ester / hydration suffixes to try stripping
# These appear at the end of the cleaned string as a separate word.
_SALT_RE = re.compile(
    r"\s+(?:"
    r"hydrochloride|hcl|"
    r"sulfate|sulphate|"
    r"sodium|potassium|calcium|magnesium|"
    r"citrate|maleate|fumarate|mesylate|tosylate|"
    r"tartrate|acetate|phosphate|bromide|chloride|"
    r"gluconate|succinate|besylate|benzoate|"
    r"monohydrate|dihydrate|trihydrate|anhydrous|"
    r"hydrochloride|hydrobromide|nitrate|oxalate"
    r")\s*$",
    re.IGNORECASE,
)

_MULTI_SPACE_RE = re.compile(r"\s+")


# ---------------------------------------------------------------------------
# Individual transformation steps (exported for unit testing)
# ---------------------------------------------------------------------------

def unicode_clean(s: str) -> str:
    """NFKD-normalize and drop combining characters (accents, diacritics)."""
    nfkd = unicodedata.normalize("NFKD", s)
    return "".join(c for c in nfkd if not unicodedata.combining(c))


def remove_dosages(s: str) -> str:
    return _DOSAGE_RE.sub(" ", s)


def remove_form_tokens(s: str) -> str:
    return _FORM_RE.sub(" ", s)


def collapse_whitespace(s: str) -> str:
    return _MULTI_SPACE_RE.sub(" ", s).strip()


def clean(raw: str) -> str:
    """
    Full Tier-0 cleaning in one call.
    Does NOT strip salt suffixes — that is a separate variant.
    """
    s = unicode_clean(str(raw) if raw is not None else "")
    s = s.lower().strip()
    s = remove_dosages(s)
    s = remove_form_tokens(s)
    s = collapse_whitespace(s)
    return s


def strip_salts(cleaned: str) -> Optional[str]:
    """
    Return a salt-stripped variant, or None if the string does not end
    with a known suffix.  Both variants are tried in Tiers 1–2.
    """
    stripped = _SALT_RE.sub("", cleaned).strip()
    return stripped if stripped != cleaned and stripped else None


def is_combo(cleaned: str) -> bool:
    """Return True when the cleaned string encodes a combination product."""
    return bool(_COMBO_SEP_RE.search(cleaned))


def split_components(cleaned: str) -> List[str]:
    """
    Split a combination string into its individual components.
    Each component is itself cleaned but not recursively split.
    """
    return [p.strip() for p in _COMBO_SEP_RE.split(cleaned) if p.strip()]


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def preprocess(raw: str) -> dict:
    """
    Full Tier-0 preprocessing for one raw FAERS drug string.

    Returns
    -------
    {
      "raw_string":    original value,
      "cleaned_string": cleaned form (empty string if input is garbage),
      "salt_stripped":  salt-stripped variant or None,
      "is_combo":       bool,
      "components":     list[str] — individual component strings for combos,
    }
    """
    raw = str(raw).strip() if raw is not None else ""
    cleaned = clean(raw)
    salt = strip_salts(cleaned)
    combo = is_combo(cleaned)
    components = split_components(cleaned) if combo else []

    return {
        "raw_string":     raw,
        "cleaned_string": cleaned,
        "salt_stripped":  salt,
        "is_combo":       combo,
        "components":     components,
    }
