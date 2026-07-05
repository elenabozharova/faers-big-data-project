"""
Drug name normalization utilities for FAERS -> ChEMBL mapping.

Public API
----------
normalize_drug_string(s) -> str
    Converts a raw FAERS drug name to a canonical normalized form.
    Combo drugs are represented with "|" as a segment separator.

split_combo(normalized) -> List[str]
    Splits a normalized drug name into its constituent components.
    Returns a single-element list for simple drugs.

is_combo(normalized) -> bool
    Returns True if the normalized name contains a "|" separator.

Design
------
Normalization steps (applied in order):
  1. Lowercase
  2. Remove dosage patterns (numbers + optional units) BEFORE separator
     substitution so that "10/500 mg" is treated as a dosage ratio, not
     a combo separator.
  3. Replace combo separators (/ + and with plus) with "|".
  4. Strip remaining punctuation (preserve "|" and word-internal hyphens).
  5. Per "|"-segment: remove standalone hyphens, then remove noise tokens
     (dosage forms, routes, release modifiers, pack descriptors).
  6. Drop empty segments; join with "|".

Examples
--------
  "aspirin 81mg"                        -> "aspirin"
  "AMOXICILLIN/CLAVULANATE 875-125 MG"  -> "amoxicillin|clavulanate"
  "METOPROLOL SUCCINATE ER 50 MG TAB"   -> "metoprolol succinate"
  "METFORMIN HCL 500 MG TABLET"         -> "metformin hcl"
  "lisinopril 10 mg oral tablet"         -> "lisinopril"
  "ASPIRIN + DIPYRIDAMOLE"               -> "aspirin|dipyridamole"
"""

import re
from typing import List


# ---------------------------------------------------------------------------
# Noise tokens: dosage forms, routes of administration, release modifiers,
# pack/presentation descriptors.
# NOTE: chemical salt descriptors (hydrochloride, succinate, sodium, …) are
# intentionally NOT in this list — they are part of the drug name as ChEMBL
# records it, so removing them would hurt synonym matching.
# ---------------------------------------------------------------------------
_NOISE_TOKENS: frozenset = frozenset({
    # Measurement units (residual after dosage-pattern removal)
    "mg", "mcg", "ug", "ml", "g", "l", "iu", "meq",
    # Dosage forms
    "tab", "tablet", "tablets", "caps", "cap", "capsule", "capsules",
    "inj", "injection", "injectable", "infusion",
    "solution", "suspension", "syrup", "elixir", "emulsion",
    "patch", "cream", "ointment", "gel", "lotion", "foam",
    "spray", "drops", "inhaler", "aerosol", "powder",
    "granule", "granules", "pellet", "pellets",
    "suppository", "suppositories", "enema",
    "lozenge", "lozenges", "pastille",
    # Routes of administration
    "iv", "po", "sc", "im", "sq",
    "oral", "topical", "sublingual", "transdermal",
    "intravenous", "subcutaneous", "intramuscular",
    "ophthalmic", "otic", "nasal", "rectal", "vaginal",
    # Release / formulation modifiers
    "sr", "xr", "er", "dr", "cr", "la", "xl", "ir",
    "extended", "release", "delayed", "controlled",
    "modified", "immediate", "sustained", "slow",
    # Pack / presentation noise
    "film", "coated", "scored", "chewable", "effervescent",
    "pack", "dose", "vial", "ampule", "ampoule", "kit",
    # Common abbreviations that add no identity information
    "nos",   # Not Otherwise Specified
    "nec",   # Not Elsewhere Classified
    "hcpcs",
})

# ---------------------------------------------------------------------------
# Compiled regex patterns
# ---------------------------------------------------------------------------

# Combo separators -> "|"
_SEP_RE = re.compile(
    r"[/+]"
    r"|\band\b"
    r"|\bwith\b"
    r"|\bplus\b"
    r"|\bwithout\b",
    re.IGNORECASE,
)

# Dosage numeric patterns with optional unit.
# Applied BEFORE separator substitution so "10/500 mg" (a ratio) is removed
# rather than converted to a combo separator.
# Matches: 81, 875-125, 0.5, 10/500, 20mg, 1.5mcg, 50%, etc.
_DOSAGE_RE = re.compile(
    r"\b\d+(?:[.,/-]\d+)*\s*(?:mg|mcg|ug|ml|g|l|iu|meq|%)?\b",
    re.IGNORECASE,
)

# Punctuation to remove (preserve "|" which we use as a combo marker,
# and "-" which may be part of drug names like "anti-infective").
_PUNCT_RE = re.compile(r"[^\w\s|]")

# Standalone hyphen: a "-" that is NOT surrounded by word characters.
# Removes artifacts like " - " or leading/trailing hyphens inside segments.
_STANDALONE_HYPHEN_RE = re.compile(r"(?<!\w)-|-(?!\w)")

# Collapse runs of whitespace
_WS_RE = re.compile(r"\s+")


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

def normalize_drug_string(s: str) -> str:
    """
    Normalize a raw FAERS drug name to a canonical, lowercase form.

    Combo drugs yield "|"-separated segments:
        "AMOXICILLIN/CLAVULANATE 875-125 MG TAB" -> "amoxicillin|clavulanate"
        "ASPIRIN + DIPYRIDAMOLE 25-200MG"         -> "aspirin|dipyridamole"

    Simple drugs yield a single segment (no "|"):
        "aspirin 81mg"                  -> "aspirin"
        "METOPROLOL SUCCINATE ER 50 MG" -> "metoprolol succinate"
        "METFORMIN HCL 500 MG TABLET"   -> "metformin hcl"

    Empty or non-string inputs return "".
    """
    if not s or not isinstance(s, str):
        return ""

    s = s.lower().strip()

    # Step 1: remove dosage patterns BEFORE converting separators.
    # This ensures "10/500 mg" is deleted rather than turned into "10|500".
    s = _DOSAGE_RE.sub(" ", s)

    # Step 2: replace combo separators with "|".
    s = _SEP_RE.sub("|", s)

    # Step 3: strip punctuation, keeping "|" and "-".
    s = _PUNCT_RE.sub(" ", s)

    # Step 4: process each "|"-delimited segment independently.
    segments: List[str] = []
    for seg in s.split("|"):
        # Remove standalone hyphens (e.g. " - " artifacts from punctuation removal).
        seg = _STANDALONE_HYPHEN_RE.sub(" ", seg)
        # Tokenize and drop noise tokens and bare digits.
        tokens = _WS_RE.split(seg.strip())
        tokens = [t for t in tokens if t and t not in _NOISE_TOKENS and not t.isdigit()]
        seg_clean = " ".join(tokens).strip()
        if seg_clean:
            segments.append(seg_clean)

    return "|".join(segments)


def split_combo(normalized: str) -> List[str]:
    """
    Split a normalized drug name into its constituent component strings.

    For simple drugs, returns a single-element list.
    For combo drugs (containing "|"), returns each component as a string.
    Returns an empty list for empty input.

    Args:
        normalized: Output of normalize_drug_string().

    Examples:
        split_combo("aspirin")                    -> ["aspirin"]
        split_combo("amoxicillin|clavulanate")    -> ["amoxicillin", "clavulanate"]
        split_combo("aspirin|dipyridamole")        -> ["aspirin", "dipyridamole"]
        split_combo("")                            -> []
    """
    if not normalized:
        return []
    return [p.strip() for p in normalized.split("|") if p.strip()]


def is_combo(normalized: str) -> bool:
    """Return True if the normalized name represents a combination drug."""
    return "|" in normalized and len(split_combo(normalized)) > 1
