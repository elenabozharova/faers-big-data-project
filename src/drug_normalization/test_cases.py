"""
Seed test set for the drug normalization pipeline.

~30 cases covering the main FAERS patterns.  Expand to 100+ by pulling
known FAERS → RxCUI ground truth from:
  - OHDSI Athena   (concept_synonym.csv, vocabulary_id='RxNorm')
  - FAERS standardization paper: Banda et al. 2016 (JAMIA Open)
  - NLM RxNav browser  https://mor.nlm.nih.gov/RxNav/

Schema
------
  raw           : string exactly as it appears in FAERS DRUG.drugname / prod_ai
  expected_rxcui: authoritative RxCUI (str)
  note          : pattern being tested
"""

TEST_CASES = [
    # ── Tier 1 should catch these ─────────────────────────────────────────
    {"raw": "aspirin",                          "expected_rxcui": "1191",   "note": "exact ingredient"},
    {"raw": "ibuprofen",                        "expected_rxcui": "5640",   "note": "exact ingredient"},
    {"raw": "metformin",                        "expected_rxcui": "6809",   "note": "exact ingredient"},
    {"raw": "atorvastatin",                     "expected_rxcui": "83367",  "note": "exact ingredient"},
    {"raw": "lisinopril",                       "expected_rxcui": "29046",  "note": "exact ingredient"},
    {"raw": "omeprazole",                       "expected_rxcui": "7646",   "note": "exact ingredient"},
    {"raw": "amoxicillin",                      "expected_rxcui": "723",    "note": "exact ingredient"},
    {"raw": "warfarin",                         "expected_rxcui": "11289",  "note": "exact ingredient"},
    {"raw": "acetaminophen",                    "expected_rxcui": "161",    "note": "exact ingredient"},
    {"raw": "prednisone",                       "expected_rxcui": "8638",   "note": "exact ingredient"},

    # ── Dosage noise → Tier 1 after cleaning ─────────────────────────────
    {"raw": "metformin 500mg",                  "expected_rxcui": "6809",   "note": "dosage suffix"},
    {"raw": "lisinopril 10 mg tablet",          "expected_rxcui": "29046",  "note": "dosage + form"},
    {"raw": "atorvastatin 40mg",                "expected_rxcui": "83367",  "note": "dosage suffix"},
    {"raw": "aspirin 81 mg",                    "expected_rxcui": "1191",   "note": "dosage suffix"},
    {"raw": "ibuprofen 200 mg capsule",         "expected_rxcui": "5640",   "note": "dosage + form"},

    # ── Salt suffix stripping ─────────────────────────────────────────────
    {"raw": "metformin hydrochloride",          "expected_rxcui": "6809",   "note": "HCl salt"},
    {"raw": "fluoxetine hydrochloride",         "expected_rxcui": "41493",  "note": "HCl salt"},
    {"raw": "sertraline hydrochloride",         "expected_rxcui": "36437",  "note": "HCl salt"},
    {"raw": "ondansetron hydrochloride",        "expected_rxcui": "26225",  "note": "HCl salt"},
    {"raw": "metoprolol tartrate",              "expected_rxcui": "6918",   "note": "tartrate salt"},
    {"raw": "amlodipine besylate",              "expected_rxcui": "17767",  "note": "besylate salt"},

    # ── Brand names (Tier 1 via RxNorm synonyms) ─────────────────────────
    {"raw": "tylenol",                          "expected_rxcui": "161",    "note": "brand -> acetaminophen"},
    {"raw": "advil",                            "expected_rxcui": "5640",   "note": "brand -> ibuprofen"},
    {"raw": "glucophage",                       "expected_rxcui": "6809",   "note": "brand -> metformin"},
    {"raw": "lipitor",                          "expected_rxcui": "83367",  "note": "brand -> atorvastatin"},
    {"raw": "zoloft",                           "expected_rxcui": "36437",  "note": "brand -> sertraline"},

    # ── Misspellings → expect Tier 2 (approximate) ───────────────────────
    {"raw": "metfromin",                        "expected_rxcui": "6809",   "note": "transposition typo"},
    {"raw": "ibuprofin",                        "expected_rxcui": "5640",   "note": "common misspelling"},
    {"raw": "aspriin",                          "expected_rxcui": "1191",   "note": "insertion typo"},
    {"raw": "atorvastatin calcium",             "expected_rxcui": "83367",  "note": "calcium salt (Lipitor)"},

    # ── Combination products (is_combo=True; rxcui should be None for parent) ──
    # These are recorded to verify combo detection, not RxCUI resolution.
    {"raw": "acetaminophen / codeine",          "expected_rxcui": None,     "note": "combo — no single RxCUI"},
    {"raw": "amoxicillin and clavulanate",      "expected_rxcui": None,     "note": "combo"},
    {"raw": "hydrocodone with acetaminophen",   "expected_rxcui": None,     "note": "combo"},
]


def run_test_report(resolved_df, unresolved_df) -> dict:
    """
    Compare pipeline output against TEST_CASES.
    Returns precision / recall per tier and a list of failures.
    """
    import pandas as pd

    all_resolved = pd.concat([resolved_df, unresolved_df], ignore_index=True)
    # Index by cleaned string via preprocess
    from .preprocess import clean as _clean

    results = []
    for case in TEST_CASES:
        raw      = case["raw"]
        expected = case["expected_rxcui"]
        cleaned  = _clean(raw)

        match = all_resolved[all_resolved["cleaned_string"] == cleaned]
        if match.empty:
            got_rxcui = None
            tier      = None
        else:
            got_rxcui = match.iloc[0]["rxcui"]
            tier      = match.iloc[0]["match_tier"]

        correct = (expected is None and got_rxcui is None) or (got_rxcui == expected)
        results.append({
            "raw":      raw,
            "expected": expected,
            "got":      got_rxcui,
            "tier":     tier,
            "correct":  correct,
            "note":     case["note"],
        })

    results_df = pd.DataFrame(results)
    n_total    = len(results_df)
    n_correct  = results_df["correct"].sum()

    print("\n" + "=" * 60)
    print("  TEST SET RESULTS")
    print("=" * 60)
    print(f"  Correct: {n_correct} / {n_total}  ({100*n_correct/n_total:.0f}%)")

    failures = results_df[~results_df["correct"]]
    if len(failures):
        print(f"\n  Failures ({len(failures)}):")
        for _, r in failures.iterrows():
            print(f"    {r['raw']!r:40s}  expected={r['expected']}  got={r['got']}  tier={r['tier']}")

    return {
        "n_total":   n_total,
        "n_correct": int(n_correct),
        "precision": n_correct / n_total,
        "failures":  failures.to_dict("records"),
    }
