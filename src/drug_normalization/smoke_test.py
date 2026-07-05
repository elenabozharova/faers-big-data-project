"""Quick smoke test for Tier 0 preprocess logic."""
import sys
import os

_SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from drug_normalization.preprocess import preprocess, clean, strip_salts, is_combo, split_components

cases = [
    ("aspirin",                   "aspirin",        None,      False),
    ("metformin 500mg",           "metformin",      None,      False),
    ("lisinopril 10 mg tablet",   "lisinopril",     None,      False),
    ("metformin hydrochloride",   "metformin hydrochloride", "metformin", False),
    ("amlodipine besylate",       "amlodipine besylate",     "amlodipine", False),
    ("acetaminophen / codeine",   "acetaminophen / codeine",    None,  True),
    ("amoxicillin and clavulanate","amoxicillin and clavulanate", None, True),
    ("amoxicillin\\clavulanic acid","amoxicillin\\clavulanic acid", None, True),
    ("budesonide\\formoterol fumarate", "budesonide\\formoterol fumarate", "budesonide\\formoterol", True),
    ("ASPIRIN",                   "aspirin",        None,      False),
    ("Ibuprofen 200 MG Capsule",  "ibuprofen",      None,      False),
    ("atorvastatin calcium",      "atorvastatin calcium",    "atorvastatin", False),
]

print("Tier 0 Preprocess Smoke Test")
print("=" * 70)
passed = 0
failed = 0

for raw, exp_clean, exp_salt, exp_combo in cases:
    p = preprocess(raw)
    got_clean = p["cleaned_string"]
    got_salt  = p["salt_stripped"]
    got_combo = p["is_combo"]

    ok_clean = got_clean == exp_clean
    ok_salt  = got_salt  == exp_salt
    ok_combo = got_combo == exp_combo
    ok = ok_clean and ok_salt and ok_combo

    status = "PASS" if ok else "FAIL"
    if ok:
        passed += 1
    else:
        failed += 1

    print(f"  [{status}] {raw!r}")
    if not ok_clean:
        print(f"         clean:  got={got_clean!r}  exp={exp_clean!r}")
    if not ok_salt:
        print(f"         salt:   got={got_salt!r}  exp={exp_salt!r}")
    if not ok_combo:
        print(f"         combo:  got={got_combo!r}  exp={exp_combo!r}")

print("=" * 70)
print(f"  {passed} passed, {failed} failed")
sys.exit(0 if failed == 0 else 1)
