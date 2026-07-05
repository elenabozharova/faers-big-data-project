"""Inspect the pipeline output parquet — unresolved strings and RxCUI levels."""
import sys
import os
_SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import pandas as pd

path = os.path.join(os.path.dirname(_SRC), "data", "processed", "drug_mapping_rxnorm.parquet")
df = pd.read_parquet(path)

unresolved = df[df["rxcui"].isna()].sort_values("cleaned_string")
print(f"\nUnresolved ({len(unresolved)}) — sample:")
print(unresolved["cleaned_string"].head(40).to_string())

# Check RxCUI levels for known cases
print("\n\nRxCUI lookup for known exact names:")
known = ["omeprazole", "amoxicillin", "prednisone", "ibuprofen", "aspirin",
         "atorvastatin", "metformin hydrochloride", "amlodipine besylate"]
for name in known:
    row = df[df["cleaned_string"] == name]
    if len(row):
        rxcui = row.iloc[0]["rxcui"]
        tier  = row.iloc[0]["match_tier"]
        print(f"  {name:<35s} rxcui={rxcui}  tier={tier}")
    else:
        print(f"  {name:<35s} NOT IN 10k SAMPLE")
