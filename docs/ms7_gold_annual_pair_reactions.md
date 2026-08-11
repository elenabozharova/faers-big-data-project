# MS7A canonical annual pair-reaction Gold dataset

MS7A publishes the stable analytical boundary consumed by future temporal
feature engineering and modeling. It does not create lag features, labels,
train/test splits, or models.

## Contract

The dataset has exactly five columns:

| Column | Type | Meaning |
|---|---|---|
| `year` | int32 | Globally assigned analytical year |
| `drug_a` | string | Lower-sorting ChEMBL molecule identifier |
| `drug_b` | string | Higher-sorting ChEMBL molecule identifier |
| `reaction_pt` | string | FAERS reaction preferred term |
| `n_cases` | int64 | Distinct contributing CASEIDs |

It is horizontally partitioned using Hive-style year directories:

```text
data/gold/annual_pair_reactions/
├── year=2023/
├── year=2024/
└── year=2025/
```

The builder supports arbitrary years and inclusive ranges. The production
repository currently has validated MS6 inputs only for 2023-2025; extending
the Gold dataset to 2013-2022 first requires those years to pass MS2-MS6.

## Invariants

- `drug_a < drug_b` and therefore `drug_a != drug_b`;
- `n_cases >= 1`;
- exactly one row per `(year, drug_a, drug_b, reaction_pt)`;
- every drug identifier matches `^CHEMBL[0-9]+$`;
- every included CASEID has at most one globally assigned analytical year;
- every row carries the year of its physical partition; and
- no analytical column is null.

Input rows must also be deterministically ordered by the complete Gold key.
This makes cross-batch duplicate detection bounded-memory and ensures stable
output ordering.

## Build

```powershell
python src/build_gold_annual_pair_reactions.py `
  --years 2023-2025
```

Existing output is protected unless `--overwrite` is supplied. Publication is
staged beside the final directory and renamed only after all requested years
pass validation. A successful dataset contains `_SUCCESS`, `_metadata.json`,
and one Parquet part per year. The validation report is written to
`reports/data_quality/ms7_gold_validation.json`.

## Production result

The 2023-2025 production build completed successfully with 77,475,041 rows:

| Year | Rows |
|---:|---:|
| 2023 | 23,256,862 |
| 2024 | 25,346,924 |
| 2025 | 28,871,255 |

All declared invariants passed. The Gold row counts exactly match the
corresponding validated MS6 unique pair-reaction pattern counts.
