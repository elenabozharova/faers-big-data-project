# MS5 multi-year normalized drug vocabulary

## Purpose and boundary

`src/build_multi_year_drug_vocabulary.py` builds the reusable normalized-name
population for later ChEMBL mapping. It reads only globally corrected MS4B
`resolved_drug_rows.parquet` files. It does not read raw or legacy DRUG rows,
perform ChEMBL/RxNorm matching, generate manual mappings, or create drug pairs.

## Input contract and analytical years

Each requested year resolves exactly one file at
`<annual-input-root>/<year>/<input-filename>`. Required projected columns are
`caseid`, `primaryid`, `prod_ai`, and `drugname`. Analytical year comes only
from the explicit annual input context; `source_quarter` and clinical dates do
not assign it.

Every resolved `(caseid, primaryid)` is validated against MS4A rows satisfying
`assignment_status == "included"` for that analytical year. Follow-up audit
flags are not used. Identifiers are trimmed strings and are never numerically
coerced.

Year syntax supports ranges and explicit nonconsecutive values:

```powershell
python src/build_multi_year_drug_vocabulary.py --years 2023-2025
python src/build_multi_year_drug_vocabulary.py --years 2018 2020 2023-2025
```

Overlapping/duplicate years, descending ranges, missing inputs, and invalid
tokens fail rather than being silently skipped.

## Frozen selected-name and normalization rules

The builder imports `_choose_raw_drug_name()` from the production case-level
pipeline. A trimmed `prod_ai` longer than two characters and not numeric after
literal periods are removed wins; otherwise trimmed `drugname` is used.

It then imports `normalize_drug_string()` from `normalize_drug_names.py`.
MS5 vocabulary construction intentionally does not alter its dosage,
punctuation, hyphen, Unicode, parenthetical, salt, `without`, or combination
behavior. The open policy questions are measured through collision,
combination, and review reports before any future normalization change.

## Vocabulary and blank policy

The Parquet artifact has one lexically sorted row per nonblank
`normalized_drug_name`:

- `first_seen_year`, `last_seen_year`, sorted list `years_seen`, `year_count`;
- `raw_row_count`, `distinct_case_count`, `distinct_primaryid_count`;
- `distinct_raw_name_count`;
- deterministic, lexically sorted, bounded list `sample_raw_names`.

Actual null/blank raw fields, blank chosen names, and names that normalize to
blank are counted in reports. Blank normalized values are not vocabulary keys.
Normalized values equal to the common serialized-null placeholders `nan`,
`none`, or `null` (case-insensitive) are also treated as blank eligibility
outcomes, without changing the frozen selector or normalizer.

## Determinism and scalability

Inputs are sorted by analytical year and read with projected PyArrow batches.
The builder never concatenates all annual rows. Batch counts and exact distinct
`(normalized name, selected PRIMARYID)` and `(normalized name, raw name)` keys
spill into a controlled temporary SQLite database. Because MS4A validates a
one-to-one selected PRIMARYID/CASEID contract and assigns each CASEID once,
distinct PRIMARYID and CASEID counts are equal for this population. SQLite is
closed and its temporary directory is removed on success or failure.

Vocabulary rows, year lists, raw-name samples, collision ordering, combination
ordering, and the stable-hash review sample are deterministic and independent
of requested input order, source row order, filesystem order, and Python hash
randomization.

## Reports and metadata

- `ms5_drug_vocabulary_by_year.csv`: annual input, blank, and vocabulary counts;
- `ms5_drug_vocabulary_summary.json`: scope, hashes, versions, commit/dirty
  state, runtime, paths, and invariant results;
- `ms5_normalization_collisions.csv`: normalized keys produced by more than one
  chosen raw name;
- `ms5_combination_names.csv`: names containing `|`, without mapping/splitting;
- `ms5_normalization_review_sample.csv`: deterministic union of frequent,
  collision, combination, stable-hash, and large raw/normalized-difference
  samples, with blank human-review fields.

Artifacts record `artifact_version`, `normalization_version`, requested years,
input paths and SHA-256 hashes, creation time, Git commit, and dirty-worktree
state. No ChEMBL version is recorded because this stage performs no mapping.
Existing outputs are protected unless `--overwrite` is explicit.

## First production run: 2023–2025

Command:

```powershell
python src/build_multi_year_drug_vocabulary.py `
  --years 2023-2025 `
  --annual-input-root data/processed/faers_case_level_annual `
  --input-filename resolved_drug_rows.parquet `
  --output data/interim/drug_mapping/multi_year_drug_vocabulary_2023_2025_norm_v1.parquet `
  --report-dir reports/data_quality `
  --normalization-version norm_v1 `
  --artifact-version vocab_v1 `
  --sample-raw-names-limit 10
```

| Year | Rows | CASEIDs / selected PRIMARYIDs | Unique normalized names | Blank normalized rows | Combination names |
|---:|---:|---:|---:|---:|---:|
| 2023 | 4,763,980 | 1,361,395 | 45,293 | 33 | 3,384 |
| 2024 | 5,008,214 | 1,308,773 | 41,365 | 33 | 2,970 |
| 2025 | 4,942,013 | 1,295,857 | 44,875 | 39 | 3,862 |
| Total/scope | 14,714,207 | 3,966,025 | 94,725 | 105 | 8,419 unique |

The scope contains 14,972 normalized-name collision groups. No row had a
blank chosen name; 104 excluded rows had nonblank selected strings that the
frozen normalizer reduced to blank, and one additional `None!` value normalized
to the prohibited placeholder `none`. Total builder runtime was 345.54
seconds. Observed process working memory was approximately 1.13 GB while exact
distinct keys spilled to temporary SQLite; no staging file remained.

Prominent high-frequency collisions include `tirzepatide` (3 raw names,
427,100 rows), `prednisone` (60; 184,366), `methotrexate` (16; 163,302),
`acetaminophen` (82; 161,346), and `aspirin` (89; 95,117). These are audit
signals, not evidence that mappings should be merged or changed.

## Future 2013–2025 run

After globally corrected annual inputs exist for every requested year:

```powershell
python src/build_multi_year_drug_vocabulary.py --years 2013-2025 `
  --output data/interim/drug_mapping/multi_year_drug_vocabulary_2013_2025_norm_v1.parquet
```

The artifact must be rebuilt over the expanded scope so first-seen years and
all counts are recomputed; current outputs must not simply be appended.

## Validation and limitations

Standard validation is:

```powershell
python -m compileall src tests
python -m unittest discover -s tests -p "test_*.py"
python src/build_multi_year_drug_vocabulary.py --help
python -m pip check
git diff --check
git status --short
```

This artifact freezes existing normalization behavior but does not approve it.
Hyphen removal, `without` as a separator, Unicode handling, dosage-related
collisions, salt fragmentation, brand-name fragmentation, and combination
representation remain explicit review questions before ChEMBL mapping.
