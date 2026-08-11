# FAERS temporal reporting-pattern project

This course project processes FDA Adverse Event Reporting System (FAERS)
reports to study whether an annual reporting pattern for a tuple
`(drug_a, drug_b, reaction_pt)` recurs in a future year. The intended data
range is 2012-2025, and standardized drug identifiers are ChEMBL molecule IDs
such as `CHEMBL25`.

## Milestone status

| Milestone | Status |
|---|---|
| MS1 | Complete |
| MS2 | Complete: reproducible raw FAERS acquisition and provenance |
| MS3 | Complete: raw schema and quality profiling |
| MS4A | Complete and production-validated: global cross-year CASEID resolution |
| MS4B | Complete and production-validated end to end |
| MS5 | Complete and production-validated for 2023-2025: reusable multi-year drug vocabulary and conservative versioned exact mapping to ChEMBL 36 |
| MS6 | Complete and production-validated for 2023-2025: final annual drug-pair-reaction support datasets |
| MS7 | Complete and production-validated for 2023-2025: canonical partitioned annual pair-reaction Gold dataset |
| MS8 | Complete and production-validated for the 2023-2025 Stage 3 pilot: leakage-safe temporal features and recurrence targets |
| MS9 and later | Pending; not started |

## Pipeline overview

The active production path is:

```text
FAERS quarterly archives
-> raw schema and quality profiles
-> global cross-year CASEID assignment and winning-version selection
-> annual MS4B resolved DRUG and REAC rows
-> multi-year normalized drug vocabulary
-> verified ChEMBL 36 preferred-name and synonym index
-> versioned exact FAERS-to-ChEMBL mapping
-> annual mapping coverage reports
-> final globally corrected annual pair-reaction support datasets
-> canonical partitioned MS7 Gold dataset
-> temporal feature engineering and modeling (pending)
```

MS4A assigns each CASEID once and selects one eligible version globally.
MS4B routes that selection back to its physical source partition and produces
the corrected annual DRUG/REAC boundary artifacts. The complete resolution
policy is documented in
[Cross-year FAERS case resolution](docs/cross_year_case_resolution.md).

## Environment and configuration

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

FAERS and ChEMBL runtime settings are read from environment variables described
in `.env.example`. Database credentials must be supplied through the local
execution environment; they are not stored in source, configuration,
documentation, or generated metadata. Raw and processed data are excluded
from Git under the repository data policy.

## MS2 - Reproducible FAERS acquisition

Official quarterly ASCII archive URLs and required members are declared in
`config/faers_sources.csv`. Local attempts and checksums are appended to the
untracked provenance log at
`data/interim/acquisition/faers_download_log.csv`.

```powershell
python src/download_faers.py --years 2023 2024 2025
python src/download_faers.py --years 2023 2024 2025 --verify-only
```

Downloads stream to same-directory `.part` files, validate ZIP structure, CRC,
and required DEMO/DRUG/REAC members, and are installed atomically. The 2012
data remains a later compatibility stage because Q1-Q3 use legacy AERS files
and Q3 is a partial quarter.

## MS3 - Raw schema and quality profiles

The read-only profiler scans FAERS tables inside each ZIP in bounded chunks
and writes `reports/data_quality/raw_schema_profile.csv` and
`reports/data_quality/raw_schema_summary.json`.

```powershell
python src/profile_faers_raw.py --years 2023 2024 2025
```

Raw ZIPs are not extracted or changed by profiling.

## MS4 - Global case resolution boundary

```powershell
python src/build_global_case_index.py --years 2023 2024 2025
python src/process_faers_case_level_annual.py --years 2023 2024 2025
```

MS4A and MS4B are complete and production-validated end to end. The deferred
downstream closure reused the corrected annual `retained_cases.parquet`,
`resolved_drug_rows.parquet`, and `resolved_reac_rows.parquet` files without
rerunning raw integration, then proved that every downstream CASEID remained
inside its globally assigned annual population.

## MS5 - Multi-Year Drug Vocabulary and Exact ChEMBL Mapping

FAERS drug names are free text. The same drug can appear with different
capitalization, spacing, strengths, formulations, brand names, and spelling
variants, so raw names are not reliable identifiers for temporal drug-pair
analysis. MS5 creates one consistent multi-year vocabulary and maps recognized
names to standardized ChEMBL identifiers.

Normalization and mapping are separate operations:

- **Normalization** cleans and standardizes text into a deterministic key.
- **ChEMBL mapping** links that normalized key to a recognized ChEMBL molecule
  ID under an explicit matching policy.

Fuzzy matching is not part of normalization, and no fuzzy mapping was accepted
in the production artifact.

### Processing flow and provenance

```text
Globally resolved DRUG rows
-> preferred raw-name selection
-> frozen drug-name normalization
-> reusable multi-year normalized vocabulary
-> ChEMBL 36 preferred-name and synonym index
-> exact preferred-name mapping
-> exact synonym mapping
-> explicit ambiguous and combination outcomes
-> versioned mapping artifact
-> annual coverage reports
```

Only the globally corrected MS4B annual `resolved_drug_rows.parquet` files were
used. Raw unfiltered quarterly DRUG files were not used to build the
vocabulary, and later nonselected follow-up PRIMARYIDs did not re-enter MS5.
Analytical year comes from the validated annual MS4B output context, not the
physical source quarter.

### Production inputs

| Year | Resolved DRUG rows | Globally assigned CASEIDs |
|---:|---:|---:|
| 2023 | 4,763,980 | 1,361,395 |
| 2024 | 5,008,214 | 1,308,773 |
| 2025 | 4,942,013 | 1,295,857 |
| **Total** | **14,714,207** | **3,966,025** |

Of these rows, 14,714,102 were eligible for mapping. The remaining 105 rows
had normalized values that were blank or invalid vocabulary placeholders;
they were not missing raw DRUG rows.

### Vocabulary results

| Metric | Result |
|---|---:|
| Global normalized vocabulary names | 94,725 |
| Normalization collisions | 14,972 |
| Combination names | 8,419 |
| Blank or invalid normalized rows | 105 |

A normalization collision means that multiple selected raw names produced the
same normalized key. This can be desirable, for example when dosage-form
variants reduce to the same name, but collisions remain quality-review signals
because aggressive normalization can also merge unrelated values. Explicit
combination names contain multiple components separated by the normalized `|`
separator.

### ChEMBL 36 index

| Metric | Result |
|---|---:|
| ChEMBL release | 36 |
| Preferred-name and synonym relationships | 150,008 |
| Unique normalized ChEMBL names | 102,527 |
| Unique ChEMBL IDs | 50,614 |
| Ambiguous normalized ChEMBL names | 3,160 |

The local database name, installed version record, and required tables were
verified before artifact generation. Preferred names and synonyms use the same
frozen normalizer as FAERS names. Multiple ChEMBL IDs associated with one
normalized name are preserved so ambiguity can be detected rather than hidden
by selecting the first record.

### Exact mapping results

| Metric | Result |
|---|---:|
| FAERS vocabulary names | 94,725 |
| Single-ID mapped names | 8,624 |
| Single-ID mapped-name coverage | 9.1042% |
| Eligible DRUG rows | 14,714,102 |
| Single-ID mapped rows | 13,138,103 |
| Single-ID row coverage | 89.2892% |
| Exact preferred-name mappings | 5,278 |
| Exact synonym mappings | 3,346 |
| Ordinary ambiguous names | 325 |
| Fully recognized combinations | 1,345 |
| Partial combinations | 6,248 |
| Ambiguous combinations | 826 |
| Ordinary unmatched names | 77,357 |

The low unique-name coverage and high row-level coverage show a long tail of
rare, noisy, or product-specific FAERS names. Frequently occurring names map
much more successfully than the many rare vocabulary variants. Specifically,
89.2892% of eligible DRUG rows received one unambiguous exact ChEMBL mapping.
This is coverage, not mapping accuracy; accuracy would require manual
validation against a reviewed sample.

### Annual coverage

| Year | Eligible rows | Single-ID row coverage | Unique names | Single-ID name coverage |
|---:|---:|---:|---:|---:|
| 2023 | 4,763,947 | 87.7468% | 45,293 | 14.2053% |
| 2024 | 5,008,181 | 89.0215% | 41,365 | 15.0949% |
| 2025 | 4,941,974 | 91.0473% | 44,875 | 13.9699% |
| **ALL** | **14,714,102** | **89.2892%** | **94,725** | **9.1042%** |

The `ALL` name denominator is the union of names, not the sum of annual unique
counts or an average of annual percentages. Common mapped names recur in
several years but count once in the union, while rare unmatched year-specific
variants accumulate. This is why combined name coverage is lower than each
annual percentage.

### Mapping policy

Ordinary names use this accepted order:

1. Exact preferred ChEMBL name.
2. Exact ChEMBL synonym.
3. Unmatched.

One unique exact candidate produces `mapping_status = mapped`. Multiple IDs
produce an explicit ambiguous status and the first candidate is never selected
arbitrarily. No fuzzy match contributes to mapped coverage.

For explicit combinations, each ordered component is mapped independently.
Fully recognized combinations retain their component ChEMBL IDs but do not
expose one single `chembl_id`. Partial and ambiguous combinations remain
explicit non-single-ID outcomes; they are not counted as ordinary mapped
drugs.

### Reproduction commands

The vocabulary builder accepts arbitrary requested analytical years:

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

$env:PGPASSWORD = "<local ChEMBL database password>"
python src/build_chembl_name_index.py

python src/build_faers_chembl_mapping.py `
  --years 2023-2025 `
  --vocabulary data/interim/drug_mapping/multi_year_drug_vocabulary_2023_2025_norm_v1.parquet `
  --name-index data/interim/drug_mapping/chembl36_name_index_norm_v1_index_v1.parquet `
  --output data/processed/mappings/faers_to_chembl_2023_2025_chembl36_norm_v1_exact_v1.parquet
```

An expanded vocabulary can be reproduced after globally corrected annual
inputs exist for every year in the requested range:

```powershell
python src/build_multi_year_drug_vocabulary.py `
  --years 2013-2025 `
  --annual-input-root data/processed/faers_case_level_annual `
  --input-filename resolved_drug_rows.parquet `
  --output data/interim/drug_mapping/multi_year_drug_vocabulary_2013_2025_norm_v1.parquet `
  --report-dir reports/data_quality `
  --normalization-version norm_v1 `
  --artifact-version vocab_v2 `
  --sample-raw-names-limit 10
```

The 2013-2025 command must not be used until globally corrected annual inputs
exist for 2013-2022. Existing outputs require the documented `--overwrite`
option when intentional replacement is needed.

### Principal artifacts and quality reports

Versioned data artifacts:

- `data/interim/drug_mapping/multi_year_drug_vocabulary_2023_2025_norm_v1.parquet`
- `data/interim/drug_mapping/chembl36_name_index_norm_v1_index_v1.parquet`
- `data/processed/mappings/faers_to_chembl_2023_2025_chembl36_norm_v1_exact_v1.parquet`

Quality reports:

- `reports/data_quality/ms5_drug_vocabulary_by_year.csv`
- `reports/data_quality/ms5_drug_vocabulary_summary.json`
- `reports/data_quality/ms5_normalization_collisions.csv`
- `reports/data_quality/ms5_combination_names.csv`
- `reports/data_quality/ms5_chembl36_name_index_summary.json`
- `reports/data_quality/ms5_chembl_mapping_summary.json`
- `reports/data_quality/ms5_chembl_mapping_coverage_by_year.csv`
- `reports/data_quality/ms5_chembl_mapping_status_by_year.csv`
- `reports/data_quality/ms5_chembl_unmatched_names_top.csv`
- `reports/data_quality/ms5_chembl_ambiguous_mappings.csv`
- `reports/data_quality/ms5_chembl_combination_mapping.csv`

Generated Parquet files and other large data artifacts are normally excluded
from Git and must be reproduced locally.

### Known limitations

- Fuzzy mappings are not accepted, RxNorm is not integrated, and no manual
  mapping campaign was performed.
- Salt forms are not unified with parent ingredients.
- Unicode normalization remains unchanged.
- The frozen normalizer removes hyphens according to its implemented regex,
  and its existing `without` separator behavior remains unchanged.
- Combination products are not converted into one ChEMBL product entity.
- Space-separated implicit combinations may remain unmatched.
- Exact mapping coverage is not manually validated mapping accuracy.

Examples such as `acetaminophen oxycodone hydrochloride`, `carbidopa
levodopa`, `sulfamethoxazole trimethoprim`, and `acetaminophen hydrocodone
bitartrate` likely contain multiple active ingredients but have no explicit
separator. The exact-only mapper deliberately does not split on whitespace
because valid single drugs and salt forms can also contain several words.
These examples are not assertions that the source records are incorrect.

### Detailed MS5 documentation

- [Current mapping audit](docs/ms5_current_mapping_audit.md): the pre-MS5
  mapping, normalization, consumer, and provenance audit.
- [Reusable multi-year vocabulary](docs/ms5_multi_year_drug_vocabulary.md):
  vocabulary design, deterministic processing, reports, and production
  results.
- [Versioned ChEMBL exact mapping](docs/ms5_chembl_exact_mapping.md): ChEMBL
  index verification, exact mapping policy, schemas, consumer contract,
  coverage, and limitations.

### Validation

Final MS5 validation recorded:

- 72 unit tests passed;
- 12 normalization smoke tests passed;
- Python compilation passed;
- CLI help checks passed;
- `pip check` reported no broken requirements;
- `git diff --check` passed with existing line-ending warnings only;
- resolved input hashes and raw ZIP hashes remained unchanged;
- no temporary SQLite databases or `.part` files remained; and
- no absolute Windows paths were persisted in generated reports.

`pytest` is not required by the project test workflow. When it is unavailable,
that is an environment condition rather than a code failure; the authoritative
suite uses `unittest`.

## Milestone development boundary

Each milestone is committed separately. MS5 source code, tests, and
documentation should be committed together only after production validation.
Generated large Parquet artifacts remain excluded according to the repository
data policy.

## MS6 - Final annual drug-pair-reaction datasets

MS6 is complete and production-validated for 2023-2025. It reused the immutable
MS4B boundary, applied the versioned MS5 exact mapping with fuzzy mappings
disabled, retained only ordinary single-ID mappings, generated deterministic
unordered pairs under the existing `all_pairs` policy, joined pairs to
case-level reactions in bounded chunks, and counted distinct contributing
CASEIDs.

| Year | Cases with >=2 mapped drugs | Pair rows | Pair-reaction rows | Unique patterns | Min3 | Min5 | Min10 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2023 | 449,412 | 11,134,568 | 133,383,154 | 23,256,862 | 6,393,185 | 4,234,914 | 2,390,729 |
| 2024 | 425,838 | 13,313,358 | 259,927,804 | 25,346,924 | 7,067,742 | 4,789,551 | 3,012,592 |
| 2025 | 427,407 | 14,122,469 | 209,952,427 | 28,871,255 | 8,434,213 | 5,658,578 | 3,326,862 |

Authoritative outputs are under
`data/processed/faers_case_level_annual_final/<year>/`. The production run used
`data/processed/mappings/faers_to_chembl_2023_2025_chembl36_norm_v1_exact_v1.parquet`
with SHA-256
`e4e28c86dea312c262ce9ab92c006430817d6c1dbfb0a348a33e1dab39d408a5`.
All raw ZIP, global-assignment, mapping, and reused boundary hashes remained
unchanged. See
[MS6 annual pair-reaction rebuild](docs/ms6_annual_pair_reaction_rebuild.md)
for usability counts, performance, validations, legacy availability, and
reproduction details. Final repository validation passed 77 unit tests, Python
compilation, both production CLI help checks, dependency checks, and Git diff
hygiene.

## MS7 - Canonical annual pair-reaction Gold dataset

MS7 is complete and production-validated for 2023-2025. It publishes exactly
`year`, `drug_a`, `drug_b`, `reaction_pt`, and `n_cases` under
`data/gold/annual_pair_reactions/year=<year>/`. The three partitions contain
77,475,041 rows in total. All pair ordering, support, key uniqueness, ChEMBL
identifier, partition-year, nullability, and global CASEID assignment
invariants passed. See [MS7 Gold dataset](docs/ms7_gold_annual_pair_reactions.md).

## MS8 - Temporal recurrence feature pilot

MS8 is complete and production-validated for the Stage 3 pilot. It generated
62,245,721 candidate rows under
`data/gold/temporal_pair_reaction_features/`, partitioned by target year. The
candidate population contains tuples observed at least once before the target
year; tuples first appearing in the target year are excluded as cold starts.
All features use history strictly before the target year, and unavailable pilot
lags are null rather than false zero observations. See
`reports/data_quality/ms8_temporal_features_validation.json`.

## Next milestone

MS9 and later have not started. The required post-MS7 sequence is documented
in [Post-MS7 temporal modeling and Spark requirements](docs/post_ms7_temporal_modeling_requirements.md):

1. Implement an always-negative baseline, previous-year persistence rule,
   support-threshold rule, logistic regression, and optionally
   gradient-boosted trees, in that order.
2. Complete the Stage 4 historical expansion for 2013-2022 after the pilot
   pipeline works end to end.
3. Port the expensive data-engineering stages to Spark and benchmark the local
   pipeline against Spark over increasing year ranges.

The expanded 2013-2022 history remains required before declaring the full
three-lag modeling dataset production-complete for the intended historical
study.
