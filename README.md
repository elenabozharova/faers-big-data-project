# FAERS temporal reporting-pattern project

This course project processes FDA Adverse Event Reporting System (FAERS)
reports to study whether an annual reporting pattern for a tuple
`(drug_a, drug_b, reaction_pt)` recurs in a future year.

The intended data range is 2012–2025. Drug identifiers are ChEMBL molecule
identifiers such as `CHEMBL25`.

## Active pipeline

### Map FAERS drug names to ChEMBL identifiers

`src/map_faers_to_chembl.py` normalizes names with
`src/normalize_drug_names.py`, queries ChEMBL names and synonyms, and writes
`data/processed/faers/faers_drug_to_chembl.parquet`.

It requires FAERS DRUG data and either a local ChEMBL PostgreSQL database or
an existing `chembl_name_index.parquet` cache.

```powershell
$env:PGPASSWORD = "your_password"
python src/map_faers_to_chembl.py
```

### Run case-level preprocessing

`src/process_faers_case_level.py` processes one quarterly archive or an
extracted FAERS directory. It resolves `CASEID`/`CASEVERSION`, deduplicates
case reactions and mapped ChEMBL drugs, creates canonical unordered drug
pairs, and counts distinct cases for each pair–reaction tuple.

```powershell
python src/process_faers_case_level.py --quarter 2025Q4
```

### Run annual preprocessing

`src/process_faers_case_level_annual.py` concatenates all four quarters for
each requested year before resolving case versions and producing annual
distinct-case counts.

```powershell
python src/process_faers_case_level_annual.py --years 2023 2024 2025
```

Use `--help` on each command for input paths, mapping options, role filtering,
and smoke-test row limits.

## Configuration

ChEMBL connection and single-quarter FAERS defaults are read from environment
variables documented in `.env.example`.

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Raw and processed data are excluded from Git.

## Milestone 2: reproducible FAERS acquisition

Official quarterly ASCII archive URLs and required members are declared in
the tracked `config/faers_sources.csv` manifest. Local attempts and checksums
are appended to the untracked provenance log at
`data/interim/acquisition/faers_download_log.csv`.

```powershell
# Pilot download
python src/download_faers.py --years 2024 --quarters 1

# One full year
python src/download_faers.py --years 2024

# Three-year pilot
python src/download_faers.py --years 2023 2024 2025

# Verify existing files without network requests
python src/download_faers.py --years 2023 2024 2025 --verify-only

# Validate and show intended actions without writing anything
python src/download_faers.py --years 2024 --dry-run
```

Downloads stream to a same-directory `.part` file, validate as ZIPs (including
CRC and required DEMO/DRUG/REAC members), and are atomically installed only
after all checks pass. Repeated normal runs validate and skip an unchanged
archive. A checksum differing from the latest successful local log entry is a
failure; `--overwrite` deliberately downloads a replacement while preserving
the current archive unless the replacement passes validation. Raw archives,
partial downloads, runtime logs, and generated data are not committed.

The 2012 data remains a later compatibility stage: Q1-Q3 use legacy AERS
files, and Q3 is a partial quarter.

## Milestone 3: raw schema and quality profiles

The read-only profiler scans FAERS tables inside each ZIP in bounded pandas
chunks and writes `raw_schema_profile.csv` plus `raw_schema_summary.json`.

```powershell
python src/profile_faers_raw.py --years 2024 --quarters 1
python src/profile_faers_raw.py --years 2024
python src/profile_faers_raw.py --years 2023 2024 2025
python src/profile_faers_raw.py --years 2024 --quarters 1 --max-rows-per-table 5000 --out-dir reports/data_quality/sample
```

Profiles preserve ordered raw column names and use trimmed, lowercase names
only for logical comparisons and deterministic schema signatures. Zero means
an applicable check found no issue; a blank/null metric is not applicable or
depends on an absent column. Repeated `PRIMARYID` values are expected in DRUG
and REAC; their duplicate keys are `PRIMARYID + DRUG_SEQ` and trimmed
`PRIMARYID + PT`. DEMO checks repeated `PRIMARYID` and `CASEID + CASEVERSION`.
`unique_raw_drug_names` is the count of nonblank, whitespace-trimmed DRUGNAME
values without case folding or drug normalization. Sample runs are marked and
cannot overwrite default full reports. Raw ZIPs are never extracted or changed.

## Current status

Implemented:

- FAERS table loading and diagnostics;
- within-run `CASEID`/`CASEVERSION` resolution;
- reaction normalization and case-level deduplication;
- FAERS drug-name normalization and ChEMBL mapping;
- canonical ChEMBL drug-pair generation;
- optional role-aware pair filtering;
- chunked pair–reaction joining;
- distinct-case support counts and annual output wrappers.

Not yet implemented:

- a consolidated multi-year Gold table;
- temporal labels and lag features;
- leakage-free train, validation, and test splits;
- persistence baselines and logistic regression;
- a PySpark implementation and scalability benchmark;
- final model evaluation and temporal visualizations.

The current annual wrapper resolves report versions within each year. A
cross-year policy for cases updated in later years must be defined before the
temporal modeling dataset is finalized.
