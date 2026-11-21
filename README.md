# FAERS temporal reporting-pattern project

This course project processes FDA Adverse Event Reporting System (FAERS)
reports to study whether an annual reporting pattern for a tuple
`(drug_a, drug_b, reaction_pt)` recurs in a future year. The study data
range is 2013-2025 (2013 is the first complete calendar year in the FAERS
format; 2012 Q1-Q3 are legacy AERS files and Q3 is partial), and standardized
drug identifiers are ChEMBL molecule IDs such as `CHEMBL25`.

## Pipeline overview

```text
FAERS quarterly archives
-> raw schema and quality profiles
-> global cross-year CASEID assignment and winning-version selection
-> annual resolved DRUG and REAC rows
-> multi-year normalized drug vocabulary
-> ChEMBL 36 preferred-name and synonym index
-> versioned exact FAERS-to-ChEMBL mapping
-> annual drug-pair-reaction support datasets
-> canonical partitioned Gold dataset
-> temporal recurrence features
-> baselines and logistic regression (Historical Baseline v1)
-> PySpark implementation, parity check and scaling benchmark
```

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
from Git under the repository data policy and must be reproduced locally.

## Stages

Exact production commands and arguments are in the linked documents.

### Acquisition

`src/download_faers.py` downloads the official quarterly ASCII archives
declared in `config/faers_sources.csv`, validates ZIP structure, CRC and the
required DEMO/DRUG/REAC members, and installs them atomically. Checksums are
logged to `data/interim/acquisition/faers_download_log.csv`.

### Raw profiling

`src/profile_faers_raw.py` scans the tables inside each ZIP read-only and
writes schema and quality profiles to `reports/data_quality/`.

### Case resolution

`src/build_global_case_index.py` assigns each CASEID to one report year and
selects one eligible version globally; `src/process_faers_case_level_annual.py`
produces the corrected annual DRUG/REAC rows. See
[Cross-year FAERS case resolution](docs/cross_year_case_resolution.md) and
[Historical expansion](docs/stage4_historical_expansion.md).

### Drug normalization and ChEMBL mapping

`src/build_multi_year_drug_vocabulary.py`, `src/build_chembl_name_index.py`
and `src/build_faers_chembl_mapping.py` normalize drug names and map them to
ChEMBL IDs by exact preferred name, then exact synonym. No fuzzy mapping is
accepted. See
[Multi-year vocabulary](docs/ms5_multi_year_drug_vocabulary.md),
[ChEMBL exact mapping](docs/ms5_chembl_exact_mapping.md) and the
[mapping audit](docs/ms5_current_mapping_audit.md) (known limitations).

### Pair-reaction datasets and Gold

`src/build_gold_annual_pair_reactions.py` publishes
`(year, drug_a, drug_b, reaction_pt, n_cases)` under `data/gold/`. See
[Annual pair-reaction rebuild](docs/ms6_annual_pair_reaction_rebuild.md) and
[Gold dataset](docs/ms7_gold_annual_pair_reactions.md).

### Temporal features and models

`src/build_temporal_features.py` builds features from history strictly before
each target year; `src/run_ms7c_models.py` runs the baselines and logistic
regression. Results are frozen as
[Historical Baseline v1](docs/historical_baseline_v1.md) in
`reports/modeling/historical_baseline_v1/`. See also the
[modeling specification](docs/ms7c_modeling_confluence.md).

### PySpark

`src/spark_pair_reactions.py` rebuilds the Gold contract in PySpark;
`src/run_spark_benchmark.py` and `src/analyze_spark_benchmark.py` run and
summarize the scaling benchmark. Parity and benchmark results are in
`reports/spark/`. See [PySpark implementation](docs/ms8a_spark_implementation.md)
and the [PySpark specification](docs/ms8_pyspark_confluence.md).

## Tests

```powershell
.venv\Scripts\python.exe -m unittest discover -s tests
```
