# Local data layout

FAERS ZIP files under `data/raw/faers/` are downloaded locally from the
official FDA URLs recorded in `config/faers_sources.csv`. Raw archives are
immutable inputs: do not edit, extract over, or rewrite them. Raw, interim,
processed, and Gold data are generated artifacts and are not committed.

Download and validate a year:

```powershell
python src/download_faers.py --years 2024
```

Verify already-downloaded archives without network access:

```powershell
python src/download_faers.py --years 2024 --verify-only
```

The tracked manifest is the reproducible declaration of official source URLs,
filenames, and required ZIP members. The untracked runtime log at
`data/interim/acquisition/faers_download_log.csv` is append-only local
provenance: it records attempts, checksums, HTTP metadata, and failures.

Raw archives can be profiled in place without extraction or modification:

```powershell
python src/profile_faers_raw.py --years 2024
```

Milestone 4A writes ignored Parquet indexes under
`data/interim/case_resolution/`. Temporary SQLite work files are cleaned by
default; `--keep-work-files` retains them for debugging.

Milestone 7 publishes the canonical five-column analytical dataset under
`data/gold/annual_pair_reactions/year=<year>/`. Gold partitions are generated
artifacts and should be rebuilt from validated MS6 annual aggregates rather
than edited directly.

Milestone 8 publishes the Stage 3 temporal feature pilot under
`data/gold/temporal_pair_reaction_features/target_year=<year>/`. These ignored
Parquet partitions are rebuilt from the canonical MS7 Gold inputs.
