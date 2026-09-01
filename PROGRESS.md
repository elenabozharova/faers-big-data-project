# Project progress

Snapshot as of 2026-08-20. Base commit `3f5c4a0` ("Historical data
expansion") plus a large **uncommitted** session that completed
**Historical Baseline v1** (2013-2022 expansion + rolling 2016-2025 temporal
experiment). Nothing from this session has been committed yet.

## Milestone status

| Milestone | Status |
|---|---|
| MS1-MS3 | Complete |
| MS4A | Complete and production-validated for 2013-2025 (promoted from the Stage 4 sandbox this session) |
| MS4B | Complete and production-validated for 2013-2025 |
| MS5 | Complete and production-validated for 2023-2025 (untouched); additive artifacts for 2013-2022 |
| MS6 | Complete and production-validated for 2023-2025 (untouched); additive artifacts for 2013-2022 |
| MS7A | Complete and production-validated for 2023-2025 (untouched); additive Gold partitions for 2013-2022 (2013-2025 complete) |
| MS7B | Complete for the 2023-2025 pilot (untouched); full 2013-2025 rebuild done for target years 2016-2025 |
| MS7C | Complete for the Stage 3 pilot; **Historical Baseline v1 complete** — 30/30 model runs (2016-2025 x M1/M2/M3) |
| MS8A-MS8C | Pending (unchanged, out of scope for this session) |

## Historical Baseline v1 — complete

Full results, methodology, and interpretation: [docs/historical_baseline_v1.md](docs/historical_baseline_v1.md).
Plan file: `C:\Users\bozha\.claude\plans\melodic-kindling-hammock.md`.

**Headline finding**: `count_t1` is confirmed as the single strongest
predictor in every one of the 10 target years, but it does not dominate to
the exclusion of deeper history — adding `count_t2`/`count_t3` (M2 vs M1)
raises mean F1 from 0.442 to 0.520, mean AUROC from 0.729 to 0.785, and mean
AUPRC from 0.438 to 0.528, consistently across all 10 years. The pilot's
apparent "count_t1 completely dominates" reading was an artifact of its
training data having only one year of real history, making every deeper
feature a constant (zero-variance) column that the model was structurally
unable to learn from.

**2023-2025 production artifacts are untouched** — verified by SHA-256
(mapping file's checksum still matches the one published in the README).

**Three genuine pre-existing pipeline bugs found and fixed** (each with a
regression test), plus two shared-output-path bugs that would have destroyed
already-published years on any additive rerun:
1. Self-referential "later followup" false-positive assertion (`integrate_global_case_assignments.py`).
2. Non-UTF-8 (`cp1252`) bytes in the 2019 Q3 archive (`faers_utils.py`).
3. Hard-required `prod_ai` column missing from 2013's FAERS DRUG table (`build_multi_year_drug_vocabulary.py`).
4. `rebuild_global_downstream.py`'s shared `annual_summary_all_years.json` blindly overwritten instead of merged.
5. `build_gold_annual_pair_reactions.py` staged and swapped the *entire* Gold output root — would have deleted 2023-2025's Gold partitions the first time `--overwrite` was used for an additive year range.

Full test suite: **95/95 passing.**

### Changed/new files this session (uncommitted)

`src/integrate_global_case_assignments.py`, `src/faers_utils.py`,
`src/build_multi_year_drug_vocabulary.py`, `src/run_ms7c_models.py`,
`src/rebuild_global_downstream.py`, `src/build_gold_annual_pair_reactions.py`,
`src/run_historical_baseline.py` (new), `docs/historical_baseline_v1.md` (new),
`docs/stage4_historical_expansion.md`, `README.md`,
`tests/test_faers_utils.py` (new), `tests/test_run_historical_baseline.py` (new),
`tests/test_ms4b_global_assignment_integration.py`,
`tests/test_build_multi_year_drug_vocabulary.py`,
`tests/test_full_global_downstream_rebuild.py`,
`tests/test_build_gold_annual_pair_reactions.py`.

New data/report artifacts (all additive, gitignored except where noted):
`data/interim/case_resolution/` (promoted), `data/interim/case_resolution_pre_stage4_backup/`,
`data/processed/faers_case_level_annual/2013..2022/`,
`data/interim/drug_mapping/multi_year_drug_vocabulary_2013_2022_norm_v1.parquet`,
`data/processed/mappings/faers_to_chembl_2013_2022_chembl36_norm_v1_exact_v1.parquet`,
`data/processed/faers_case_level_annual_final/2013..2022/`,
`data/gold/annual_pair_reactions/year=2013..2022/`,
`data/gold/temporal_pair_reaction_features_historical_v1/`,
`data/processed/modeling/historical_baseline_v1/`,
`reports/data_quality/historical_2013_2022/`, `reports/modeling/historical_baseline_v1/`.

### Known pre-existing issue (unrelated to this session's code)

`data/processed/faers_case_level_annual_final/2024/` and `2025/`'s
`drug_pair_reaction_counts.parquet`, and `data/gold/temporal_pair_reaction_features/`
(the original pilot's feature directory), are currently Access-Denied at the
OS level (even `icacls` fails) — same symptom class as the
`data/gold/annual_pair_reactions/` ACL issue fixed earlier this session, but
this one is still open. Not caused by anything run this session (verified by
code inspection — no command this session ever targeted those years/paths)
and does not block any completed work, but worth a look.

## Next up

1. MS8A: port the expensive data-engineering stages to PySpark.
2. MS8B: Spark/local correctness and output-parity validation.
3. MS8C: scaling benchmark and engine recommendation.

Full spec: [docs/post_ms7_temporal_modeling_requirements.md](docs/post_ms7_temporal_modeling_requirements.md)
and [docs/ms8_pyspark_confluence.md](docs/ms8_pyspark_confluence.md).
