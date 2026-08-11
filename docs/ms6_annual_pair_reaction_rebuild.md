# MS6 final annual pair-reaction rebuild

## Scope and milestone boundary

The full downstream rebuild was deferred at the MS4B retained/resolved
boundary until the final MS5 mapping existed. One production run now serves
two distinct conclusions:

- **MS4B closure:** prove that the globally selected CASEID/PRIMARYID
  population propagates correctly through every downstream annual artifact.
- **MS6 production:** create final 2023-2025 ChEMBL-mapped annual
  drug-pair-reaction support datasets.

MS4B is complete and production-validated end to end. MS6 is complete and
production-validated for 2023-2025. MS7A was not started: this work creates no
lag features, temporal labels, models, or Spark implementation.

## Operational strategy

`src/rebuild_global_downstream.py` reuses these immutable MS4B inputs:

```text
data/processed/faers_case_level_annual/<year>/retained_cases.parquet
data/processed/faers_case_level_annual/<year>/resolved_drug_rows.parquet
data/processed/faers_case_level_annual/<year>/resolved_reac_rows.parquet
```

The run does not reread or reintegrate raw FAERS tables. It validates each
boundary against `data/interim/case_resolution/global_case_assignments.parquet`,
creates same-volume hard links for the three reused files, stages the complete
year in a hidden sibling directory, and publishes it atomically only after all
year-level invariants pass. Existing annual outputs require `--overwrite`.

Authoritative outputs are:

```text
data/processed/faers_case_level_annual_final/2023/
data/processed/faers_case_level_annual_final/2024/
data/processed/faers_case_level_annual_final/2025/
```

## Mapping and downstream flow

The production mapping is:

```text
data/processed/mappings/faers_to_chembl_2023_2025_chembl36_norm_v1_exact_v1.parquet
SHA-256 e4e28c86dea312c262ce9ab92c006430817d6c1dbfb0a348a33e1dab39d408a5
```

Only `mapping_status = mapped` with a nonblank single ChEMBL ID contributes to
`case_drugs`. Ordinary ambiguity, unmatched names, fully recognized component
combinations, partial combinations, and ambiguous combinations remain visible
in `case_drugs_long` but contribute no ordinary case drug or pair. Fuzzy
mappings are disabled.

The unchanged downstream methodology is:

```text
resolved REAC -> normalized/excluded/deduplicated case reactions
resolved DRUG -> frozen name selection/normalization -> MS5 exact mapping
mapped rows -> one row per CASEID/ChEMBL ID -> deterministic unordered pairs
pairs -> role audit -> all_pairs selection
selected pairs + case reactions -> bounded pair-reaction join
unique case-pair-reaction rows -> distinct-CASEID support -> min3/min5/min10
```

The pair-reaction join used the existing 2,000,000-row chunk cap. Support is
the number of distinct logical CASEIDs for `(drug_a, drug_b, reaction_pt)`, not
raw DRUG or REAC multiplicity.

## Production counts

| Year | Retained | Resolved DRUG | Resolved REAC | Case reactions | Case drugs | Pair rows | Pair-reaction rows |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2023 | 1,361,395 | 4,763,980 | 3,860,778 | 3,510,978 | 3,176,711 | 11,134,568 | 133,383,154 |
| 2024 | 1,308,773 | 5,008,214 | 3,791,696 | 3,450,588 | 3,168,844 | 13,313,358 | 259,927,804 |
| 2025 | 1,295,857 | 4,942,013 | 3,711,750 | 3,399,331 | 3,235,435 | 14,122,469 | 209,952,427 |

`case_drugs_long` has one row per resolved DRUG row. Pair-role audit and
selected-pair row counts equal the pair counts above because the frozen
production policy is `all_pairs`.

| Year | Unique patterns | Min3 | Min5 | Min10 |
|---:|---:|---:|---:|---:|
| 2023 | 23,256,862 | 6,393,185 | 4,234,914 | 2,390,729 |
| 2024 | 25,346,924 | 7,067,742 | 4,789,551 | 3,012,592 |
| 2025 | 28,871,255 | 8,434,213 | 5,658,578 | 3,326,862 |

## Case mapping usability

| Year | Zero mapped | One mapped | At least two mapped | Combination-only | Ambiguous-only |
|---:|---:|---:|---:|---:|---:|
| 2023 | 108,081 | 803,902 | 449,412 | 7 | 3,736 |
| 2024 | 78,081 | 804,854 | 425,838 | 11 | 3,054 |
| 2025 | 64,701 | 803,749 | 427,407 | 7 | 4,460 |

Counts use unique ChEMBL IDs per case. Duplicate DRUG rows never create extra
drugs. “Combination-only” and “ambiguous-only” mean that the case has no
accepted single-ID mapping and has at least one corresponding explicit audit
status; these categories are audit flags and need not be mutually exclusive
with other nonmapped statuses. Every retained case had at least one DRUG row,
and no case with DRUG rows had only blank/invalid normalized names.

## MS4B validation

The retained counts match global assignments exactly. Resolved DRUG/REAC and
all downstream case populations are annual retained-case subsets. No
nonselected PRIMARYID reached a resolved file, no unexpected CASEID or
PRIMARYID was observed, and transitivity from disjoint retained populations
proves downstream cross-year uniqueness.

The audited 2024 CASEIDs `24346610`, `24346611`, and `24346612` remain retained
and have resolved DRUG rows. They are absent from resolved REAC, case reactions,
pair reactions, and support counts, exactly as required.

## Performance and storage

| Year | Runtime seconds | Join chunks | Final logical file size |
|---:|---:|---:|---:|
| 2023 | 1,042.911 | 69 | 466.94 MB |
| 2024 | 1,381.920 | 135 | 501.10 MB |
| 2025 | 1,460.473 | 110 | 508.42 MB |

Case-drug mapping dominated each year at 653.663, 665.977, and 804.940
seconds. The pair-reaction joins took 138.367, 310.820, and 295.463 seconds.
The independently observed process peak was 17,960.96 MB. The Windows memory
helper was corrected after the run; the stage report records this conservative
run-level observed peak for every stage rather than inventing stage-specific
measurements. Case-drug deduplication is performed inside mapping, so its
separate timing row is zero and its row counts—not an additive duration—are
the useful metric.

## Reports and legacy comparison

- `reports/data_quality/ms4b_downstream_validation_by_year.csv`
- `reports/data_quality/ms4b_downstream_validation_summary.json`
- `reports/data_quality/ms6_case_mapping_usability_by_year.csv`
- `reports/data_quality/ms6_annual_pair_reaction_validation.csv`
- `reports/data_quality/ms6_annual_pair_reaction_summary.json`
- `reports/data_quality/ms6_stage_performance_by_year.csv`
- `reports/data_quality/ms6_production_input_integrity.json`
- `reports/data_quality/ms4b_ms6_legacy_global_downstream_comparison.csv`

Only earlier annual retained-case totals were available for legacy comparison.
Legacy reaction, mapped-drug, pair, pair-reaction, and support artifacts do not
exist locally, so those report fields explicitly say `unavailable`; no values
were invented. Differences from a hypothetical legacy rebuild would combine
global CASEID/PRIMARYID correction with the conservative final MS5 mapping and
cannot be separated without generating another full counterfactual pipeline.

## Integrity and validation

Before and after production, SHA-256 hashes were computed for all 12 raw ZIPs,
global assignments, the final mapping, and all nine reused boundary Parquets.
Every hash was unchanged. No raw, assignment, mapping, or resolved input was
rewritten. No staging directory, `.part`, or temporary SQLite file remained,
and persisted report paths are repository-relative POSIX strings.

The complete final suite passed 77 unit tests. Python compilation, the annual
and rebuild CLI help commands, `pip check`, and `git diff --check` also passed;
the latter reported only the repository's existing LF-to-CRLF warnings.

Production command:

```powershell
python src/rebuild_global_downstream.py --years 2023 2024 2025
```

Validation commands:

```powershell
python -m compileall src tests
python -m unittest discover -s tests -p "test_*.py"
python src/process_faers_case_level_annual.py --help
python src/rebuild_global_downstream.py --help
python -m pip check
git diff --check
git status --short
```

## Limitations and next boundary

The run uses exact ChEMBL single-ID mappings only and does not establish manual
mapping accuracy. Combination components are not split into case drugs. The
production role policy remains `all_pairs`; `non_cc` remains a sensitivity
mode. No legacy full downstream artifacts exist for an effect-separated
comparison. Peak memory is a process-level observation rather than a
stage-isolated profiler measurement.

MS7A remains pending. Temporal Gold construction, lag features, train/test
splits, models, and Spark are outside this milestone.
