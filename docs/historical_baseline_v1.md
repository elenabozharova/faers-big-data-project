# Historical Baseline v1 — 2013-2025 rolling temporal evaluation

## Purpose

The 2023-2025 pilot (MS7C) could not test whether deeper reporting history
helps predict drug-pair-reaction recurrence, because no pilot target year had
three full years of real lag history. Historical Baseline v1 promotes the
validated Stage 4 MS4A result, builds 2013-2022 additively through MS7A
(2023-2025 untouched), rebuilds MS7B across the full 2013-2025 range, and
reruns the existing MS7C baselines and logistic regression across every
target year that now has three complete lags (2016-2025), for three feature
sets (M1/M2/M3), to test that question directly.

## Scope decision

2023-2025's MS5 (vocabulary + ChEMBL mapping), MS6, and MS7A artifacts are
**untouched** -- verified by SHA-256 (the 2023-2025 mapping file's checksum
still matches the one published in the README). Only 2013-2022 got new,
additively versioned MS5/MS6/MS7A artifacts. This is safe because
`build_faers_chembl_mapping.py`'s `exact_component()` maps a normalized name
against the ChEMBL 36 index alone, independent of which years fed the
vocabulary -- adding 2013-2022 names cannot change any 2023-2025 mapping
outcome.

## Feature sets

| Set | Features |
|---|---|
| M1 | `count_t1` |
| M2 | `count_t1`, `count_t2`, `count_t3` |
| M3 | `count_t1`, `count_t2`, `count_t3`, `present_t1..3`, `historical_total`, `years_observed`, `years_since_last_seen` (existing MS7C support set) |

All logistic-regression features are standardized (`sklearn.StandardScaler`,
zero mean / unit variance) per training year before fitting
(`run_ms7c_models.train_logistic`), so coefficients are in standardized space
and directly comparable across features and across years. They are **not**
comparable to raw, unscaled feature units. For every 2016-2025 target year,
`count_t1/t2/t3` and `present_t1/t2/t3` are always observed (the full
2013-2025 history covers all three lags), so their `_unavailable` indicator
columns are constant zero and correctly receive a zero coefficient -- this is
expected, not a bug.

## Results

### Per-year metrics (M3, the full support feature set)

| Target year | Train year | n_samples | Positive rate | Support-rule F1 | Logistic F1 | Logistic AUROC | Logistic AUPRC |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2016 | 2015 | 30,420,807 | 0.2033 | 0.4523 | 0.5196 | 0.7494 | 0.5269 |
| 2017 | 2016 | 42,368,465 | 0.1747 | 0.4389 | 0.5117 | 0.7685 | 0.5338 |
| 2018 | 2017 | 52,716,536 | 0.1700 | 0.4190 | 0.5169 | 0.7830 | 0.5515 |
| 2019 | 2018 | 63,679,998 | 0.1650 | 0.4482 | 0.5210 | 0.7993 | 0.5621 |
| 2020 | 2019 | 74,305,569 | 0.1683 | 0.4537 | 0.5399 | 0.8065 | 0.5733 |
| 2021 | 2020 | 86,277,409 | 0.1282 | 0.3953 | 0.4839 | 0.8180 | 0.5370 |
| 2022 | 2021 | 95,221,750 | 0.1252 | 0.4431 | 0.4965 | 0.8266 | 0.5501 |
| 2023 | 2022 | 104,532,085 | 0.1284 | 0.4935 | 0.5190 | 0.8388 | 0.5805 |
| 2024 | 2023 | 114,369,768 | 0.1351 | 0.5276 | 0.4998 | 0.8490 | 0.6179 |
| 2025 | 2024 | 124,265,716 | 0.1411 | 0.5285 | 0.5149 | 0.8529 | 0.6385 |

Logistic regression (M3) beats the simple support-threshold rule on F1 in
**8 of 10 years**; the rule wins in 2024 and 2025. AUROC and AUPRC rise
steadily as more historical years become available for later target years
(2016: AUROC 0.749 -> 2025: AUROC 0.853), consistent with deeper history
improving the model's ranking quality.

### Ablation: M1 vs M2 vs M3 (mean +/- SD across the 10 target years)

| Feature set | Mean F1 | SD F1 | Mean AUROC | SD AUROC | Mean AUPRC | SD AUPRC |
|---|---:|---:|---:|---:|---:|---:|
| M1 (`count_t1` only) | 0.4424 | 0.0398 | 0.7291 | 0.0248 | 0.4379 | 0.0425 |
| M2 (3 annual lags) | 0.5200 | 0.0255 | 0.7850 | 0.0260 | 0.5284 | 0.0330 |
| M3 (full support set) | 0.5123 | 0.0156 | **0.8092** | 0.0346 | **0.5672** | 0.0367 |

M2 improves every metric over M1 by a clear margin (F1 +0.078, AUROC +0.056,
AUPRC +0.090) -- the three-year-lag features carry real, additional signal
beyond the most recent year alone. M3 has the highest AUROC and AUPRC (both
threshold-independent ranking metrics) but a very slightly lower F1 than M2
at the fixed 0.5 decision threshold; this looks like a threshold-calibration
artifact from the extra features shifting the probability distribution, not
worse ranking quality -- the AUROC/AUPRC gains say M3 is still the strongest
model at separating positives from negatives.

### `count_t1` coefficient dominance

M2 (standardized coefficients, mean across 10 years):

| Feature | Mean coefficient | Rank-1 years (out of 10) |
|---|---:|---:|
| `count_t1` | 5.393 | 10/10 |
| `count_t2` | 2.240 | 0/10 |
| `count_t3` | 1.472 | 0/10 |

M3 (standardized coefficients, mean across 10 years, non-zero features):

| Feature | Mean coefficient | Rank-1 years (out of 10) |
|---|---:|---:|
| `count_t1` | 2.524 | 10/10 |
| `years_observed` | 0.580 | 0/10 |
| `historical_total` | 0.509 | 0/10 |
| `count_t2` | 0.496 | 0/10 |
| `years_since_last_seen` | -0.301 | 0/10 |
| `count_t3` | 0.169 | 0/10 |
| `present_t1` | 0.026 | 0/10 |
| `present_t2` | -0.008 | 0/10 |
| `present_t3` | -0.009 | 0/10 |

The recency ordering `|coef(count_t1)| > |coef(count_t2)| > |coef(count_t3)|`
holds in **10/10** years in the M3 model. `count_t1` is the single strongest
feature in every year in both M2 and M3 -- that part of the pilot's finding
holds up. But it does not dominate to the exclusion of everything else:
`count_t2` and `count_t3` have clearly non-zero, substantial coefficients
(not noise-level), and M2's mean F1/AUROC/AUPRC are all meaningfully higher
than M1's. The pilot's apparent "count_t1 completely dominates" reading was
an artifact of its training data: with only one prior year (2023) available
for the train_year=2024 partition, `count_t2`, `count_t3`,
`years_observed`, `years_since_last_seen`, and `present_t1..3` were all
*constant* (zero-variance) columns in that specific training set -- the
model was structurally unable to learn from them, not evidence they carry no
signal.

## Scientific interpretation

**Does the full 2013-2025 history support the pilot's conclusion that the
most recent year's support dominates recurrence prediction?**
Partially. `count_t1` is confirmed as the single strongest predictor in
every one of the 10 years tested (M2 and M3 coefficient ranks). But "dominates"
in the sense of "the rest doesn't matter" is not supported: M2 clearly
outperforms M1.

**Do older temporal lags and cumulative history provide measurable predictive
value beyond `count_t1`?** Yes. M2's mean F1 is 17.5% relatively higher than
M1's (0.520 vs 0.442), with a similarly-sized AUROC and AUPRC gain, and this
holds consistently across all 10 target years (M2 beats M1 on every metric
in every single year -- see `ablation_per_year_wide.csv`). M3's additional
features (`historical_total`, `years_observed`, `years_since_last_seen`)
further improve ranking quality (AUROC, AUPRC) beyond M2, though not F1 at a
fixed threshold.

**Conclusion**: *count_t1 remains the strongest single predictor, but the
three-year pilot overstated its exclusivity -- multi-year history and
cumulative-support features consistently add real, non-trivial predictive
value once the data actually contains enough history to show it.*

## Temporal anomalies (observed, not causally asserted)

- Positive rate (prevalence) drifts from ~20% (2016) down to ~13% (2021-2023)
  before rising again to ~14% (2025) -- an observed pattern. A plausible but
  unverified explanation is that the growing candidate population (more years
  of history contribute more rarely-recurring historical pairs as the window
  widens) dilutes prevalence; this was not causally tested.
- 2021 has the lowest support-rule and logistic F1 of the 10 years
  (F1 0.395/0.484) and the lowest positive rate alongside 2022-2023. This
  coincides with the COVID-19 pandemic period, a plausible explanation given
  documented FAERS reporting-volume shifts in 2020-2021, but this analysis
  did not isolate reporting-volume effects from other confounds, so it is
  reported as a coincidence worth further investigation, not a demonstrated
  cause.
- Sample size grows roughly monotonically from 30.4M (2016) to 124.3M (2025)
  rows, tracking the growing candidate population as more historical years
  become available -- expected by construction, not an anomaly.

## Known limitations / deviations from the pilot

- **Candidate populations are not directly comparable across the pilot and
  this run.** The pilot's target_year=2025 partition used `available_years
  =[2023,2024,2025]`; this run's target_year=2025 partition uses
  `available_years=2013-2025`. A tuple observed only once, a decade ago, now
  qualifies as a candidate where it would not have under the pilot's narrower
  window. This is very likely why this run's target_year=2025 numbers
  (support F1 0.5285, logistic F1 0.5149) are lower than the original pilot's
  (support F1 0.6298, logistic F1 0.6176) for what is nominally the same
  train/test year pair -- a larger, more diverse, harder candidate
  population, not a regression in model quality. Direct pilot-vs-historical
  number comparisons should account for this.
- The known, expected +/-1-2 case MS4B discrepancy for 2023-2025 documented in
  `docs/stage4_historical_expansion.md` (from the wider 2013-2025 global
  index revealing slightly earlier `INIT_FDA_DT` values for a small number of
  cases) is **not** corrected here, since 2023-2025's MS4B/MS6/MS7A were
  intentionally left unrebuilt.
- Gradient-boosted trees remain optional and were not implemented, per the
  existing MS7C scope.
- Three genuine pre-existing pipeline bugs were found and fixed while
  building this dataset (see the final report's "Problems/fixes" section);
  each has a regression test.

## Reproduction

```powershell
# Phase 1: promote MS4A (already done; not re-run)
# Phase 2: MS5 (2013-2022 only)
.venv\Scripts\python.exe src\build_multi_year_drug_vocabulary.py `
  --years 2013-2022 --annual-input-root data\processed\faers_case_level_annual `
  --input-filename resolved_drug_rows.parquet `
  --output data\interim\drug_mapping\multi_year_drug_vocabulary_2013_2022_norm_v1.parquet `
  --report-dir reports\data_quality\historical_2013_2022 `
  --normalization-version norm_v1 --artifact-version vocab_v1_hist_2013_2022 `
  --sample-raw-names-limit 10

.venv\Scripts\python.exe src\build_faers_chembl_mapping.py `
  --years 2013-2022 `
  --vocabulary data\interim\drug_mapping\multi_year_drug_vocabulary_2013_2022_norm_v1.parquet `
  --name-index data\interim\drug_mapping\chembl36_name_index_norm_v1_index_v1.parquet `
  --output data\processed\mappings\faers_to_chembl_2013_2022_chembl36_norm_v1_exact_v1.parquet `
  --annual-input-root data\processed\faers_case_level_annual `
  --report-dir reports\data_quality\historical_2013_2022 `
  --normalization-version norm_v1 --vocabulary-artifact-version vocab_v1_hist_2013_2022 `
  --mapping-artifact-version exact_v1_hist_2013_2022 --name-index-artifact-version index_v1

# MS6 (2013-2022 only)
.venv\Scripts\python.exe src\rebuild_global_downstream.py `
  --years 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 `
  --global-case-assignments data\interim\case_resolution\global_case_assignments.parquet `
  --mapping-file data\processed\mappings\faers_to_chembl_2013_2022_chembl36_norm_v1_exact_v1.parquet `
  --report-dir reports\data_quality\historical_2013_2022

# MS7A (2013-2022 only, additive)
.venv\Scripts\python.exe src\build_gold_annual_pair_reactions.py `
  --years 2013-2022 `
  --report reports\data_quality\historical_2013_2022\ms7_gold_validation.json

# MS7B (full 2013-2025 history; target years must include one year before the
# earliest evaluation year, i.e. 2015, so it can serve as the 2016 training partition)
.venv\Scripts\python.exe src\build_temporal_features.py `
  --available-years 2013-2025 --target-years 2015-2025 `
  --output-root data\gold\temporal_pair_reaction_features_historical_v1 `
  --report reports\data_quality\historical_2013_2022\ms7b_temporal_features_validation.json `
  --memory-limit 22GB --threads 2

# Historical Baseline v1 (30 model runs: 10 target years x M1/M2/M3)
.venv\Scripts\python.exe src\run_historical_baseline.py `
  --feature-root data\gold\temporal_pair_reaction_features_historical_v1 `
  --target-years 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 `
  --batch-size 500000
```

Outputs: `reports/modeling/historical_baseline_v1/` (per-year metrics, full
detailed metrics, coefficients, coefficient dominance, ablation tables, and
`historical_baseline_v1_manifest.json`, which also records the source commit,
data/target years, feature and model definitions, and environment details for
reproduction).
