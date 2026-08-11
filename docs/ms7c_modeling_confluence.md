# MS7C - Temporal recurrence baselines and models

| Field | Value |
|---|---|
| Status | Planned |
| Owner | Data Science / Data Engineering |
| Input | MS7B temporal recurrence feature dataset |
| Pilot years | Features from 2023-2024; targets for 2024-2025 |
| Output | Reproducible model evaluations, predictions, and selected thresholds |
| Out of scope | Historical 2013-2022 expansion and MS8 PySpark implementation |

## Objective

Determine how well prior FAERS reporting history predicts whether a known
`(drug_a, drug_b, reaction_pt)` tuple will be observed again in a target year.
MS7C must establish simple, interpretable baselines before fitting statistical
or machine-learning models.

## Business question

For a tuple observed at least once before target year `t`, can information
available strictly before `t` predict whether that tuple appears in `t`?

This is recurrence prediction for previously observed tuples. It is not a
cold-start model for tuples first appearing in the target year, and it does
not establish causality or clinical risk.

## Input contract

The model phase consumes the partitioned MS7B feature dataset under:

```text
data/gold/temporal_pair_reaction_features/
├── target_year=2024/
└── target_year=2025/
```

Required columns are:

```text
target_year
drug_a
drug_b
reaction_pt
count_t1
count_t2
count_t3
present_t1
present_t2
present_t3
historical_total
years_observed
years_since_last_seen
target
```

Unavailable lags are null and must remain distinguishable from an available
year with zero support. Model preprocessing must be fitted using training data
only and recorded as part of each model artifact.

## Required implementation order

Models must be implemented and evaluated in this order. A later item must not
replace or obscure the results of an earlier baseline.

### 1. Always-negative baseline

Predict `target = 0` for every eligible row. This establishes the accuracy and
specificity obtainable without identifying any recurrence.

### 2. Previous-year persistence rule

Predict recurrence when `present_t1 = 1`. Rows without an available `t-1` lag
are not silently converted to negative; their handling must be reported.

### 3. Support-threshold rule

Predict recurrence when `count_t1 >= k`. Select `k` using the validation
period only. Report results for the evaluated threshold grid and preserve the
selected value in machine-readable metadata.

### 4. Logistic regression

Fit an interpretable binary logistic-regression model using the numeric and
binary MS7B features. Document null handling, scaling, class weighting, solver,
regularization, convergence, coefficient direction, and decision threshold.
Hyperparameters and the classification threshold must be selected without
using the final test labels.

### 5. Optional gradient-boosted trees

Implement only after the preceding four approaches are complete. It must use
the same candidate population, chronological split, and primary metrics so its
performance is directly comparable. The additional complexity is justified
only if it produces a meaningful improvement over logistic regression.

## Temporal evaluation design

Random row splitting is prohibited for the primary evaluation because rows
from different target years represent a time sequence.

The Stage 3 pilot has only target years 2024 and 2025. Therefore:

- 2024 may be used for fitting and threshold development;
- 2025 is the held-out forward test year; and
- results are preliminary because the pilot has no independent middle
  validation year.

Any threshold or hyperparameter selection performed on 2024 must be explicit.
The 2025 labels may be used only for final pilot evaluation. Stronger
train/validation/test separation becomes mandatory after Stage 4 adds
historical target years.

## Metrics

The recurrence target may be imbalanced, so accuracy alone is insufficient.
Every approach must report:

- eligible row count and positive prevalence;
- confusion matrix: TP, FP, TN, FN;
- precision;
- recall;
- F1 score;
- specificity;
- balanced accuracy; and
- predicted-positive rate.

Probability-producing models must additionally report ROC AUC, average
precision / PR AUC, and Brier score. Metrics must be reported by target year
and for any explicitly defined aggregate view.

## Required artifacts

```text
data/processed/modeling/ms7c/
├── predictions/
│   └── target_year=<year>/
├── models/
└── metadata/

reports/modeling/ms7c/
├── model_comparison.csv
├── metrics_by_target_year.csv
├── confusion_matrices.csv
├── support_threshold_search.csv
└── ms7c_summary.json
```

Large predictions and serialized model files remain ignored generated
artifacts. Small validation reports and documentation may be committed.

## Validation requirements

- Input schema exactly matches the documented MS7B contract.
- One prediction exists for every eligible input row for each evaluated model.
- Prediction keys are unique by
  `(model_id, target_year, drug_a, drug_b, reaction_pt)`.
- Predictions and targets are binary for rule-based classifiers.
- Probabilities are finite and within `[0, 1]` where applicable.
- No feature uses the target year or any later year.
- Preprocessing and threshold selection do not use held-out test labels.
- All models use the same documented candidate population.
- Recomputed metrics match the persisted predictions.
- Runs record code version, parameters, input hashes, row counts, and runtime.

## Acceptance criteria

MS7C is complete for the Stage 3 pilot when:

1. The first four required approaches are implemented in order.
2. Each approach has reproducible predictions and validation reports.
3. The 2025 forward-test results are reported without training leakage.
4. Results are compared on the same eligible population and metrics.
5. Logistic-regression coefficients and preprocessing are documented.
6. All validation requirements pass through automated tests.
7. Limitations of the two-target-year pilot are stated prominently.

Gradient-boosted trees are optional and do not block MS7C completion.

## Risks and limitations

- The Stage 3 pilot offers only one fitting year and one forward-test year.
- Pilot rows have incomplete three-year lag histories.
- Candidate generation excludes cold-start tuples first observed in the target
  year, so results apply only to recurrence among previously observed tuples.
- FAERS reporting patterns do not measure incidence or prove drug-event
  causality.
- Strong performance may reflect persistence in reporting behavior rather
  than a clinically meaningful safety signal.

## Milestone boundary and next steps

MS7C covers local baseline and model evaluation only. After the pilot pipeline
works end to end, Stage 4 expands the data through 2013-2022 and enables a
stronger multi-year temporal evaluation. MS8 then ports expensive data
engineering stages to PySpark and benchmarks the local pipeline against Spark.
