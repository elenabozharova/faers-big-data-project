# MS7B/MS7C temporal modeling and MS8A/MS8B/MS8C PySpark requirements

This document defines the required work after the canonical MS7A Gold
dataset. MS7B and MS7C are implemented and production-validated for the Stage
3 pilot. Stage 4 historical expansion and MS8A, MS8B, and MS8C remain pending.

## Staged data rollout

The project deliberately validates the pipeline in progressively larger data
stages:

### Stage 2 - One complete year

Use all four quarters of 2024 (`2024 Q1-Q4`) to exercise the existing annual
wrapper and validate a complete single-year run.

### Stage 3 - Three-year pilot

Use all quarters from 2023-2025. This pilot is intended to exercise:

- multiple annual outputs;
- repeated CASEIDs across years;
- temporal recurrence;
- drug-mapping coverage; and
- preliminary feature generation.

MS2 through MS7A are production-validated for this pilot window. MS7B uses
it to implement and validate preliminary temporal logic. Within this window,
2024 can have one prior-year lag and 2025 can have two prior-year lags. A
three-prior-year row is not available unless an earlier year is added or 2026
becomes available, so unavailable lags must be represented and documented
explicitly rather than silently treated as observed zero counts.

### Stage 4 - Historical expansion

After the 2023-2025 pipeline and preliminary features work end to end, ingest
and process 2013-2022 through the same validated milestone contracts. Rebuild
the multi-year mapping, annual MS6 aggregates, and MS7A Gold partitions for the
expanded range before producing the full historical feature and modeling
datasets.

## MS7B - Temporal features and recurrence target

For every target year `t`, features must use observations strictly before `t`.
The binary target is derived only from whether the tuple is observed in `t`:

```text
history before t -> features
observed in t     -> target
```

The initial feature contract is:

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

`count_t1` is the support in `t-1`, with `count_t2` and `count_t3` defined
analogously. The `present_*` columns are their binary presence indicators.
`historical_total`, `years_observed`, and `years_since_last_seen` must be
calculated only from years before the target year. Feature construction and
dataset splitting must prevent future-year leakage.

The current 2023-2025 Gold dataset is the Stage 3 pilot for preliminary feature
generation. A complete three-lag feature row requires either a 2026 target or
Gold data through at least `t-3`. Stage 4 therefore expands the history through
2013-2022 before the full historical three-lag feature dataset is declared
production-complete. Pilot rows with incomplete history must distinguish an
unavailable lag from an observed year with zero support.

### Stage 3 production result

The MS7B feature builder uses DuckDB for bounded-memory Parquet aggregation and
publishes Hive-style partitions under
`data/gold/temporal_pair_reaction_features/`.

| Target year | Candidate rows | Positive targets | Available history |
|---:|---:|---:|---|
| 2024 | 23,256,862 | 9,614,927 | 2023 |
| 2025 | 38,988,859 | 13,070,201 | 2023-2024 |
| **Total** | **62,245,721** | **22,685,128** | |

All schema, uniqueness, temporal-boundary, pair-ordering, identifier, history,
and binary-target invariants passed. The candidate population is the union of
tuples observed before each target year. Target-only cold starts are excluded
because they have no historical evidence from which to predict recurrence.

Reproduce the pilot build with:

```powershell
python src/build_temporal_features.py `
  --available-years 2023-2025 `
  --target-years 2024-2025
```

## MS7C - Baselines and statistical model

The complete milestone definition, acceptance criteria, evaluation protocol,
and artifact contract are maintained in the
[MS7C Confluence-ready specification](ms7c_modeling_confluence.md).

Models must be implemented and evaluated in this order:

1. Always-negative baseline.
2. Previous-year persistence rule.
3. Support-threshold rule.
4. Logistic regression.
5. Optional gradient-boosted trees.

Evaluation must use chronological train, validation, and test boundaries.
Random row splitting is not acceptable for the primary temporal evaluation.
Model comparisons must use the same eligible population and report metrics
appropriate for a likely imbalanced recurrence target.

## MS8A/MS8B/MS8C - PySpark implementation and benchmark

The complete scope, parity contract, benchmark design, artifacts, and
acceptance criteria are maintained in the
[MS8A-MS8C Confluence-ready specification](ms8_pyspark_confluence.md).

Spark work must port the expensive data-engineering stages, not merely train a
model from an already-small feature table. The required Spark scope is:

- multi-file ingestion;
- global case filtering;
- pair-reaction joins;
- `countDistinct(caseid)` support calculation;
- annual aggregation; and
- lag-feature generation.

Pandas/PyArrow and Spark implementations must be checked for output parity and
benchmarked over increasing year ranges. Benchmarks must record input size,
runtime, peak memory where available, Spark configuration, partition counts,
and output row counts. The comparison should identify the scale at which
Spark's distributed execution offsets its startup and shuffle overhead.
