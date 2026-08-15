# Milestone 8 - PySpark implementation and benchmark

Milestone 8 is delivered in three ordered sub-milestones:

| Sub-milestone | Name | Outcome |
|---|---|---|
| MS8A | PySpark pipeline implementation | Expensive FAERS transformations run end to end in PySpark |
| MS8B | Correctness and parity validation | Spark outputs match the validated local pipeline |
| MS8C | Performance benchmark and recommendation | Increasing data ranges are measured and an engine recommendation is documented |

| Field | Value |
|---|---|
| Status | Planned |
| Prerequisites | MS7A Gold, MS7B temporal features, and MS7C pilot models |
| Pilot data | FAERS 2023-2025 |
| Historical benchmark data | FAERS 2013-2025 after Stage 4 expansion |
| Primary comparison | Existing local Pandas/PyArrow/DuckDB pipeline vs PySpark |
| Output | Validated Spark datasets, parity reports, and scaling benchmarks |

## Objective

Implement the computationally expensive FAERS data-engineering stages in
PySpark and determine when distributed processing becomes more effective than
the existing local pipeline.

This milestone must demonstrate both correctness and scalability. A faster
Spark result is not acceptable if it changes case selection, drug pairs,
reaction counts, temporal features, or other analytical semantics.

## Why this milestone exists

The local pipeline is already capable of producing validated datasets for the
2023-2025 pilot. PySpark is introduced to evaluate larger historical ranges,
parallelize expensive joins and aggregations, and document the scale at which
Spark's distributed execution offsets its startup, serialization, shuffle,
and coordination overhead.

The goal is not to rewrite every small script in Spark. The goal is to port
the stages where data volume, joins, grouping, or distinct counting make a
distributed engine useful.

## Scope

The PySpark implementation must include:

1. Multi-file FAERS ingestion.
2. Global case filtering and selected-version application.
3. Drug-pair generation.
4. Pair-reaction join.
5. Distinct contributing-case calculation with `countDistinct(caseid)`.
6. Annual pair-reaction aggregation.
7. Lag-feature generation from annual Gold data.
8. Output-parity validation against the existing local implementation.
9. Performance benchmarks over increasing year ranges.

## Out of scope

- Replacing the validated ChEMBL mapping policy.
- Changing the global CASEID assignment rules.
- Redefining the MS7A Gold schema.
- Redefining the MS7B feature or target semantics.
- Training only the final small model table in Spark and calling that a Spark
  pipeline implementation.
- Distributed model training unless separately justified after the data
  engineering benchmark.
- Streaming ingestion or real-time FAERS processing.
- Production cloud deployment, cluster autoscaling, or managed platform
  selection unless added as a separate milestone.

## Input contracts

The Spark pipeline must consume the same logical inputs as the local pipeline:

- official FAERS quarterly archives or their reproducibly converted tabular
  representation;
- global CASEID assignments and selected PRIMARYIDs;
- versioned FAERS-to-ChEMBL mapping;
- MS7A annual Gold partitions for feature-generation parity; and
- the same configuration values, exclusion rules, role-filter policy, and
  analytical year boundaries.

Input versions and hashes must be recorded for every benchmark run.

## Required Spark stages

### MS8A - PySpark pipeline implementation

MS8A ports the expensive data-engineering stages. It is complete only when the
2023-2025 pilot can run end to end in PySpark and publish outputs conforming to
the existing analytical contracts.

#### 1. Multi-file ingestion

Read multiple quarterly and annual partitions with an explicit schema. Avoid
schema inference during timed benchmark runs. Preserve source year, quarter,
table, and physical-file provenance needed for validation.

#### 2. Global case filtering

Apply the globally selected CASEID/PRIMARYID boundary. A case must contribute
to no more than one assigned analytical year, and nonselected follow-up
versions must not re-enter downstream processing.

#### 3. Drug-pair generation

Build deterministic unordered pairs from uniquely mapped ChEMBL identifiers.
Require `drug_a < drug_b`, reject self-pairs, and preserve the configured pair
role policy.

#### 4. Pair-reaction join

Join selected drug pairs to deduplicated case-level reaction terms. Document
the chosen partitioning and join strategy, including whether either side is
broadcast and why.

#### 5. Distinct case support

Aggregate by:

```text
year, drug_a, drug_b, reaction_pt
```

and calculate:

```text
countDistinct(caseid) AS n_cases
```

Raw joined row counts must never be substituted for distinct case support.

#### 6. Annual aggregation

Publish the same five-column MS7A contract:

```text
year
drug_a
drug_b
reaction_pt
n_cases
```

Write horizontally partitioned Parquet by year and preserve deterministic
analytical results regardless of Spark partition count.

#### 7. Lag-feature generation

Reproduce the MS7B feature contract using only history before each target
year. Unavailable history must remain distinct from an observed year with zero
support. Candidate-population and cold-start rules must match the local
implementation exactly.

## MS8B - Correctness and parity validation

For every compared year range, Spark and the local pipeline must agree on:

- assigned CASEID population;
- retained selected PRIMARYIDs;
- mapped case-drug rows;
- unique case-level drug pairs;
- unique case-pair-reaction rows;
- annual `(year, drug_a, drug_b, reaction_pt)` keys;
- `n_cases` for every annual key;
- temporal feature keys and values;
- recurrence targets; and
- total, per-year, and stage-level row counts.

Validation must compare complete outputs where practical. If a full comparison
is temporarily infeasible, use deterministic hashes plus targeted row-level
samples, and document the limitation. Zero unexplained mismatches are allowed.

## MS8C - Performance benchmark and recommendation

MS8C begins only after MS8B establishes analytical parity for the pilot.

### Benchmark plan

Run both implementations over increasing ranges after the corresponding local
outputs exist:

| Benchmark | Range | Purpose |
|---|---|---|
| B1 | 2024 | Single complete year and Spark overhead baseline |
| B2 | 2023-2025 | Existing three-year pilot and cross-year behavior |
| B3 | 2020-2025 | Medium historical range after Stage 4 |
| B4 | 2016-2025 | Larger historical range after Stage 4 |
| B5 | 2013-2025 | Full intended historical benchmark |

The exact intermediate ranges may be adjusted to available hardware, but the
benchmark must contain at least three increasing data sizes and include the
2023-2025 pilot.

### Metrics to capture

For every engine, range, and stage, record:

- input files and bytes;
- input and output row counts;
- wall-clock duration;
- peak process or executor memory where available;
- temporary disk and shuffle bytes;
- CPU or executor utilization where available;
- Spark executor, core, memory, and partition configuration;
- number and size of output files;
- failed or retried tasks;
- code version and input hashes; and
- parity-validation status.

Do not compare only total runtime. Stage-level measurements are required to
identify where Spark helps and where it adds overhead.

### Fair-comparison rules

- Both engines must implement the same analytical rules.
- Both must read equivalent input data and write equivalent logical outputs.
- Timed runs must exclude one-time dependency installation and source-data
  download.
- Clearly label cold-cache and warm-cache measurements.
- Run each benchmark repeatedly and report the median; retain individual run
  measurements.
- Record hardware and storage characteristics.
- Spark configuration must be fixed or explicitly versioned between runs.
- Performance tuning must not weaken validation or change output semantics.

## Required artifacts

```text
src/spark/
├── ingest_faers.py
├── build_case_boundary.py
├── build_pair_reactions.py
├── build_annual_gold.py
└── build_temporal_features.py

config/spark/
└── benchmark_profiles.yaml

reports/benchmark/ms8/
├── benchmark_runs.csv
├── stage_performance.csv
├── output_parity.csv
├── spark_configuration.json
└── ms8_summary.json
```

Large Spark outputs, event logs, shuffle data, and temporary benchmark files
must remain ignored generated artifacts.

## Validation requirements

- Explicit schemas are used for production reads.
- No CASEID contributes to more than one assigned year.
- Pair ordering and identifier invariants match MS7A.
- Annual support uses distinct CASEIDs.
- Spark Gold keys and counts match local Gold outputs.
- Spark temporal features and targets match MS7B outputs.
- Output schemas and nullability match their contracts.
- Repeated runs with different safe Spark partition counts produce equivalent
  analytical results.
- Every benchmark has complete provenance and configuration metadata.
- Temporary files are cleaned after successful runs.

## Acceptance criteria

### MS8A acceptance

MS8A is complete when:

1. All required expensive stages have working PySpark implementations.
2. The 2023-2025 pilot completes end to end in Spark.
3. Explicit schemas, configuration, provenance, and output contracts are
   implemented and tested.

### MS8B acceptance

MS8B is complete when:

1. Spark and local outputs have zero unexplained analytical mismatches.
2. Gold and temporal-feature parity checks pass for the pilot.
3. Repeated safe partition configurations produce equivalent analytical
   results.
4. Automated tests cover Spark transformations and parity contracts.

### MS8C acceptance

MS8C is complete when:

1. At least three increasing data ranges are benchmarked, subject to Stage 4
   data availability.
2. Benchmark runs are reproducible from checked-in code and configuration.
3. Stage-level runtime and resource results are documented.
4. The final report states where Spark becomes beneficial, remains slower, or
   is constrained by the available hardware.

Milestone 8 is complete only when MS8A, MS8B, and MS8C are complete.

The combined acceptance conditions are:

1. All required expensive stages have working PySpark implementations.
2. The 2023-2025 pilot completes end to end in Spark.
3. Spark and local outputs have zero unexplained analytical mismatches.
4. At least three increasing data ranges are benchmarked, subject to Stage 4
   data availability.
5. Benchmark runs are reproducible from checked-in code and configuration.
6. Stage-level runtime and resource results are documented.
7. The final report states where Spark becomes beneficial, remains slower, or
   is constrained by the available hardware.
8. Automated tests cover transformation logic and output contracts.

Spark does not need to outperform the local pipeline at every scale. A valid
outcome may show that local execution is preferable for the pilot while Spark
becomes useful only for larger historical ranges.

## Risks and limitations

- A single-machine Spark benchmark is not equivalent to a distributed cluster
  benchmark.
- Spark startup and shuffle overhead may dominate small ranges.
- Distinct counting and pair-reaction joins may require substantial shuffle
  storage.
- Different Parquet file layouts can affect performance even when logical data
  is identical.
- Historical FAERS schema differences may affect Stage 4 ingestion before they
  affect Spark itself.
- Hardware limits may prevent the full 2013-2025 comparison; any reduced scope
  must be documented rather than silently omitted.

## Deliverable summary

MS8A-MS8C deliver a correctness-equivalent PySpark implementation of the
expensive FAERS stages, validated parity, and an evidence-based comparison with
the existing local pipeline. The conclusion should recommend which engine to
use at each tested scale and explain that recommendation using measured results.
