# MS8A/MS8B: PySpark implementation and parity

Implementation of the spec in [ms8_pyspark_confluence.md](ms8_pyspark_confluence.md).
Status on 2026-10-12: MS8A implemented, MS8B parity validated for 2013-2025,
MS8C (scaling benchmark) pending.

## What runs in Spark

`src/spark_pair_reactions.py` rebuilds the MS7A Gold contract
`(year, drug_a, drug_b, reaction_pt, n_cases)` from the same inputs as
`rebuild_global_downstream.py`: the MS4B resolved DRUG/REAC rows and the
per-period MS5 exact ChEMBL mapping (2013-2022 and 2023-2025 artifacts).

| Stage | Spark implementation |
|---|---|
| Reaction terms | `normalize_reaction_pt` (shared with the local pipeline) on distinct `pt` values via `mapInArrow`; drop blanks and `REACTION_PT_EXCLUDE`; distinct `(caseid, reaction_pt)` |
| Drug identity | `_choose_raw_drug_name` + `normalize_drug_string` on distinct `(prod_ai, drugname)` values via `mapInArrow`; broadcast join to the mapping; only `mapping_status == "mapped"`; distinct `(caseid, chembl_id)` |
| Pairs | self-join on `caseid` with `drug_a < drug_b` (string order, as in the local pipeline) |
| Support | join to reactions on `caseid`, `countDistinct(caseid)` per `(drug_a, drug_b, reaction_pt)` |
| Output | each task writes its own Parquet part file with pyarrow (`year=YYYY/part-NNNNN.parquet`) |

Only the default `all_pairs` role filter is implemented; it is the filter used
by every production year.

Design decisions:

- **Shared string semantics.** With pandas 3, the local pipeline's string
  methods run on Arrow kernels, where regex `\s` does not match `\xa0`. A
  hand-written Python copy would disagree on rare terms, so Spark calls the
  local functions themselves. `tests/test_spark_pair_reactions.py` pins this.
- **Normalize distinct values, not rows.** The drug normalizer runs once per
  distinct name pair instead of once per DRUG row. Output is identical, but
  part of Spark's speedup over the local pipeline comes from this change, not
  from parallelism.
- **Per-year jobs.** Every case belongs to one analytical year and every key
  contains the year, so per-year processing equals one multi-year job while
  bounding shuffle size on a single machine.
- **No Hadoop writer.** Spark's local-filesystem writer needs `winutils.exe` on
  Windows, and collecting a year to the driver exceeds
  `spark.driver.maxResultSize` from 2016 on. Distributed pyarrow writes avoid both.
- **pandas UDFs are not used.** They crash with pandas 3 under PySpark 4.2;
  `mapInArrow` does not go through PySpark's pandas serializer.

## Running

Requires Java 17+ (`JAVA_HOME`) and `pip install -r requirements.txt`.

```powershell
# parity against production Gold (never writes under data/gold)
.venv\Scripts\python.exe src\spark_pair_reactions.py --years 2013-2025 --cores 8 `
  --check-against data\gold\annual_pair_reactions

# benchmark rows (one CSV row per year); --noop times compute without writing
.venv\Scripts\python.exe src\spark_pair_reactions.py --years 2013-2025 --cores 1 `
  --out-dir data\interim\spark\bench --benchmark-csv reports\spark\benchmark.csv
```

Default output: `data/interim/spark/annual_pair_reactions/` (gitignored).
Defaults: `--driver-memory 6g`, shuffle partitions `max(8, 4 * cores)`.

## MS8B parity result

[`reports/spark/parity_2013_2025.json`](../reports/spark/parity_2013_2025.json):
all 13 years are identical to the production Gold partitions:
256,643,925 rows, 0 keys only in Spark, 0 keys only in Gold, 0 `n_cases`
mismatches, 0 duplicate keys. The run used `--cores 8` and took 1,294 s.
2013 was also identical at `--cores 1`, so output does not depend on
parallelism or partition count.

Machine and versions: [`reports/spark/environment.json`](../reports/spark/environment.json).
About 7 GB of commit memory was free during the runs, so a 12g heap failed to
allocate. Keep the heap at 6g on this machine, or close other applications.

## MS8C scaling benchmark

```powershell
.venv\Scripts\python.exe src\run_spark_benchmark.py --dry-run   # show the 42 configs
.venv\Scripts\python.exe src\run_spark_benchmark.py             # about 3-4 h unattended
.venv\Scripts\python.exe src\analyze_spark_benchmark.py         # tables + figures
```

| Phase | Configs | Output |
|---|---|---|
| A, scaling grid | 2013-2025 at `--cores` 1, 2, 4, 8, 16, 24 | `reports/spark/benchmark_scaling.csv` |
| B, repeats | 2023 x 3 repeats per core count, alternating order | `reports/spark/benchmark_repeats.csv` |
| C, compute only | same as B with `--noop` | `reports/spark/benchmark_noop.csv` |

Every config runs in a fresh Spark session. The runner keeps Windows awake,
continues past a failed config, and logs exit code, wall time and free memory
per config to `reports/spark/benchmark_runs.csv`. Full logs go to
`data/interim/spark/bench_logs/`. Close memory-heavy applications first;
free memory affects timings.

`analyze_spark_benchmark.py` produces:
- `benchmark_summary.json`: speedup, efficiency, Amdahl serial fraction,
  repeat mean ± SD, write share
- `benchmark_by_year.csv`: local pipeline vs Spark at 1 core vs the best core
  count, with join size per year
- `fig_spark_speedup` and `fig_spark_runtime_by_year` in `reports/figures/paper/`

Local pipeline reference timings come from
`reports/data_quality/**/ms6_stage_performance_by_year.csv`. They include
validation and intermediate artifact writes, so the local-vs-Spark ratio mixes
engine, algorithm (distinct-value normalization) and I/O effects. The
parallelism effect alone is Spark at 1 core vs Spark at N cores.
