"""MS8A: PySpark build of the annual (year, drug_a, drug_b, reaction_pt, n_cases) Gold contract.

Reads the same immutable MS4B boundary (resolved DRUG/REAC rows) and the same
versioned MS5 ChEMBL mapping as rebuild_global_downstream.py, and reproduces its
semantics with the default `all_pairs` role filter:

  reactions: pt -> lowercase/strip/collapse whitespace, drop blanks and
             REACTION_PT_EXCLUDE, one row per (caseid, reaction_pt)
  drugs:     prod_ai if meaningful else drugname -> frozen normalizer ->
             exact mapping; only mapping_status == "mapped" survives;
             one row per (caseid, chembl_id)
  pairs:     every unordered pair per case with drug_a < drug_b (string order)
  support:   countDistinct(caseid) per (drug_a, drug_b, reaction_pt)

Years are processed one at a time: every case belongs to exactly one analytical
year and every key contains the year, so per-year processing is semantically
identical to one multi-year job while bounding shuffle size on a single node.

String normalization reuses the local pipeline's own pandas code inside
mapInArrow (pandas 3 runs these on Arrow string kernels, whose regex semantics
differ from Python's re). It runs once per distinct (prod_ai, drugname) and
distinct pt value instead of once per row; the result is identical, but part of
the speedup over the local pipeline is this deduplication rather than parallelism.

Each Spark task writes its own Parquet part file with pyarrow (mapInArrow),
because Hadoop's local-filesystem writer needs winutils.exe on Windows. Output
never passes through the driver, so it is not bounded by driver memory.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

SRC = Path(__file__).resolve().parent
ROOT = SRC.parent
# Spark's Python workers must import the frozen normalizer and use this interpreter.
os.environ["PYTHONPATH"] = os.pathsep.join(filter(None, [str(SRC), os.environ.get("PYTHONPATH")]))
os.environ.setdefault("PYSPARK_PYTHON", sys.executable)

import pandas as pd  # noqa: E402

import process_faers_case_level as pipeline  # noqa: E402
from faers_utils import REACTION_PT_EXCLUDE  # noqa: E402
from normalize_drug_names import normalize_drug_string  # noqa: E402

DEFAULT_BOUNDARY_ROOT = ROOT / "data/processed/faers_case_level_annual"
DEFAULT_OUT_DIR = ROOT / "data/interim/spark/annual_pair_reactions"
MAPPINGS = ROOT / "data/processed/mappings"
GOLD_SCHEMA = pa.schema([
    pa.field("year", pa.int32(), nullable=False),
    pa.field("drug_a", pa.string(), nullable=False),
    pa.field("drug_b", pa.string(), nullable=False),
    pa.field("reaction_pt", pa.string(), nullable=False),
    pa.field("n_cases", pa.int64(), nullable=False),
])


def mapping_for_year(year: int) -> Path:
    """Same per-period artifacts that produced the production Gold partitions."""
    period = "2013_2022" if year <= 2022 else "2023_2025"
    return MAPPINGS / f"faers_to_chembl_{period}_chembl36_norm_v1_exact_v1.parquet"


def drug_name_batches(batches):
    """mapInArrow: (prod_ai, drugname) -> normalized_drug_name via the local pipeline's own code."""
    for batch in batches:
        frame = batch.to_pandas().astype("str")  # same dtype pd.read_parquet gives the local pipeline
        chosen = pipeline._choose_raw_drug_name(frame)
        blank = chosen.str.strip() == ""
        normalized = pd.Series(None, index=frame.index, dtype=object)
        normalized[~blank] = chosen[~blank].apply(normalize_drug_string)
        yield pa.RecordBatch.from_arrays(
            [batch.column(0), batch.column(1), pa.array(normalized, type=pa.string(), from_pandas=True)],
            names=["prod_ai", "drugname", "normalized_drug_name"])


def reaction_pt_batches(batches):
    """mapInArrow: pt -> reaction_pt via process_faers_case_level.normalize_reaction_pt."""
    for batch in batches:
        reaction_pt = pipeline.normalize_reaction_pt(batch.column(0).to_pandas().astype("str"))
        yield pa.RecordBatch.from_arrays(
            [batch.column(0), pa.array(reaction_pt, type=pa.string(), from_pandas=True)], names=["pt", "reaction_pt"])


def build_spark(cores: int, shuffle_partitions: int, driver_memory: str):
    from pyspark.sql import SparkSession

    return (
        SparkSession.builder.master(f"local[{cores}]").appName("faers-ms8a")
        .config("spark.driver.memory", driver_memory)
        .config("spark.sql.shuffle.partitions", str(shuffle_partitions))
        .config("spark.ui.enabled", "false")
        .getOrCreate()
    )


def year_counts(spark, year: int, boundary_root: Path):
    """Lazy DataFrame of one year's Gold rows."""
    from pyspark.sql import functions as F

    caseid = F.nullif(F.trim(F.col("caseid")), F.lit(""))
    source = boundary_root / str(year)

    drug = spark.read.parquet(str(source / "resolved_drug_rows.parquet"))
    prod_ai = F.col("prod_ai") if "prod_ai" in drug.columns else F.lit(None).cast("string")
    drug = drug.select(
        caseid.alias("caseid"),
        F.coalesce(prod_ai, F.lit("")).alias("prod_ai"),
        F.coalesce(F.col("drugname"), F.lit("")).alias("drugname"),
    ).where(F.col("caseid").isNotNull())

    names = drug.select("prod_ai", "drugname").distinct().mapInArrow(
        drug_name_batches, "prod_ai string, drugname string, normalized_drug_name string")
    mapping = (
        spark.read.parquet(str(mapping_for_year(year)))
        .where(F.col("mapping_status") == "mapped")
        .select("normalized_drug_name", F.col("chembl_id").alias("drug"))
    )
    case_drugs = (
        drug.join(names, ["prod_ai", "drugname"])
        .join(F.broadcast(mapping), "normalized_drug_name")
        .select("caseid", "drug").distinct()
    )
    a, b = case_drugs.alias("a"), case_drugs.alias("b")
    pairs = a.join(b, (F.col("a.caseid") == F.col("b.caseid")) & (F.col("a.drug") < F.col("b.drug"))).select(
        F.col("a.caseid").alias("caseid"), F.col("a.drug").alias("drug_a"), F.col("b.drug").alias("drug_b"))

    reac = spark.read.parquet(str(source / "resolved_reac_rows.parquet")).select(
        caseid.alias("caseid"), F.coalesce(F.col("pt"), F.lit("")).alias("pt"))
    pts = reac.select("pt").distinct().mapInArrow(reaction_pt_batches, "pt string, reaction_pt string")
    reactions = (
        reac.where(F.col("caseid").isNotNull())
        .join(F.broadcast(pts), "pt")
        .where(F.col("reaction_pt").isNotNull() & ~F.col("reaction_pt").isin(*sorted(REACTION_PT_EXCLUDE)))
        .select("caseid", "reaction_pt").distinct()
    )
    return (
        pairs.join(reactions, "caseid")
        .groupBy("drug_a", "drug_b", "reaction_pt")
        .agg(F.countDistinct("caseid").alias("n_cases"))
        .select(F.lit(year).cast("int").alias("year"), "drug_a", "drug_b", "reaction_pt", "n_cases")
    )


def write_year(counts, out_dir: Path, year: int) -> int:
    """Each Spark task writes its own part file with pyarrow (distributed write, no driver collect)."""
    year_dir = out_dir / f"year={year}"
    if year_dir.exists():
        shutil.rmtree(year_dir)  # stale part files from a run with a different partition count
    year_dir.mkdir(parents=True)

    def write_partition(batches):
        from pyspark import TaskContext

        context = TaskContext.get()
        destination = year_dir / f"part-{context.partitionId():05d}.parquet"
        staged = destination.with_name(f"{destination.name}.{context.attemptNumber()}.tmp")
        writer, rows = None, 0
        for batch in batches:
            if writer is None:
                writer = pq.ParquetWriter(staged, GOLD_SCHEMA)
            writer.write_batch(batch.cast(GOLD_SCHEMA))
            rows += batch.num_rows
        if writer is not None:
            writer.close()
            os.replace(staged, destination)
        yield pa.RecordBatch.from_pydict({"rows": [rows]})

    return counts.mapInArrow(write_partition, "rows long").groupBy().sum("rows").first()[0] or 0


def parity(spark_dir: Path, gold_root: Path, year: int) -> dict:
    """Exact key and n_cases comparison against the local Gold partition (DuckDB, a third engine)."""
    import duckdb

    gold_file = gold_root / f"year={year}" / "part-00000.parquet"
    key = "year, drug_a, drug_b, reaction_pt"
    s, g = f"read_parquet('{(spark_dir / '*.parquet').as_posix()}')", f"read_parquet('{gold_file.as_posix()}')"
    con = duckdb.connect()
    one = lambda sql: con.sql(sql).fetchone()[0]  # noqa: E731
    result = {
        "year": year,
        "spark_rows": one(f"select count(*) from {s}"),
        "gold_rows": one(f"select count(*) from {g}"),
        "spark_duplicate_keys": one(f"select count(*) - count(distinct ({key})) from {s}"),
        "keys_only_in_spark": one(f"select count(*) from (select {key} from {s} except select {key} from {g})"),
        "keys_only_in_gold": one(f"select count(*) from (select {key} from {g} except select {key} from {s})"),
        "n_cases_mismatches": one(f"select count(*) from {s} x join {g} y using ({key}) where x.n_cases <> y.n_cases"),
    }
    result["identical"] = (result["spark_rows"] == result["gold_rows"] and result["spark_duplicate_keys"] == 0
                           and result["keys_only_in_spark"] == result["keys_only_in_gold"]
                           == result["n_cases_mismatches"] == 0)
    return result


def append_csv(path: Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        if new:
            writer.writeheader()
        writer.writerow(row)


def parse_years(values: list[str]) -> list[int]:
    years = set()
    for value in values:
        start, _, end = value.partition("-")
        years.update(range(int(start), int(end or start) + 1))
    return sorted(years)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--years", nargs="+", required=True, help="e.g. 2023 or 2013-2025")
    parser.add_argument("--cores", type=int, default=1, help="local[N]; start at 1 per the benchmark design")
    parser.add_argument("--shuffle-partitions", type=int, help="default max(8, 4 * cores)")
    parser.add_argument("--driver-memory", default="6g")
    parser.add_argument("--boundary-root", type=Path, default=DEFAULT_BOUNDARY_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--noop", action="store_true", help="execute the full plan without writing output (compute-only timing)")
    parser.add_argument("--check-against", type=Path, help="Gold root to compare each written year against")
    parser.add_argument("--benchmark-csv", type=Path, help="append one timing row per year")
    args = parser.parse_args()

    years = parse_years(args.years)
    gold = (ROOT / "data/gold").resolve()
    if args.out_dir.resolve().is_relative_to(gold):
        raise ValueError("refusing to write under data/gold: year folders in --out-dir are replaced")
    shuffle_partitions = args.shuffle_partitions or max(8, 4 * args.cores)

    started = time.perf_counter()
    spark = build_spark(args.cores, shuffle_partitions, args.driver_memory)
    session_seconds = time.perf_counter() - started
    results = []
    try:
        for year in years:
            year_started = time.perf_counter()
            counts = year_counts(spark, year, args.boundary_root)
            row = {"year": year, "cores": args.cores, "shuffle_partitions": shuffle_partitions,
                   "mode": "noop" if args.noop else "write"}
            if args.noop:
                counts.write.format("noop").mode("overwrite").save()
                row["output_rows"] = ""
            else:
                row["output_rows"] = write_year(counts, args.out_dir, year)
            row["seconds"] = round(time.perf_counter() - year_started, 3)
            if args.check_against and not args.noop:
                row["parity"] = parity(args.out_dir / f"year={year}", args.check_against, year)
            row["session_start_seconds"] = round(session_seconds, 3)
            row["finished_at"] = datetime.now(timezone.utc).isoformat()
            print(json.dumps(row), flush=True)
            if args.benchmark_csv:
                append_csv(args.benchmark_csv, {k: v for k, v in row.items() if k != "parity"})
            results.append(row)
    finally:
        spark.stop()
    failed = [r["year"] for r in results if "parity" in r and not r["parity"]["identical"]]
    if failed:
        raise SystemExit(f"parity FAILED for years {failed}")


if __name__ == "__main__":
    main()
