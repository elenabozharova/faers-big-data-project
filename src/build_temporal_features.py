"""Build leakage-safe MS8 temporal recurrence features from MS7 Gold Parquet."""

from __future__ import annotations

import argparse
import json
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import duckdb


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GOLD_ROOT = ROOT / "data/gold/annual_pair_reactions"
DEFAULT_OUTPUT_ROOT = ROOT / "data/gold/temporal_pair_reaction_features"
DEFAULT_REPORT = ROOT / "reports/data_quality/ms8_temporal_features_validation.json"
KEYS = ["drug_a", "drug_b", "reaction_pt"]
FEATURE_COLUMNS = [
    "target_year", *KEYS, "count_t1", "count_t2", "count_t3",
    "present_t1", "present_t2", "present_t3", "historical_total",
    "years_observed", "years_since_last_seen", "target",
]


def parse_years(values: list[str]) -> list[int]:
    years: set[int] = set()
    for value in values:
        if "-" in value:
            start, end = (int(part) for part in value.split("-", 1))
            if start > end:
                raise ValueError(f"invalid year range: {value}")
            years.update(range(start, end + 1))
        else:
            years.add(int(value))
    return sorted(years)


def _sql_path(path: Path) -> str:
    return path.resolve().as_posix().replace("'", "''")


def _source_file(gold_root: Path, year: int) -> Path:
    path = gold_root / f"year={year}" / "part-00000.parquet"
    if not path.is_file():
        raise FileNotFoundError(f"MS7 Gold partition missing for {year}: {path}")
    return path


def feature_query(history_files: list[Path], target_file: Path, target_year: int, first_year: int) -> str:
    history = ", ".join(f"'{_sql_path(path)}'" for path in history_files)
    target = _sql_path(target_file)

    def lag(name: str, offset: int, present: bool = False) -> str:
        year = target_year - offset
        if year < first_year:
            return f"CAST(NULL AS TINYINT) AS {name}" if present else f"CAST(NULL AS BIGINT) AS {name}"
        expression = f"COALESCE(MAX(CASE WHEN year = {year} THEN n_cases END), 0)"
        if present:
            expression = f"CAST(({expression}) > 0 AS TINYINT)"
        return f"{expression} AS {name}"

    return f"""
        WITH history AS (
            SELECT year, drug_a, drug_b, reaction_pt, n_cases
            FROM read_parquet([{history}], union_by_name = true)
            WHERE year < {target_year}
        ),
        aggregated AS (
            SELECT
                CAST({target_year} AS INTEGER) AS target_year,
                drug_a, drug_b, reaction_pt,
                {lag('count_t1', 1)},
                {lag('count_t2', 2)},
                {lag('count_t3', 3)},
                {lag('present_t1', 1, True)},
                {lag('present_t2', 2, True)},
                {lag('present_t3', 3, True)},
                CAST(SUM(n_cases) AS BIGINT) AS historical_total,
                CAST(COUNT(*) AS INTEGER) AS years_observed,
                CAST({target_year} - MAX(year) AS INTEGER) AS years_since_last_seen
            FROM history
            GROUP BY drug_a, drug_b, reaction_pt
        ),
        observed_target AS (
            SELECT drug_a, drug_b, reaction_pt
            FROM read_parquet('{target}')
        )
        SELECT
            a.*,
            CAST(t.drug_a IS NOT NULL AS TINYINT) AS target
        FROM aggregated a
        LEFT JOIN observed_target t USING (drug_a, drug_b, reaction_pt)
        ORDER BY drug_a, drug_b, reaction_pt
    """


def validate_partition(connection: duckdb.DuckDBPyConnection, path: Path, target_year: int) -> dict:
    source = _sql_path(path)
    columns = [row[0] for row in connection.execute(f"DESCRIBE SELECT * FROM read_parquet('{source}')").fetchall()]
    if columns != FEATURE_COLUMNS:
        raise ValueError(f"{target_year}: unexpected feature schema: {columns}")
    result = connection.execute(f"""
        SELECT
            COUNT(*) AS rows,
            COUNT(*) - COUNT(DISTINCT (drug_a, drug_b, reaction_pt)) AS duplicate_keys,
            COUNT(*) FILTER (WHERE target_year <> {target_year}) AS wrong_year,
            COUNT(*) FILTER (WHERE drug_a >= drug_b OR drug_a !~ '^CHEMBL[0-9]+$' OR drug_b !~ '^CHEMBL[0-9]+$') AS invalid_pairs,
            COUNT(*) FILTER (WHERE historical_total < 1 OR years_observed < 1 OR years_since_last_seen < 1) AS invalid_history,
            COUNT(*) FILTER (WHERE target NOT IN (0, 1)) AS invalid_targets,
            SUM(target) AS positive_targets
        FROM read_parquet('{source}')
    """).fetchone()
    names = [item[0] for item in connection.description]
    metrics = dict(zip(names, result))
    for name in ("duplicate_keys", "wrong_year", "invalid_pairs", "invalid_history", "invalid_targets"):
        if metrics[name] != 0:
            raise ValueError(f"{target_year}: {name}={metrics[name]}")
    if metrics["rows"] == 0:
        raise ValueError(f"{target_year}: empty temporal feature partition")
    return {key: int(value or 0) for key, value in metrics.items()}


def build_features(
    available_years: list[int], target_years: list[int], gold_root: Path,
    output_root: Path, report_path: Path, overwrite: bool,
    memory_limit: str = "8GB", threads: int = 4,
) -> dict:
    available_years = sorted(set(available_years))
    target_years = sorted(set(target_years))
    if not available_years or not target_years:
        raise ValueError("available and target years are required")
    if any(target not in available_years for target in target_years):
        raise ValueError("every target year must have an available Gold partition")
    if output_root.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {output_root}; pass --overwrite to replace it")
    sources = {year: _source_file(gold_root, year) for year in available_years}
    for target in target_years:
        if not any(year < target for year in available_years):
            raise ValueError(f"{target}: no history is available before target year")

    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_root.name}-", dir=output_root.parent))
    database = staging / "ms8_work.duckdb"
    connection = duckdb.connect(str(database))
    connection.execute(f"SET memory_limit = '{memory_limit}'")
    connection.execute(f"SET threads = {int(threads)}")
    results = []
    try:
        for target in target_years:
            history = [sources[year] for year in available_years if year < target]
            partition = staging / f"target_year={target}"
            partition.mkdir()
            destination = partition / "part-00000.parquet"
            query = feature_query(history, sources[target], target, min(available_years))
            connection.execute(
                f"COPY ({query}) TO '{_sql_path(destination)}' "
                "(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880)"
            )
            metrics = validate_partition(connection, destination, target)
            metrics["target_year"] = target
            metrics["history_years"] = [year for year in available_years if year < target]
            results.append(metrics)
        connection.close()
        database.unlink(missing_ok=True)
        metadata = {
            "milestone": "MS8",
            "status": "stage_3_pilot",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "candidate_population": "tuples observed at least once before target_year",
            "unavailable_lag_policy": "null; distinct from observed zero support",
            "available_years": available_years,
            "target_years": target_years,
            "schema": FEATURE_COLUMNS,
            "partitions": results,
            "total_rows": sum(item["rows"] for item in results),
            "all_invariants_passed": True,
        }
        (staging / "_SUCCESS").touch()
        (staging / "_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        if output_root.exists():
            shutil.rmtree(output_root)
        staging.replace(output_root)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        return metadata
    except Exception:
        connection.close()
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description="Build MS8 leakage-safe temporal recurrence features from MS7 Gold data.")
    parser.add_argument("--available-years", nargs="+", required=True)
    parser.add_argument("--target-years", nargs="+", required=True)
    parser.add_argument("--gold-root", type=Path, default=DEFAULT_GOLD_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--memory-limit", default="8GB")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    result = build_features(
        parse_years(args.available_years), parse_years(args.target_years),
        args.gold_root, args.output_root, args.report, args.overwrite,
        args.memory_limit, args.threads,
    )
    print(f"MS8 pilot features published: {result['total_rows']:,} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
