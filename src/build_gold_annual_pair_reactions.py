"""Publish the canonical MS7A annual pair-reaction Gold dataset.

The builder consumes validated MS6 annual aggregates, validates their strict
analytical contract in bounded batches, and atomically publishes Hive-style
year partitions. Feature engineering is deliberately outside this milestone.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_ROOT = ROOT / "data/processed/faers_case_level_annual_final"
DEFAULT_OUTPUT_ROOT = ROOT / "data/gold/annual_pair_reactions"
DEFAULT_ASSIGNMENTS = ROOT / "data/interim/case_resolution/global_case_assignments.parquet"
DEFAULT_REPORT = ROOT / "reports/data_quality/ms7_gold_validation.json"
GOLD_COLUMNS = ["year", "drug_a", "drug_b", "reaction_pt", "n_cases"]
CHEMBL_PATTERN = re.compile(r"^CHEMBL[0-9]+$")
GOLD_SCHEMA = pa.schema([
    pa.field("year", pa.int32(), nullable=False),
    pa.field("drug_a", pa.string(), nullable=False),
    pa.field("drug_b", pa.string(), nullable=False),
    pa.field("reaction_pt", pa.string(), nullable=False),
    pa.field("n_cases", pa.int64(), nullable=False),
])


def repository_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return path.name


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
    if not years:
        raise ValueError("at least one year is required")
    return sorted(years)


def validate_global_assignments(path: Path, years: list[int]) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"global assignment artifact missing: {path}")
    table = pq.read_table(path, columns=["caseid", "assigned_report_year", "assignment_status"])
    frame = table.to_pandas()
    included = frame.loc[frame["assignment_status"].eq("included")].copy()
    if included["caseid"].isna().any() or included["assigned_report_year"].isna().any():
        raise ValueError("included global assignments contain null CASEID or assigned year")
    conflicts = included.groupby("caseid", sort=False)["assigned_report_year"].nunique().gt(1)
    if conflicts.any():
        raise ValueError(f"{int(conflicts.sum())} CASEIDs contribute to more than one assigned year")
    requested = included["assigned_report_year"].astype(int).isin(years)
    return {
        "included_assignment_rows": int(len(included)),
        "included_assignment_rows_for_requested_years": int(requested.sum()),
        "caseids_assigned_to_multiple_years": 0,
    }


def _validate_batch(frame: pd.DataFrame, year: int) -> None:
    missing = set(GOLD_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(f"{year}: source is missing columns {sorted(missing)}")
    if frame[GOLD_COLUMNS].isna().any().any():
        raise ValueError(f"{year}: Gold columns contain null values")
    if not frame["year"].eq(year).all():
        raise ValueError(f"{year}: source contains a different analytical year")
    if not frame["drug_a"].lt(frame["drug_b"]).all():
        raise ValueError(f"{year}: drug pairs are not strictly ordered")
    if not frame["n_cases"].ge(1).all():
        raise ValueError(f"{year}: n_cases must be at least 1")
    valid_a = frame["drug_a"].astype(str).str.fullmatch(CHEMBL_PATTERN)
    valid_b = frame["drug_b"].astype(str).str.fullmatch(CHEMBL_PATTERN)
    if not (valid_a & valid_b).all():
        raise ValueError(f"{year}: invalid ChEMBL identifier")


def publish_year(source: Path, destination: Path, year: int, batch_size: int) -> dict:
    if not source.is_file():
        raise FileNotFoundError(f"{year}: MS6 aggregate missing: {source}")
    parquet = pq.ParquetFile(source)
    if not set(GOLD_COLUMNS).issubset(parquet.schema_arrow.names):
        missing = set(GOLD_COLUMNS) - set(parquet.schema_arrow.names)
        raise ValueError(f"{year}: source is missing columns {sorted(missing)}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = pq.ParquetWriter(destination, GOLD_SCHEMA, compression="snappy")
    rows = 0
    previous_key: tuple | None = None
    try:
        for batch in parquet.iter_batches(batch_size=batch_size, columns=GOLD_COLUMNS):
            frame = batch.to_pandas()
            _validate_batch(frame, year)
            keys = list(frame[["year", "drug_a", "drug_b", "reaction_pt"]].itertuples(index=False, name=None))
            if len(keys) != len(set(keys)):
                raise ValueError(f"{year}: duplicate Gold key within a batch")
            if keys:
                if previous_key is not None and keys[0] <= previous_key:
                    raise ValueError(f"{year}: duplicate or non-deterministically ordered Gold key")
                if any(left >= right for left, right in zip(keys, keys[1:])):
                    raise ValueError(f"{year}: duplicate or non-deterministically ordered Gold key")
                previous_key = keys[-1]
            normalized = pa.Table.from_pandas(frame[GOLD_COLUMNS], schema=GOLD_SCHEMA, preserve_index=False)
            writer.write_table(normalized)
            rows += len(frame)
    finally:
        writer.close()
    if rows == 0:
        raise ValueError(f"{year}: Gold partition would be empty")
    return {"year": year, "rows": rows, "source": repository_path(source), "partition": f"year={year}"}


def build_gold(
    years: list[int], input_root: Path, output_root: Path, assignments: Path,
    report_path: Path, batch_size: int, overwrite: bool,
) -> dict:
    if output_root.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {output_root}; pass --overwrite to replace it")
    assignment_result = validate_global_assignments(assignments, years)
    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_root.name}-", dir=output_root.parent))
    results = []
    try:
        for year in years:
            partition = staging / f"year={year}"
            partition.mkdir()
            result = publish_year(
                input_root / str(year) / "drug_pair_reaction_counts.parquet",
                partition / "part-00000.parquet", year, batch_size,
            )
            results.append(result)
        metadata = {
            "milestone": "MS7A",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "schema": GOLD_COLUMNS,
            "partitioning": "year",
            "years": years,
            "total_rows": sum(item["rows"] for item in results),
            "partitions": results,
            "global_assignment_validation": assignment_result,
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
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the validated MS7A partitioned annual pair-reaction Gold dataset.")
    parser.add_argument("--years", nargs="+", required=True, help="Years or inclusive ranges, for example 2023-2025")
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--assignments", type=Path, default=DEFAULT_ASSIGNMENTS)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--batch-size", type=int, default=500_000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error("--batch-size must be positive")
    metadata = build_gold(
        parse_years(args.years), args.input_root, args.output_root,
        args.assignments, args.report, args.batch_size, args.overwrite,
    )
    print(f"MS7A Gold dataset published: {metadata['total_rows']:,} rows across {len(metadata['years'])} year(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
