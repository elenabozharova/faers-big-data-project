"""Build a deterministic vocabulary from globally corrected annual DRUG rows."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sqlite3
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from download_faers import sha256_file
from integrate_global_case_assignments import normalize_identifier, repository_path
from normalize_drug_names import normalize_drug_string
from process_faers_case_level import _choose_raw_drug_name


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_ROOT = ROOT / "data/processed/faers_case_level_annual"
DEFAULT_ASSIGNMENTS = ROOT / "data/interim/case_resolution/global_case_assignments.parquet"
DEFAULT_OUTPUT_DIR = ROOT / "data/interim/drug_mapping"
DEFAULT_REPORT_DIR = ROOT / "reports/data_quality"
REQUIRED_COLUMNS = ("caseid", "primaryid", "prod_ai", "drugname")
INVALID_NORMALIZED_SENTINELS = frozenset({"nan", "none", "null"})


@dataclass(frozen=True)
class AnnualDrugInput:
    analytical_year: int
    path: Path


def parse_years(values: list[str]) -> list[int]:
    if not values:
        raise ValueError("at least one analytical year is required")
    years: list[int] = []
    for token in values:
        token = token.strip()
        if not token:
            raise ValueError("year tokens must not be blank")
        if "-" in token:
            parts = token.split("-")
            if len(parts) != 2 or not all(part.isdigit() for part in parts):
                raise ValueError(f"invalid year range: {token!r}")
            start, end = map(int, parts)
            if end < start:
                raise ValueError(f"descending year range is not allowed: {token!r}")
            expanded = list(range(start, end + 1))
        elif token.isdigit():
            expanded = [int(token)]
        else:
            raise ValueError(f"invalid year token: {token!r}")
        for year in expanded:
            if not 1900 <= year <= 2100:
                raise ValueError(f"year outside supported range 1900-2100: {year}")
            if year in years:
                raise ValueError(f"duplicate year after expansion: {year}")
            years.append(year)
    return sorted(years)


def resolve_annual_drug_inputs(input_root: Path, years: list[int], filename: str) -> list[AnnualDrugInput]:
    if not years or len(years) != len(set(years)):
        raise ValueError("requested years must be nonempty and unique")
    resolved = []
    seen_paths = set()
    for year in sorted(years):
        directory = input_root / str(year)
        path = directory / filename
        if not directory.is_dir():
            raise FileNotFoundError(f"missing annual input directory: {directory}")
        if not path.is_file():
            raise FileNotFoundError(f"missing annual input parquet: {path}")
        real = path.resolve()
        if real in seen_paths:
            raise ValueError(f"annual input resolved more than once: {path}")
        seen_paths.add(real)
        try:
            schema = pq.ParquetFile(path).schema_arrow.names
        except Exception as exc:
            raise ValueError(f"cannot read annual input parquet {path}: {exc}") from exc
        missing = sorted(set(REQUIRED_COLUMNS) - set(schema))
        if missing:
            raise ValueError(f"annual input {path} missing required columns: {missing}")
        resolved.append(AnnualDrugInput(year, path))
    return resolved


def _load_assignment_maps(path: Path, years: list[int]) -> tuple[dict[int, dict[str, str]], int]:
    columns = ["caseid", "selected_primaryid", "assigned_report_year", "assignment_status"]
    assignments = pd.read_parquet(path, columns=columns)
    assignments["caseid"] = normalize_identifier(assignments["caseid"])
    assignments["selected_primaryid"] = normalize_identifier(assignments["selected_primaryid"])
    included = assignments[
        assignments["assignment_status"].eq("included")
        & assignments["assigned_report_year"].isin(years)
    ].copy()
    if included[["caseid", "selected_primaryid"]].isna().any().any():
        raise ValueError("included assignments contain blank identifiers")
    if included["caseid"].duplicated().any():
        raise ValueError("an included CASEID appears in more than one requested year")
    if included.groupby("selected_primaryid")["caseid"].nunique().gt(1).any():
        raise ValueError("a selected PRIMARYID maps to multiple CASEIDs")
    maps = {}
    for year in years:
        annual = included[included["assigned_report_year"].eq(year)]
        maps[year] = dict(zip(annual["selected_primaryid"], annual["caseid"]))
    return maps, len(included)


def _open_store(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.executescript("""
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        PRAGMA temp_store=FILE;
        CREATE TABLE stats (
            normalized TEXT NOT NULL, year INTEGER NOT NULL, row_count INTEGER NOT NULL,
            PRIMARY KEY (normalized, year)
        ) WITHOUT ROWID;
        CREATE TABLE norm_pid (
            normalized TEXT NOT NULL, primaryid TEXT NOT NULL,
            PRIMARY KEY (normalized, primaryid)
        ) WITHOUT ROWID;
        CREATE TABLE norm_raw (
            normalized TEXT NOT NULL, raw_name TEXT NOT NULL,
            PRIMARY KEY (normalized, raw_name)
        ) WITHOUT ROWID;
    """)
    return connection


def _upsert_batch(connection: sqlite3.Connection, frame: pd.DataFrame, year: int) -> None:
    if frame.empty:
        return
    counts = frame.groupby("normalized", sort=False).size()
    connection.executemany(
        "INSERT INTO stats VALUES (?, ?, ?) ON CONFLICT(normalized, year) "
        "DO UPDATE SET row_count=row_count+excluded.row_count",
        ((name, year, int(count)) for name, count in counts.items()),
    )
    connection.executemany(
        "INSERT OR IGNORE INTO norm_pid VALUES (?, ?)",
        frame[["normalized", "primaryid"]].drop_duplicates().itertuples(index=False, name=None),
    )
    connection.executemany(
        "INSERT OR IGNORE INTO norm_raw VALUES (?, ?)",
        frame[["normalized", "raw_name"]].drop_duplicates().itertuples(index=False, name=None),
    )


def _read_vocabulary(connection: sqlite3.Connection, sample_limit: int) -> pd.DataFrame:
    stats: dict[str, dict] = {}
    for name, year, count in connection.execute(
        "SELECT normalized, year, row_count FROM stats ORDER BY normalized, year"
    ):
        item = stats.setdefault(name, {"years": [], "rows": 0})
        item["years"].append(int(year)); item["rows"] += int(count)
    pid_counts = dict(connection.execute(
        "SELECT normalized, COUNT(*) FROM norm_pid GROUP BY normalized"
    ))
    raw_counts: Counter = Counter()
    raw_samples: dict[str, list[str]] = defaultdict(list)
    for name, raw in connection.execute("SELECT normalized, raw_name FROM norm_raw ORDER BY normalized, raw_name"):
        raw_counts[name] += 1
        if len(raw_samples[name]) < sample_limit:
            raw_samples[name].append(raw)
    rows = []
    for name in sorted(stats):
        years = stats[name]["years"]
        distinct_ids = int(pid_counts[name])
        rows.append({
            "normalized_drug_name": name, "first_seen_year": min(years),
            "last_seen_year": max(years), "years_seen": years, "year_count": len(years),
            "raw_row_count": stats[name]["rows"], "distinct_case_count": distinct_ids,
            "distinct_primaryid_count": distinct_ids, "distinct_raw_name_count": raw_counts[name],
            "sample_raw_names": raw_samples[name],
        })
    return pd.DataFrame(rows)


def build_multi_year_vocabulary(
    annual_inputs: list[AnnualDrugInput], *, assignments_path: Path = DEFAULT_ASSIGNMENTS,
    sample_raw_names_limit: int = 10, batch_size: int = 250_000,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    if sample_raw_names_limit < 1:
        raise ValueError("sample_raw_names_limit must be positive")
    inputs = sorted(annual_inputs, key=lambda item: item.analytical_year)
    years = [item.analytical_year for item in inputs]
    if not inputs or len(years) != len(set(years)):
        raise ValueError("annual inputs must be nonempty with unique years")
    assignment_maps, included_count = _load_assignment_maps(assignments_path, years)
    input_hashes_before = {item.path: sha256_file(item.path) for item in inputs}
    annual_rows = []
    total_blank_normalized = total_rows = 0
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="ms5_vocab_") as temp_dir:
        db_path = Path(temp_dir) / "vocabulary.sqlite"
        connection = _open_store(db_path)
        try:
            for item in inputs:
                year_started = time.perf_counter()
                expected = assignment_maps[item.analytical_year]
                seen_cases, seen_primaryids = set(), set()
                metrics = Counter()
                annual_raw_names, annual_norms, annual_combos = set(), set(), set()
                parquet = pq.ParquetFile(item.path)
                for record_batch in parquet.iter_batches(batch_size=batch_size, columns=list(REQUIRED_COLUMNS)):
                    frame = record_batch.to_pandas()
                    frame["caseid"] = normalize_identifier(frame["caseid"])
                    frame["primaryid"] = normalize_identifier(frame["primaryid"])
                    if frame[["caseid", "primaryid"]].isna().any().any():
                        raise ValueError(f"{item.path}: blank caseid or primaryid")
                    authoritative = frame["primaryid"].map(expected)
                    if authoritative.isna().any():
                        raise ValueError(f"{item.path}: unexpected selected PRIMARYID")
                    if authoritative.ne(frame["caseid"]).any():
                        raise ValueError(f"{item.path}: CASEID/selected PRIMARYID assignment mismatch")
                    seen_cases.update(frame["caseid"]); seen_primaryids.update(frame["primaryid"])
                    prod = frame["prod_ai"].fillna("").astype(str).str.strip()
                    drug = frame["drugname"].fillna("").astype(str).str.strip()
                    metrics["rows"] += len(frame)
                    metrics["blank_prod"] += int(prod.eq("").sum())
                    metrics["blank_drug"] += int(drug.eq("").sum())
                    metrics["both_blank"] += int((prod.eq("") & drug.eq("")).sum())
                    chosen = _choose_raw_drug_name(frame)
                    metrics["blank_chosen"] += int(chosen.eq("").sum())
                    normalized = chosen.map(normalize_drug_string)
                    # Do not promote common serialized-null placeholders to
                    # vocabulary keys; the frozen selector/normalizer outputs
                    # remain unchanged and this is an eligibility guard only.
                    blank_norm = normalized.eq("") | normalized.str.casefold().isin(INVALID_NORMALIZED_SENTINELS)
                    metrics["blank_norm"] += int(blank_norm.sum())
                    usable = pd.DataFrame({
                        "normalized": normalized[~blank_norm], "raw_name": chosen[~blank_norm],
                        "primaryid": frame.loc[~blank_norm, "primaryid"],
                    })
                    annual_raw_names.update(usable["raw_name"])
                    annual_norms.update(usable["normalized"])
                    annual_combos.update(name for name in usable["normalized"].unique() if "|" in name)
                    _upsert_batch(connection, usable, item.analytical_year)
                    connection.commit()
                if not seen_primaryids.issubset(expected):
                    raise AssertionError("resolved PRIMARYIDs are not a selected-assignment subset")
                annual_rows.append({
                    "year": item.analytical_year, "resolved_drug_rows": metrics["rows"],
                    "distinct_caseids": len(seen_cases), "distinct_primaryids": len(seen_primaryids),
                    "unique_selected_raw_names": len(annual_raw_names),
                    "unique_normalized_names": len(annual_norms), "blank_prod_ai_rows": metrics["blank_prod"],
                    "blank_drugname_rows": metrics["blank_drug"],
                    "both_raw_fields_blank_rows": metrics["both_blank"],
                    "blank_chosen_name_rows": metrics["blank_chosen"],
                    "blank_normalized_name_rows": metrics["blank_norm"],
                    "combination_normalized_names": len(annual_combos),
                    "runtime_seconds": round(time.perf_counter() - year_started, 3),
                })
                total_rows += metrics["rows"]; total_blank_normalized += metrics["blank_norm"]
            vocabulary = _read_vocabulary(connection, sample_raw_names_limit)
        finally:
            connection.close()
    input_hashes_after = {item.path: sha256_file(item.path) for item in inputs}
    if input_hashes_before != input_hashes_after:
        raise AssertionError("an annual resolved DRUG input changed during processing")
    if int(vocabulary["raw_row_count"].sum()) + total_blank_normalized != total_rows:
        raise AssertionError("vocabulary rows plus blank-normalized rows do not reconcile")
    if not vocabulary["normalized_drug_name"].is_unique:
        raise AssertionError("normalized vocabulary key is not unique")
    if not (vocabulary["raw_row_count"] >= vocabulary["distinct_primaryid_count"]).all():
        raise AssertionError("raw_row_count is smaller than distinct_primaryid_count")
    if not (vocabulary["distinct_primaryid_count"] == vocabulary["distinct_case_count"]).all():
        raise AssertionError("selected PRIMARYID/CASEID one-to-one contract failed")
    if not all(len(sample) <= sample_raw_names_limit for sample in vocabulary["sample_raw_names"]):
        raise AssertionError("raw-name sample limit exceeded")
    requested = set(years)
    if not all(row.year_count == len(row.years_seen) and set(row.years_seen) <= requested
               for row in vocabulary.itertuples(index=False)):
        raise AssertionError("years_seen invariant failed")
    for row in vocabulary.itertuples(index=False):
        if any(normalize_drug_string(value) != row.normalized_drug_name for value in row.sample_raw_names):
            raise AssertionError("a sample raw name does not normalize to its vocabulary key")
    metadata = {
        "requested_years": years, "input_row_count": total_rows,
        "included_assignments_in_scope": included_count,
        "blank_normalized_name_rows": total_blank_normalized,
        "input_hashes": {repository_path(path): digest for path, digest in input_hashes_before.items()},
        "total_runtime_seconds": round(time.perf_counter() - started, 3),
        "temporary_files_remaining": 0, "all_invariants_passed": True,
    }
    return vocabulary, pd.DataFrame(annual_rows), metadata


def _review_sample(vocabulary: pd.DataFrame) -> pd.DataFrame:
    selected: dict[str, set[str]] = defaultdict(set)
    for row in vocabulary.nlargest(50, "raw_row_count").itertuples(): selected[row.normalized_drug_name].add("top_frequent")
    collisions = vocabulary[vocabulary["distinct_raw_name_count"] > 1].nlargest(50, "raw_row_count")
    for row in collisions.itertuples(): selected[row.normalized_drug_name].add("high_frequency_collision")
    combos = vocabulary[vocabulary["normalized_drug_name"].str.contains("|", regex=False)].nlargest(50, "raw_row_count")
    for row in combos.itertuples(): selected[row.normalized_drug_name].add("combination")
    ordinary = vocabulary[~vocabulary["normalized_drug_name"].isin(selected)].copy()
    ordinary["_hash"] = ordinary["normalized_drug_name"].map(lambda value: hashlib.sha256(value.encode()).hexdigest())
    for row in ordinary.sort_values("_hash").head(50).itertuples():
        selected[row.normalized_drug_name].add("stable_hash")
    difference = vocabulary.assign(_difference=vocabulary.apply(
        lambda row: max((abs(len(raw) - len(row["normalized_drug_name"])) for raw in row["sample_raw_names"]), default=0), axis=1
    )).nlargest(50, "_difference")
    for row in difference.itertuples(): selected[row.normalized_drug_name].add("large_raw_normalized_difference")
    lookup = vocabulary.set_index("normalized_drug_name")
    rows = []
    for name in sorted(selected):
        row = lookup.loc[name]
        rows.append({"sample_category": "|".join(sorted(selected[name])),
                     "chosen_raw_name": row["sample_raw_names"][0] if len(row["sample_raw_names"]) else "",
                     "normalized_drug_name": name, "raw_row_count": row["raw_row_count"],
                     "distinct_case_count": row["distinct_case_count"],
                     "years_seen": json.dumps([
                         value.item() if hasattr(value, "item") else value
                         for value in row["years_seen"]
                     ]), "review_decision": "",
                     "review_comment": "", "expected_normalized_name": ""})
    return pd.DataFrame(rows)


def _git_metadata() -> tuple[str, bool]:
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True, check=True).stdout)
    return commit, dirty


def write_vocabulary_outputs(vocabulary: pd.DataFrame, annual_summary: pd.DataFrame, metadata: dict, *,
                             output: Path, report_dir: Path, normalization_version: str,
                             artifact_version: str, sample_raw_names_limit: int, overwrite: bool) -> list[Path]:
    paths = [output, report_dir / "ms5_drug_vocabulary_by_year.csv",
             report_dir / "ms5_drug_vocabulary_summary.json",
             report_dir / "ms5_normalization_collisions.csv",
             report_dir / "ms5_combination_names.csv",
             report_dir / "ms5_normalization_review_sample.csv"]
    existing = [path for path in paths if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"output exists; use --overwrite: {existing[0]}")
    output.parent.mkdir(parents=True, exist_ok=True); report_dir.mkdir(parents=True, exist_ok=True)
    staging = {path: path.with_name(path.name + ".part") for path in paths}
    for path in staging.values():
        if path.exists():
            path.unlink()
    collisions = vocabulary[vocabulary["distinct_raw_name_count"] > 1].copy()
    collisions = collisions.sort_values(["raw_row_count", "normalized_drug_name"], ascending=[False, True])
    combos = vocabulary[vocabulary["normalized_drug_name"].str.contains("|", regex=False)].copy()
    combos.insert(1, "component_count", combos["normalized_drug_name"].str.count(r"\|") + 1)
    combos = combos[["normalized_drug_name", "component_count", "raw_row_count", "distinct_case_count",
                    "first_seen_year", "last_seen_year", "years_seen", "sample_raw_names"]]
    commit, dirty = _git_metadata()
    summary = dict(metadata)
    summary.update({
        "study_start_year": min(metadata["requested_years"]), "study_end_year": max(metadata["requested_years"]),
        "year_count": len(metadata["requested_years"]), "input_files": list(metadata["input_hashes"]),
        "input_file_count": len(metadata["input_hashes"]),
        "distinct_case_count": int(annual_summary["distinct_caseids"].sum()),
        "distinct_primaryid_count": int(annual_summary["distinct_primaryids"].sum()),
        "unique_normalized_names": len(vocabulary),
        "blank_selected_name_rows": int(annual_summary["blank_chosen_name_rows"].sum()),
        "combination_name_count": len(combos), "normalization_collision_count": len(collisions),
        "normalization_version": normalization_version, "artifact_version": artifact_version,
        "sample_raw_names_limit": sample_raw_names_limit, "output_path": repository_path(output),
        "created_at": datetime.now(timezone.utc).isoformat(), "code_commit": commit,
        "worktree_dirty": dirty, "absolute_paths_persisted": False,
    })
    serializable = vocabulary.copy()
    # PyArrow supports deterministic list columns directly.
    try:
        serializable.to_parquet(staging[output], index=False)
        annual_summary.to_csv(staging[paths[1]], index=False)
        with staging[paths[2]].open("w", encoding="utf-8") as stream:
            json.dump(summary, stream, indent=2, sort_keys=True); stream.write("\n")
        collision_csv = collisions.copy()
        combo_csv = combos.copy()
        for frame in (collision_csv, combo_csv):
            for column in ("years_seen", "sample_raw_names"):
                frame[column] = frame[column].map(lambda values: json.dumps(
                    [value.item() if hasattr(value, "item") else value for value in values],
                    ensure_ascii=False,
                ))
        collision_csv.to_csv(staging[paths[3]], index=False)
        combo_csv.to_csv(staging[paths[4]], index=False)
        _review_sample(vocabulary).to_csv(staging[paths[5]], index=False)
        for path in paths:
            os.replace(staging[path], path)
    finally:
        for path in staging.values():
            if path.exists():
                path.unlink()
    return paths


def parse_args():
    parser = argparse.ArgumentParser(description="Build an MS5 multi-year vocabulary from globally corrected DRUG rows.")
    parser.add_argument("--years", nargs="+", required=True, help="Explicit years/ranges, e.g. 2018 2020 2023-2025")
    parser.add_argument("--annual-input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--input-filename", default="resolved_drug_rows.parquet")
    parser.add_argument("--global-case-assignments", type=Path, default=DEFAULT_ASSIGNMENTS)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--normalization-version", default="norm_v1")
    parser.add_argument("--artifact-version", default="vocab_v1")
    parser.add_argument("--sample-raw-names-limit", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=250_000)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args(); years = parse_years(args.years)
    scope = "_".join(map(str, years))
    output = args.output or DEFAULT_OUTPUT_DIR / f"multi_year_drug_vocabulary_{scope}_{args.normalization_version}.parquet"
    inputs = resolve_annual_drug_inputs(args.annual_input_root, years, args.input_filename)
    print(f"Requested years: {years}")
    for item in inputs: print(f"  {item.analytical_year}: {item.path}")
    vocabulary, annual, metadata = build_multi_year_vocabulary(
        inputs, assignments_path=args.global_case_assignments,
        sample_raw_names_limit=args.sample_raw_names_limit, batch_size=args.batch_size,
    )
    paths = write_vocabulary_outputs(
        vocabulary, annual, metadata, output=output, report_dir=args.report_dir,
        normalization_version=args.normalization_version, artifact_version=args.artifact_version,
        sample_raw_names_limit=args.sample_raw_names_limit, overwrite=args.overwrite,
    )
    print(f"Vocabulary rows: {len(vocabulary):,}")
    for path in paths: print(f"Saved: {path}")


if __name__ == "__main__":
    main()
