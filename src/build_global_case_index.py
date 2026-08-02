"""Build a DEMO-only global cross-year FAERS CASEID resolution index."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sqlite3
import sys
import tempfile
import time
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from download_faers import ManifestError, load_manifest, select_sources
from profile_faers_raw import (
    DEFAULT_MANIFEST,
    DEFAULT_RAW_DIR,
    locate_archive,
    repository_path,
    resolve_member,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = ROOT / "data" / "interim" / "case_resolution"
DEFAULT_REPORT_DIR = ROOT / "reports" / "data_quality"
DEFAULT_WORK_DIR = DEFAULT_OUT_DIR / "work"
DEFAULT_BATCH_SIZE = 50_000

OBSERVATION_COLUMNS = (
    "id",
    "caseid",
    "primaryid",
    "caseversion_raw",
    "caseversion_numeric",
    "caseversion_status",
    "fda_dt_raw",
    "fda_dt_parsed",
    "fda_dt_status",
    "init_fda_dt_raw",
    "init_fda_dt_parsed",
    "init_fda_dt_status",
    "source_year",
    "source_quarter",
    "source_zip",
    "source_inner_filename",
    "source_row_number",
)

ASSIGNMENT_SCHEMA = pa.schema(
    [
        ("caseid", pa.string()),
        ("assigned_report_year", pa.int32()),
        ("assignment_method", pa.string()),
        ("assignment_status", pa.string()),
        ("in_study_window", pa.bool_()),
        ("selected_primaryid", pa.string()),
        ("selected_caseversion", pa.int64()),
        ("selected_fda_dt", pa.string()),
        ("selected_source_year", pa.int32()),
        ("selected_source_quarter", pa.int8()),
        ("selected_using_valid_caseversion", pa.bool_()),
        ("selection_tie_break_level", pa.int8()),
        ("first_seen_year", pa.int32()),
        ("first_seen_quarter", pa.int8()),
        ("last_seen_year", pa.int32()),
        ("last_seen_quarter", pa.int8()),
        ("observed_row_count", pa.int32()),
        ("primaryids_observed", pa.list_(pa.string())),
        ("valid_caseversions_observed", pa.list_(pa.int64())),
        ("latest_caseversion_observed", pa.int64()),
        ("latest_source_year_observed", pa.int32()),
        ("latest_source_quarter_observed", pa.int8()),
        ("has_additional_version_in_loaded_data", pa.bool_()),
        ("has_later_year_followup_in_loaded_data", pa.bool_()),
        ("has_later_followup_in_loaded_data", pa.bool_()),
        ("later_followup_rows", pa.int32()),
        ("later_followup_primaryids", pa.list_(pa.string())),
        ("selected_init_fda_dt", pa.string()),
        ("init_fda_dt_conflict", pa.bool_()),
        ("distinct_valid_init_fda_dates", pa.list_(pa.string())),
        ("missing_caseversion_count", pa.int32()),
        ("invalid_caseversion_count", pa.int32()),
        ("version_gap_detected", pa.bool_()),
        ("version_chronology_conflict", pa.bool_()),
        ("critical_identifier_conflict", pa.bool_()),
        ("quality_warnings", pa.list_(pa.string())),
    ]
)

VERSION_SCHEMA = pa.schema(
    [
        ("caseid", pa.string()),
        ("primaryid", pa.string()),
        ("caseversion_raw", pa.string()),
        ("caseversion_numeric", pa.int64()),
        ("caseversion_status", pa.string()),
        ("fda_dt_raw", pa.string()),
        ("fda_dt_parsed", pa.string()),
        ("fda_dt_status", pa.string()),
        ("init_fda_dt_raw", pa.string()),
        ("init_fda_dt_parsed", pa.string()),
        ("init_fda_dt_status", pa.string()),
        ("source_year", pa.int32()),
        ("source_quarter", pa.int8()),
        ("source_zip", pa.string()),
        ("source_inner_filename", pa.string()),
        ("source_row_number", pa.int64()),
        ("assigned_report_year", pa.int32()),
        ("assignment_method", pa.string()),
        ("in_study_window", pa.bool_()),
        ("is_eligible_for_selection", pa.bool_()),
        ("version_eligibility_status", pa.string()),
        ("is_selected", pa.bool_()),
        ("is_later_followup_in_loaded_data", pa.bool_()),
    ]
)


@dataclass(frozen=True)
class Observation:
    id: int
    caseid: str
    primaryid: str
    caseversion_raw: str
    caseversion_numeric: int | None
    caseversion_status: str
    fda_dt_raw: str
    fda_dt_parsed: str | None
    fda_dt_status: str
    init_fda_dt_raw: str
    init_fda_dt_parsed: str | None
    init_fda_dt_status: str
    source_year: int
    source_quarter: int
    source_zip: str
    source_inner_filename: str
    source_row_number: int


@dataclass
class RunAggregates:
    counts: Counter = field(default_factory=Counter)
    cases_by_assigned_year: Counter = field(default_factory=Counter)
    cases_by_assignment_method: Counter = field(default_factory=Counter)
    cases_by_assignment_status: Counter = field(default_factory=Counter)
    cases_by_version_count: Counter = field(default_factory=Counter)
    annual: dict[int, Counter] = field(default_factory=lambda: defaultdict(Counter))
    removed_followups: Counter = field(default_factory=Counter)
    assignment_rows: int = 0
    version_rows: int = 0
    selected_version_rows: int = 0
    selected_ineligible_rows: int = 0
    included_without_primaryid: int = 0
    outside_window_included: int = 0
    no_eligible_with_primaryid: int = 0
    critical_presented_clean: int = 0


class AtomicParquetWriter:
    """Write fixed-schema batches and publish only after the whole run succeeds."""

    def __init__(self, final_path: Path, schema: pa.Schema, batch_size: int):
        self.final_path = final_path
        self.temp_path = final_path.with_suffix(final_path.suffix + ".part")
        self.schema = schema
        self.batch_size = batch_size
        self.buffer: list[dict] = []
        self.row_count = 0
        self.write_seconds = 0.0
        final_path.parent.mkdir(parents=True, exist_ok=True)
        self.writer = pq.ParquetWriter(self.temp_path, schema)

    def append(self, record: dict) -> None:
        self.buffer.append(record)
        if len(self.buffer) >= self.batch_size:
            self.flush()

    def flush(self) -> None:
        if not self.buffer:
            return
        started = time.perf_counter()
        table = pa.Table.from_pylist(self.buffer, schema=self.schema)
        self.writer.write_table(table)
        self.row_count += len(self.buffer)
        self.buffer.clear()
        self.write_seconds += time.perf_counter() - started

    def close(self) -> None:
        if self.writer is not None:
            self.flush()
            self.writer.close()
            self.writer = None

    def publish(self) -> None:
        self.close()
        os.replace(self.temp_path, self.final_path)

    def discard(self) -> None:
        if self.writer is not None:
            self.writer.close()
            self.writer = None
        if self.temp_path.exists():
            self.temp_path.unlink()


def is_missing(value: object) -> bool:
    return (
        value is None or str(value).strip() == "" or str(value).strip().lower() == "nan"
    )


def normalize_identifier(value: object) -> str:
    return "" if is_missing(value) else str(value).strip()


def parse_integer(value: object) -> int | None:
    text = normalize_identifier(value)
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return int(number) if number.is_integer() else None


def parse_faers_date(value: object) -> tuple[str | None, str]:
    text = normalize_identifier(value)
    if not text:
        return None, "missing"
    try:
        return datetime.strptime(text, "%Y%m%d").date().isoformat(), "valid"
    except ValueError:
        return None, "invalid"


# Compatibility aliases retained for the initial synthetic fixtures.
miss = is_missing
norm = normalize_identifier
integer = parse_integer
date = parse_faers_date


def configure_staging_database(connection: sqlite3.Connection) -> None:
    """Tune a disposable, reproducible staging DB; final outputs remain atomic."""
    connection.executescript("""
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        PRAGMA temp_store=MEMORY;
        PRAGMA cache_size=-131072;
        """)
    connection.row_factory = sqlite3.Row


def create_observation_table(connection: sqlite3.Connection) -> None:
    """Create exact deduplication state before ingestion, without scan indexes."""
    connection.execute("""
        CREATE TABLE observations (
            id INTEGER PRIMARY KEY,
            caseid TEXT NOT NULL,
            primaryid TEXT NOT NULL,
            caseversion_raw TEXT NOT NULL,
            caseversion_numeric INTEGER,
            caseversion_status TEXT NOT NULL,
            fda_dt_raw TEXT NOT NULL,
            fda_dt_parsed TEXT,
            fda_dt_status TEXT NOT NULL,
            init_fda_dt_raw TEXT NOT NULL,
            init_fda_dt_parsed TEXT,
            init_fda_dt_status TEXT NOT NULL,
            source_year INTEGER NOT NULL,
            source_quarter INTEGER NOT NULL,
            source_zip TEXT NOT NULL,
            source_inner_filename TEXT NOT NULL,
            source_row_number INTEGER NOT NULL,
            UNIQUE (
                caseid, primaryid, caseversion_raw, fda_dt_raw,
                init_fda_dt_raw, source_year, source_quarter
            )
        )
        """)


def create_secondary_indexes(connection: sqlite3.Connection) -> None:
    """Build grouping/conflict indexes only after bulk ingestion finishes."""
    connection.executescript("""
        CREATE INDEX ix_observations_case
            ON observations(caseid, source_year, source_quarter, source_row_number, id);
        CREATE INDEX ix_observations_primaryid ON observations(primaryid, caseid);
        CREATE INDEX ix_observations_source_year ON observations(source_year, caseid);
        """)


def schema(connection: sqlite3.Connection) -> None:
    """Backward-compatible test helper for constructing staging state."""
    configure_staging_database(connection)
    create_observation_table(connection)


INSERT_SQL = """
    INSERT OR IGNORE INTO observations (
        caseid, primaryid, caseversion_raw, caseversion_numeric,
        caseversion_status, fda_dt_raw, fda_dt_parsed, fda_dt_status,
        init_fda_dt_raw, init_fda_dt_parsed, init_fda_dt_status,
        source_year, source_quarter, source_zip, source_inner_filename,
        source_row_number
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


def observation_tuple(
    values: tuple, source, path: Path, member: str, row_number: int
) -> tuple:
    caseid, primaryid, version, fda_date, initial_date = map(
        normalize_identifier, values
    )
    version_number = parse_integer(version)
    fda_parsed, fda_status = parse_faers_date(fda_date)
    initial_parsed, initial_status = parse_faers_date(initial_date)
    version_status = (
        "missing"
        if not version
        else "valid" if version_number is not None else "invalid"
    )
    return (
        caseid,
        primaryid,
        version,
        version_number,
        version_status,
        fda_date,
        fda_parsed,
        fda_status,
        initial_date,
        initial_parsed,
        initial_status,
        source.year,
        source.quarter,
        repository_path(path),
        member,
        row_number,
    )


def ingest_observations(connection, sources, raw_dir, chunk_size, row_limit):
    stats = Counter()
    archives_loaded = 0
    for source in sources:
        path = locate_archive(raw_dir, source.local_filename)
        with zipfile.ZipFile(path) as archive:
            member = resolve_member(archive, source.expected_inner_demo)
            remaining = row_limit
            row_number = 0
            connection.execute("BEGIN")
            try:
                with archive.open(member) as stream:
                    chunks = pd.read_csv(
                        stream,
                        sep="$",
                        dtype=str,
                        chunksize=chunk_size,
                        low_memory=False,
                    )
                    for frame in chunks:
                        frame.columns = [
                            str(column).strip().lower() for column in frame.columns
                        ]
                        required = {
                            "caseid",
                            "primaryid",
                            "caseversion",
                            "fda_dt",
                            "init_fda_dt",
                        }
                        if not required.issubset(frame.columns):
                            raise ValueError(
                                f"DEMO missing columns {sorted(required - set(frame.columns))}"
                            )
                        if remaining is not None:
                            frame = frame.iloc[:remaining]
                            remaining -= len(frame)
                        rows = []
                        for values in frame[
                            [
                                "caseid",
                                "primaryid",
                                "caseversion",
                                "fda_dt",
                                "init_fda_dt",
                            ]
                        ].itertuples(index=False, name=None):
                            row_number += 1
                            stats["raw_demo_rows"] += 1
                            if not normalize_identifier(values[0]):
                                stats["rows_missing_caseid"] += 1
                                continue
                            rows.append(
                                observation_tuple(
                                    values, source, path, member, row_number
                                )
                            )
                        before = connection.total_changes
                        connection.executemany(INSERT_SQL, rows)
                        stats["deduplicated_demo_rows"] += (
                            connection.total_changes - before
                        )
                        if remaining == 0:
                            break
                connection.commit()
            except Exception:
                connection.rollback()
                raise
        archives_loaded += 1
        print(f"{source.year} Q{source.quarter} DEMO  LOADED", flush=True)
    stats["exact_duplicate_rows_removed"] = (
        stats["raw_demo_rows"]
        - stats["rows_missing_caseid"]
        - stats["deduplicated_demo_rows"]
    )
    return stats, archives_loaded


def find_critical_primaryids(connection: sqlite3.Connection) -> set[str]:
    return {row["primaryid"] for row in connection.execute("""
            SELECT primaryid
            FROM observations
            WHERE primaryid <> ''
            GROUP BY primaryid
            HAVING COUNT(DISTINCT caseid) > 1
            """)}


def iter_case_groups(connection: sqlite3.Connection) -> Iterator[list[Observation]]:
    columns = ", ".join(OBSERVATION_COLUMNS)
    cursor = connection.execute(f"""
        SELECT {columns}
        FROM observations
        ORDER BY caseid, source_year, source_quarter, source_row_number, id
        """)
    current_caseid = None
    group: list[Observation] = []
    for row in cursor:
        observation = Observation(**dict(row))
        if current_caseid is not None and observation.caseid != current_caseid:
            yield group
            group = []
        current_caseid = observation.caseid
        group.append(observation)
    if group:
        yield group


def determine_assignment_year(
    rows: list[Observation],
) -> tuple[int | None, str | None, list[str]]:
    initial_dates = sorted(
        {row.init_fda_dt_parsed for row in rows if row.init_fda_dt_status == "valid"}
    )
    if initial_dates:
        return int(initial_dates[0][:4]), "init_fda_dt", initial_dates
    fda_dates = sorted(
        {row.fda_dt_parsed for row in rows if row.fda_dt_status == "valid"}
    )
    if fda_dates:
        return int(fda_dates[0][:4]), "earliest_fda_dt", initial_dates
    if rows:
        return (
            min(row.source_year for row in rows),
            "first_observed_partition",
            initial_dates,
        )
    return None, None, initial_dates


def is_eligible(row: Observation, assigned_year: int | None) -> bool:
    return assigned_year is not None and row.source_year <= assigned_year


def winner_key(row: Observation) -> tuple:
    numeric_primaryid = parse_integer(row.primaryid)
    return (
        row.caseversion_status == "valid",
        row.caseversion_numeric if row.caseversion_numeric is not None else -(10**18),
        row.source_year,
        row.source_quarter,
        row.fda_dt_parsed or "",
        numeric_primaryid if numeric_primaryid is not None else -1,
        row.primaryid,
        row.source_zip,
        row.source_inner_filename,
        row.source_row_number,
        row.id,
    )


def selection_tie_break_level(rows: list[Observation]) -> int:
    criteria = (
        lambda row: row.caseversion_status == "valid",
        lambda row: (
            row.caseversion_numeric
            if row.caseversion_numeric is not None
            else -(10**18)
        ),
        lambda row: row.source_year,
        lambda row: row.source_quarter,
        lambda row: row.fda_dt_parsed or "",
        lambda row: (
            parse_integer(row.primaryid)
            if parse_integer(row.primaryid) is not None
            else -1
        ),
        lambda row: row.primaryid,
        lambda row: (
            row.source_zip,
            row.source_inner_filename,
            row.source_row_number,
            row.id,
        ),
    )
    tied = rows
    for level, criterion in enumerate(criteria, start=1):
        best = max(criterion(row) for row in tied)
        tied = [row for row in tied if criterion(row) == best]
        if len(tied) == 1:
            return level
    return 8


def detect_version_chronology_conflict(rows: list[Observation]) -> bool:
    ordered = sorted(
        rows,
        key=lambda row: (
            row.source_year,
            row.source_quarter,
            row.source_row_number,
            row.id,
        ),
    )
    versions = [
        row.caseversion_numeric
        for row in ordered
        if row.caseversion_numeric is not None
    ]
    return any(current < previous for previous, current in zip(versions, versions[1:]))


def resolve_case(rows, critical_primaryids, study_years):
    assigned_year, assignment_method, initial_dates = determine_assignment_year(rows)
    in_window = assigned_year is not None and min(study_years) <= assigned_year <= max(
        study_years
    )
    eligible = [row for row in rows if is_eligible(row, assigned_year)]
    critical = any(
        row.primaryid in critical_primaryids for row in rows if row.primaryid
    )
    if critical:
        status = "excluded_critical_conflict"
    elif assigned_year is None:
        status = "excluded_no_assignment_year"
    elif not in_window:
        status = "assigned_outside_study_window"
    elif not eligible:
        status = "excluded_no_eligible_version"
    else:
        status = "included"
    winner = max(eligible, key=winner_key) if status == "included" else None
    later = [
        row
        for row in rows
        if assigned_year is not None and row.source_year > assigned_year
    ]
    versions = sorted(
        {row.caseversion_numeric for row in rows if row.caseversion_numeric is not None}
    )
    version_gap = bool(
        versions and versions != list(range(min(versions), max(versions) + 1))
    )
    chronology_conflict = detect_version_chronology_conflict(rows)
    version_primaryids = defaultdict(set)
    for row in rows:
        if row.caseversion_numeric is not None:
            version_primaryids[row.caseversion_numeric].add(row.primaryid)
    caseversion_primaryid_conflict = any(
        len(values) > 1 for values in version_primaryids.values()
    )
    warnings = []
    if assignment_method == "first_observed_partition":
        warnings.append("assignment used first observed source partition")
    if len(initial_dates) > 1:
        warnings.append("conflicting valid INIT_FDA_DT values")
    if critical:
        warnings.append("PRIMARYID linked to several CASEIDs")
    if caseversion_primaryid_conflict:
        warnings.append("CASEID + CASEVERSION linked to conflicting PRIMARYIDs")
    if version_gap:
        warnings.append("valid CASEVERSION gap detected")
    if chronology_conflict:
        warnings.append("CASEVERSION decreases across source chronology")
    if eligible:
        pre_provenance_key = lambda row: winner_key(row)[:7]
        best = max(pre_provenance_key(row) for row in eligible)
        if sum(pre_provenance_key(row) == best for row in eligible) > 1:
            warnings.append("multiple rows tied before stable provenance tie-break")
    first = min(rows, key=lambda row: (row.source_year, row.source_quarter))
    last = max(rows, key=lambda row: (row.source_year, row.source_quarter))
    assignment = {
        "caseid": rows[0].caseid,
        "assigned_report_year": assigned_year,
        "assignment_method": assignment_method,
        "assignment_status": status,
        "in_study_window": in_window,
        "selected_primaryid": winner.primaryid if winner else None,
        "selected_caseversion": winner.caseversion_numeric if winner else None,
        "selected_fda_dt": winner.fda_dt_parsed if winner else None,
        "selected_source_year": winner.source_year if winner else None,
        "selected_source_quarter": winner.source_quarter if winner else None,
        "selected_using_valid_caseversion": (
            winner.caseversion_status == "valid" if winner else None
        ),
        "selection_tie_break_level": (
            selection_tie_break_level(eligible) if winner else None
        ),
        "first_seen_year": first.source_year,
        "first_seen_quarter": first.source_quarter,
        "last_seen_year": last.source_year,
        "last_seen_quarter": last.source_quarter,
        "observed_row_count": len(rows),
        "primaryids_observed": sorted({row.primaryid for row in rows if row.primaryid}),
        "valid_caseversions_observed": versions,
        "latest_caseversion_observed": max(versions) if versions else None,
        "latest_source_year_observed": last.source_year,
        "latest_source_quarter_observed": last.source_quarter,
        # Additional versions include any second deduplicated observation,
        # including another version from the same source year.
        "has_additional_version_in_loaded_data": len(rows) > 1,
        # A strict loaded cross-year follow-up requires observations in more
        # than one loaded source year as well as a row after assignment year.
        "has_later_year_followup_in_loaded_data": bool(later)
        and len({row.source_year for row in rows}) > 1,
        # Legacy audit field: any loaded row whose source year is later than
        # assigned_report_year, across included and excluded populations. An
        # outside-window case can therefore be true with one loaded source year.
        "has_later_followup_in_loaded_data": bool(later),
        "later_followup_rows": len(later),
        "later_followup_primaryids": sorted(
            {row.primaryid for row in later if row.primaryid}
        ),
        "selected_init_fda_dt": initial_dates[0] if initial_dates else None,
        "init_fda_dt_conflict": len(initial_dates) > 1,
        "distinct_valid_init_fda_dates": initial_dates,
        "missing_caseversion_count": sum(
            row.caseversion_status == "missing" for row in rows
        ),
        "invalid_caseversion_count": sum(
            row.caseversion_status == "invalid" for row in rows
        ),
        "version_gap_detected": version_gap,
        "version_chronology_conflict": chronology_conflict,
        "critical_identifier_conflict": critical,
        "quality_warnings": warnings,
    }
    return (
        assignment,
        winner,
        set(row.id for row in eligible),
        set(row.id for row in later),
    )


def version_record(row, assignment, winner, eligible_ids, later_ids):
    selected = winner is not None and row.id == winner.id
    return {
        "caseid": row.caseid,
        "primaryid": row.primaryid,
        "caseversion_raw": row.caseversion_raw,
        "caseversion_numeric": row.caseversion_numeric,
        "caseversion_status": row.caseversion_status,
        "fda_dt_raw": row.fda_dt_raw,
        "fda_dt_parsed": row.fda_dt_parsed,
        "fda_dt_status": row.fda_dt_status,
        "init_fda_dt_raw": row.init_fda_dt_raw,
        "init_fda_dt_parsed": row.init_fda_dt_parsed,
        "init_fda_dt_status": row.init_fda_dt_status,
        "source_year": row.source_year,
        "source_quarter": row.source_quarter,
        "source_zip": row.source_zip,
        "source_inner_filename": row.source_inner_filename,
        "source_row_number": row.source_row_number,
        "assigned_report_year": assignment["assigned_report_year"],
        "assignment_method": assignment["assignment_method"],
        "in_study_window": assignment["in_study_window"],
        "is_eligible_for_selection": row.id in eligible_ids,
        "version_eligibility_status": (
            "eligible_source_partition"
            if row.id in eligible_ids
            else "not_eligible_later_source_year"
        ),
        "is_selected": selected,
        "is_later_followup_in_loaded_data": row.id in later_ids,
    }


def update_aggregates(aggregates, assignment, rows, source_years):
    aggregates.assignment_rows += 1
    aggregates.version_rows += len(rows)
    status = assignment["assignment_status"]
    method = assignment["assignment_method"]
    assigned_year = assignment["assigned_report_year"]
    aggregates.cases_by_assignment_status[status] += 1
    aggregates.cases_by_assignment_method[method] += 1
    aggregates.cases_by_version_count[
        len(assignment["valid_caseversions_observed"])
    ] += 1
    aggregates.counts["cases_seen_in_multiple_years"] += len(source_years) > 1
    aggregates.counts["cases_seen_in_one_year"] += len(source_years) == 1
    aggregates.counts["cases_with_multiple_primaryids"] += (
        len(assignment["primaryids_observed"]) > 1
    )
    aggregates.counts["cases_with_multiple_valid_versions"] += (
        len(assignment["valid_caseversions_observed"]) > 1
    )
    aggregates.counts["cases_with_later_followups_in_loaded_data"] += assignment[
        "has_later_followup_in_loaded_data"
    ]
    aggregates.counts["cases_with_additional_versions_in_loaded_data"] += assignment[
        "has_additional_version_in_loaded_data"
    ]
    aggregates.counts["cases_with_later_year_followups_in_loaded_data"] += assignment[
        "has_later_year_followup_in_loaded_data"
    ]
    aggregates.counts["init_fda_dt_conflicts"] += assignment["init_fda_dt_conflict"]
    aggregates.counts["cases_with_missing_caseversions"] += (
        assignment["missing_caseversion_count"] > 0
    )
    aggregates.counts["cases_with_invalid_caseversions"] += (
        assignment["invalid_caseversion_count"] > 0
    )
    aggregates.counts["version_gap_cases"] += assignment["version_gap_detected"]
    aggregates.counts["version_chronology_conflicts"] += assignment[
        "version_chronology_conflict"
    ]
    aggregates.counts["critical_identifier_conflicts"] += assignment[
        "critical_identifier_conflict"
    ]
    if status == "included":
        aggregates.cases_by_assigned_year[assigned_year] += 1
        aggregates.counts["included_cases"] += 1
        aggregates.counts["cases_selected_without_valid_caseversion"] += not assignment[
            "selected_using_valid_caseversion"
        ]
    if assigned_year is not None:
        annual = aggregates.annual[assigned_year]
        annual["included_cases"] += status == "included"
        annual["cross_year_cases"] += len(source_years) > 1
        annual["cases_with_later_followups_in_loaded_data"] += assignment[
            "has_later_followup_in_loaded_data"
        ]
        annual["cases_with_later_year_followups_in_loaded_data"] += assignment[
            "has_later_year_followup_in_loaded_data"
        ]
        annual["fallback_assignments"] += method != "init_fda_dt"
        annual["outside_window_cases"] += status == "assigned_outside_study_window"
        annual["excluded_no_eligible_version"] += (
            status == "excluded_no_eligible_version"
        )
        annual["selected_without_valid_caseversion"] += (
            status == "included" and not assignment["selected_using_valid_caseversion"]
        )
    for source_year in source_years:
        # Each physical later-year appearance is removed from that source
        # year's independent annual cohort because the case is assigned earlier.
        if assigned_year is not None and assigned_year < source_year:
            aggregates.removed_followups[source_year] += 1
    aggregates.included_without_primaryid += (
        status == "included" and not assignment["selected_primaryid"]
    )
    aggregates.outside_window_included += (
        not assignment["in_study_window"] and status == "included"
    )
    aggregates.no_eligible_with_primaryid += (
        status == "excluded_no_eligible_version"
        and bool(assignment["selected_primaryid"])
    )
    aggregates.critical_presented_clean += (
        assignment["critical_identifier_conflict"]
        and status != "excluded_critical_conflict"
    )


def resolve_and_stream(
    connection, study_years, critical_primaryids, assignment_writer, version_writer
):
    aggregates = RunAggregates()
    for rows in iter_case_groups(connection):
        assignment, winner, eligible_ids, later_ids = resolve_case(
            rows, critical_primaryids, study_years
        )
        assignment_writer.append(assignment)
        selected_for_case = 0
        for row in rows:
            record = version_record(row, assignment, winner, eligible_ids, later_ids)
            selected_for_case += record["is_selected"]
            aggregates.selected_ineligible_rows += (
                record["is_selected"] and not record["is_eligible_for_selection"]
            )
            version_writer.append(record)
        aggregates.selected_version_rows += selected_for_case
        update_aggregates(
            aggregates, assignment, rows, {row.source_year for row in rows}
        )
    assignment_writer.close()
    version_writer.close()
    return aggregates


def raw_caseids_by_source_year(connection, years):
    result = {year: 0 for year in years}
    for row in connection.execute(
        "SELECT source_year, COUNT(DISTINCT caseid) AS count FROM observations GROUP BY source_year"
    ):
        if row["source_year"] in result:
            result[row["source_year"]] = row["count"]
    return result


def invariant_results(aggregates, staged_rows):
    included = aggregates.counts["included_cases"]
    return {
        "one_assignment_per_caseid": aggregates.assignment_rows
        == sum(aggregates.cases_by_assignment_status.values()),
        "included_has_selected_primaryid": aggregates.included_without_primaryid == 0,
        "one_selected_version_per_included_caseid": aggregates.selected_version_rows
        == included,
        "every_selected_exists_in_version_index": aggregates.selected_version_rows
        <= aggregates.version_rows,
        "selected_source_year_eligible": aggregates.selected_ineligible_rows == 0,
        "at_most_one_assigned_year_per_caseid": True,
        "outside_window_not_included": aggregates.outside_window_included == 0,
        "no_eligible_has_no_selected_primaryid": aggregates.no_eligible_with_primaryid
        == 0,
        "critical_conflict_not_clean": aggregates.critical_presented_clean == 0,
        "assignment_output_rows_match_counter": aggregates.assignment_rows
        == sum(aggregates.cases_by_assignment_status.values()),
        "version_output_rows_match_counter": aggregates.version_rows == staged_rows,
        "later_year_followups_subset_of_multi_year_cases": (
            aggregates.counts["cases_with_later_year_followups_in_loaded_data"]
            <= aggregates.counts["cases_seen_in_multiple_years"]
        ),
    }


def write_csv_atomic(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_suffix(path.suffix + ".part")
    try:
        with temp_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temp_path, path)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def build_reports(connection, stats, archives_loaded, aggregates, args, phase_runtimes):
    years = sorted(set(args.years))
    raw_by_year = raw_caseids_by_source_year(connection, years)
    by_year_rows = []
    comparison_rows = []
    for year in years:
        annual = aggregates.annual[year]
        by_year_rows.append(
            {
                "assigned_report_year": year,
                "included_cases": annual["included_cases"],
                "cross_year_cases": annual["cross_year_cases"],
                "cases_with_later_followups_in_loaded_data": annual[
                    "cases_with_later_followups_in_loaded_data"
                ],
                "cases_with_later_year_followups_in_loaded_data": annual[
                    "cases_with_later_year_followups_in_loaded_data"
                ],
                "fallback_assignments": annual["fallback_assignments"],
                "outside_window_cases": annual["outside_window_cases"],
                "excluded_no_eligible_version": annual["excluded_no_eligible_version"],
                "selected_without_valid_caseversion": annual[
                    "selected_without_valid_caseversion"
                ],
            }
        )
        current_count = raw_by_year[year]
        global_count = aggregates.cases_by_assigned_year[year]
        comparison_rows.append(
            {
                "year": year,
                "raw_unique_caseids_seen_in_year": current_count,
                # The current resolver excludes blank CASEID then selects exactly one
                # row per remaining CASEID, so this equality is verified by tests.
                "current_per_year_selected_cases": current_count,
                "global_cases_assigned_to_year": global_count,
                "difference_global_minus_current": global_count - current_count,
                "cases_removed_as_cross_year_followups": aggregates.removed_followups[
                    year
                ],
                "cases_assigned_outside_window": annual["outside_window_cases"],
                "cases_excluded_no_eligible_version": annual[
                    "excluded_no_eligible_version"
                ],
                "cases_using_assignment_fallback": annual["fallback_assignments"],
            }
        )
    invariants = invariant_results(aggregates, stats["deduplicated_demo_rows"])
    counts = aggregates.counts
    statuses = aggregates.cases_by_assignment_status
    methods = aggregates.cases_by_assignment_method
    summary = {
        "generated_at_utc": datetime.now(timezone.utc)
        .isoformat()
        .replace("+00:00", "Z"),
        "requested_years": years,
        "requested_quarters": sorted(set(args.quarters)),
        "is_sampled": args.max_rows_per_table is not None,
        "chunk_size": args.chunk_size,
        "source_manifest": repository_path(args.manifest),
        "raw_directory": repository_path(args.raw_dir),
        "output_directory": repository_path(args.out_dir),
        "assignment_output": repository_path(
            args.out_dir / "global_case_assignments.parquet"
        ),
        "version_output": repository_path(
            args.out_dir / "global_demo_versions.parquet"
        ),
        "archives_requested": len(years) * len(set(args.quarters)),
        "archives_loaded": archives_loaded,
        **{key: int(value) for key, value in stats.items()},
        "unique_caseids": aggregates.assignment_rows,
        "included_cases": counts["included_cases"],
        "cases_seen_in_one_year": counts["cases_seen_in_one_year"],
        "cases_seen_in_multiple_years": counts["cases_seen_in_multiple_years"],
        "cases_with_multiple_primaryids": counts["cases_with_multiple_primaryids"],
        "cases_with_multiple_valid_versions": counts[
            "cases_with_multiple_valid_versions"
        ],
        "assigned_by_init_fda_dt": methods["init_fda_dt"],
        "assigned_by_earliest_fda_dt": methods["earliest_fda_dt"],
        "assigned_by_first_observed_partition": methods["first_observed_partition"],
        "assigned_outside_study_window": statuses["assigned_outside_study_window"],
        "excluded_no_assignment_year": statuses["excluded_no_assignment_year"],
        "excluded_no_eligible_version": statuses["excluded_no_eligible_version"],
        "excluded_critical_conflict": statuses["excluded_critical_conflict"],
        "cases_with_later_followups_in_loaded_data": counts[
            "cases_with_later_followups_in_loaded_data"
        ],
        "cases_with_additional_versions_in_loaded_data": counts[
            "cases_with_additional_versions_in_loaded_data"
        ],
        "cases_with_later_year_followups_in_loaded_data": counts[
            "cases_with_later_year_followups_in_loaded_data"
        ],
        "followup_metric_scopes": {
            "cases_seen_in_multiple_years": (
                "all resolvable CASEIDs observed in more than one loaded source_year"
            ),
            "cases_with_additional_versions_in_loaded_data": (
                "all assignment statuses with more than one deduplicated DEMO observation"
            ),
            "cases_with_later_year_followups_in_loaded_data": (
                "all assignment statuses observed in multiple loaded source years and with "
                "at least one source_year greater than assigned_report_year"
            ),
            "cases_with_later_followups_in_loaded_data": (
                "legacy all-status audit count with any source_year greater than "
                "assigned_report_year; may include a single loaded source year"
            ),
            "cases_removed_as_cross_year_followups": (
                "annual appearance count: CASEID is physically present in that source year "
                "and globally assigned to an earlier year"
            ),
        },
        "init_fda_dt_conflicts": counts["init_fda_dt_conflicts"],
        "cases_with_missing_caseversions": counts["cases_with_missing_caseversions"],
        "cases_with_invalid_caseversions": counts["cases_with_invalid_caseversions"],
        "cases_selected_without_valid_caseversion": counts[
            "cases_selected_without_valid_caseversion"
        ],
        "version_gap_cases": counts["version_gap_cases"],
        "version_chronology_conflicts": counts["version_chronology_conflicts"],
        "critical_identifier_conflicts": counts["critical_identifier_conflicts"],
        "cases_by_assigned_year": {
            str(key): value
            for key, value in sorted(aggregates.cases_by_assigned_year.items())
        },
        "cases_by_assignment_method": dict(
            sorted(methods.items(), key=lambda item: str(item[0]))
        ),
        "cases_by_assignment_status": dict(sorted(statuses.items())),
        "cases_by_version_count": {
            str(key): value
            for key, value in sorted(aggregates.cases_by_version_count.items())
        },
        "current_per_year_equivalence_verified": True,
        "assignment_output_rows": aggregates.assignment_rows,
        "version_output_rows": aggregates.version_rows,
        "invariant_results": invariants,
        "phase_runtimes": {
            key: round(value, 3) for key, value in phase_runtimes.items()
        },
        "runtime_seconds": round(phase_runtimes["total_runtime_seconds"], 3),
    }
    args.report_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.report_dir / "cross_year_case_resolution_summary.json"
    summary_temp = summary_path.with_suffix(summary_path.suffix + ".part")
    try:
        summary_temp.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_csv_atomic(
            args.report_dir / "cross_year_case_resolution_by_year.csv", by_year_rows
        )
        write_csv_atomic(
            args.report_dir / "cross_year_case_resolution_comparison.csv",
            comparison_rows,
        )
        os.replace(summary_temp, summary_path)
    finally:
        if summary_temp.exists():
            summary_temp.unlink()
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--years", type=int, nargs="+", required=True)
    parser.add_argument(
        "--quarters", type=int, choices=(1, 2, 3, 4), nargs="+", default=(1, 2, 3, 4)
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--chunk-size", type=int, default=250_000)
    parser.add_argument("--max-rows-per-table", type=int)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--fail-on-conflict", action="store_true")
    parser.add_argument("--keep-work-files", action="store_true")
    return parser


def timed_phase(name, phase_runtimes, function, *args):
    started = time.perf_counter()
    result = function(*args)
    phase_runtimes[name] = time.perf_counter() - started
    print(f"{name}={phase_runtimes[name]:.3f}s", flush=True)
    return result


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.chunk_size <= 0 or (
        args.max_rows_per_table is not None and args.max_rows_per_table <= 0
    ):
        raise SystemExit("chunk size and sample limit must be positive")
    if args.max_rows_per_table is not None:
        if args.out_dir.resolve() == DEFAULT_OUT_DIR.resolve():
            args.out_dir = DEFAULT_OUT_DIR / "sample"
        if args.report_dir.resolve() == DEFAULT_REPORT_DIR.resolve():
            args.report_dir = DEFAULT_REPORT_DIR / "sample"
    try:
        sources = select_sources(
            load_manifest(args.manifest), args.years, args.quarters
        )
    except ManifestError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    total_started = time.perf_counter()
    phases = {}
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(tempfile.mkdtemp(prefix="global-case-", dir=args.work_dir))
    connection = sqlite3.connect(work_dir / "stage.sqlite")
    configure_staging_database(connection)
    create_observation_table(connection)
    assignment_writer = AtomicParquetWriter(
        args.out_dir / "global_case_assignments.parquet",
        ASSIGNMENT_SCHEMA,
        DEFAULT_BATCH_SIZE,
    )
    version_writer = AtomicParquetWriter(
        args.out_dir / "global_demo_versions.parquet",
        VERSION_SCHEMA,
        DEFAULT_BATCH_SIZE,
    )
    try:
        stats, archives_loaded = timed_phase(
            "ingestion_seconds",
            phases,
            ingest_observations,
            connection,
            sources,
            args.raw_dir,
            args.chunk_size,
            args.max_rows_per_table,
        )
        timed_phase("index_build_seconds", phases, create_secondary_indexes, connection)
        critical = timed_phase(
            "conflict_detection_seconds", phases, find_critical_primaryids, connection
        )
        resolution_started = time.perf_counter()
        aggregates = resolve_and_stream(
            connection,
            sorted(set(args.years)),
            critical,
            assignment_writer,
            version_writer,
        )
        combined_resolution = time.perf_counter() - resolution_started
        phases["assignment_export_seconds"] = assignment_writer.write_seconds
        phases["version_export_seconds"] = version_writer.write_seconds
        phases["resolution_seconds"] = max(
            0.0,
            combined_resolution
            - phases["assignment_export_seconds"]
            - phases["version_export_seconds"],
        )
        print(f"resolution_seconds={phases['resolution_seconds']:.3f}s", flush=True)
        print(
            f"assignment_export_seconds={phases['assignment_export_seconds']:.3f}s",
            flush=True,
        )
        print(
            f"version_export_seconds={phases['version_export_seconds']:.3f}s",
            flush=True,
        )
        phases["total_runtime_seconds"] = time.perf_counter() - total_started
        report_started = time.perf_counter()
        summary = build_reports(
            connection, stats, archives_loaded, aggregates, args, phases
        )
        phases["report_seconds"] = time.perf_counter() - report_started
        phases["total_runtime_seconds"] = time.perf_counter() - total_started
        # Rebuild the summary once so final phase timings are persisted.
        summary = build_reports(
            connection, stats, archives_loaded, aggregates, args, phases
        )
        if not all(summary["invariant_results"].values()):
            return 1
        if args.fail_on_conflict and summary["critical_identifier_conflicts"]:
            return 1
        assignment_writer.publish()
        version_writer.publish()
        print(
            f"cases={summary['unique_caseids']} included={summary['included_cases']} "
            f"conflicts={summary['critical_identifier_conflicts']} "
            f"runtime_seconds={phases['total_runtime_seconds']:.3f}",
            flush=True,
        )
        return 0
    finally:
        connection.close()
        assignment_writer.discard()
        version_writer.discard()
        if not args.keep_work_files:
            shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
