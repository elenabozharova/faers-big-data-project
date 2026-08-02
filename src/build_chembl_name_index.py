"""Build a versioned, ambiguity-preserving ChEMBL name index."""

from __future__ import annotations

import argparse
import json
import os
import re
import sqlite3
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import psycopg2
import pyarrow as pa
import pyarrow.parquet as pq

import config
from download_faers import sha256_file
from integrate_global_case_assignments import repository_path
from normalize_drug_names import normalize_drug_string


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/interim/drug_mapping/chembl36_name_index_norm_v1_index_v1.parquet"
DEFAULT_SUMMARY = ROOT / "reports/data_quality/ms5_chembl36_name_index_summary.json"
INVALID_NAMES = frozenset({"", "nan", "none", "null"})
SOURCE_TABLES = ("molecule_dictionary", "molecule_synonyms", "version")
INDEX_SCHEMA = pa.schema([
    ("normalized_chembl_name", pa.string()), ("chembl_id", pa.string()),
    ("source_name", pa.string()), ("name_type", pa.string()),
    ("source_table", pa.string()), ("chembl_version", pa.int16()),
    ("normalization_version", pa.string()),
])


def _git_metadata() -> tuple[str, bool]:
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True, check=True).stdout)
    return commit, dirty


def verify_chembl36(connection) -> dict:
    """Prove database/table identity before labelling an artifact ChEMBL 36."""
    with connection.cursor() as cursor:
        cursor.execute("SELECT current_database(), current_schema()")
        database_name, schema_name = cursor.fetchone()
        cursor.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema=current_schema() AND table_name = ANY(%s)",
            (list(SOURCE_TABLES),),
        )
        found = {row[0] for row in cursor.fetchall()}
        missing = set(SOURCE_TABLES) - found
        if missing:
            raise RuntimeError(f"ChEMBL schema missing required table(s): {sorted(missing)}")
        cursor.execute('SELECT * FROM "version"')
        columns = [description.name for description in cursor.description]
        rows = cursor.fetchall()
    evidence = " ".join(str(value) for row in rows for value in row if value is not None)
    database_matches = bool(re.fullmatch(r"chembl_?36", str(database_name), re.IGNORECASE))
    version_matches = bool(re.search(r"(?:chembl[_ -]?)?36(?:\D|$)", evidence, re.IGNORECASE))
    if not database_matches or not version_matches:
        raise RuntimeError(
            f"cannot verify ChEMBL 36: database={database_name!r}, "
            f"version table columns={columns}, values={rows[:5]!r}"
        )
    return {"database_name": database_name, "schema_name": schema_name,
            "version_table_columns": columns, "version_evidence": evidence[:500]}


def _stage_source_rows(connection, staging: sqlite3.Connection, *, batch_size: int) -> dict:
    queries = {
        "preferred": "SELECT chembl_id, pref_name FROM molecule_dictionary WHERE pref_name IS NOT NULL",
        "synonym": "SELECT md.chembl_id, ms.synonyms FROM molecule_synonyms ms "
                   "JOIN molecule_dictionary md ON md.molregno=ms.molregno WHERE ms.synonyms IS NOT NULL",
    }
    counts = {"preferred_source_rows": 0, "synonym_source_rows": 0,
              "blank_or_invalid_names_excluded": 0, "missing_chembl_ids_excluded": 0}
    for name_type, query in queries.items():
        cursor = connection.cursor(name=f"ms5_{name_type}")
        cursor.itersize = batch_size; cursor.execute(query)
        table = "molecule_dictionary" if name_type == "preferred" else "molecule_synonyms"
        key = f"{name_type}_source_rows"
        while True:
            rows = cursor.fetchmany(batch_size)
            if not rows: break
            counts[key] += len(rows)
            staged = []
            for chembl_id, source_name in rows:
                cid = "" if chembl_id is None else str(chembl_id).strip()
                if not cid:
                    counts["missing_chembl_ids_excluded"] += 1; continue
                original = "" if source_name is None else str(source_name).strip()
                normalized = normalize_drug_string(original)
                if normalized.casefold() in INVALID_NAMES:
                    counts["blank_or_invalid_names_excluded"] += 1; continue
                staged.append((normalized, cid, original, name_type, table))
            staging.executemany("INSERT OR IGNORE INTO names VALUES (?, ?, ?, ?, ?)", staged)
            staging.commit()
        cursor.close()
    return counts


def build_name_index(connection, output: Path, *, normalization_version: str,
                     batch_size: int = 100_000, overwrite: bool = False) -> dict:
    if output.exists() and not overwrite:
        raise FileExistsError(f"index exists; use --overwrite: {output}")
    verification = verify_chembl36(connection)
    output.parent.mkdir(parents=True, exist_ok=True)
    staged_output = output.with_name(output.name + ".part")
    if staged_output.exists(): staged_output.unlink()
    with tempfile.TemporaryDirectory(prefix="ms5_chembl_index_") as directory:
        staging = sqlite3.connect(Path(directory) / "names.sqlite")
        staging.execute("CREATE TABLE names (normalized TEXT, chembl_id TEXT, source_name TEXT, "
                        "name_type TEXT, source_table TEXT, PRIMARY KEY(normalized,name_type,chembl_id,source_name)) WITHOUT ROWID")
        try:
            counts = _stage_source_rows(connection, staging, batch_size=batch_size)
            index_rows = staging.execute("SELECT COUNT(*) FROM names").fetchone()[0]
            unique_names = staging.execute("SELECT COUNT(DISTINCT normalized) FROM names").fetchone()[0]
            unique_ids = staging.execute("SELECT COUNT(DISTINCT chembl_id) FROM names").fetchone()[0]
            ambiguous = staging.execute("SELECT COUNT(*) FROM (SELECT normalized FROM names GROUP BY normalized "
                                        "HAVING COUNT(DISTINCT chembl_id)>1)").fetchone()[0]
            writer = pq.ParquetWriter(staged_output, INDEX_SCHEMA, compression="snappy")
            try:
                cursor = staging.execute("SELECT normalized,chembl_id,source_name,name_type,source_table "
                                         "FROM names ORDER BY normalized,name_type,chembl_id,source_name")
                while True:
                    rows = cursor.fetchmany(batch_size)
                    if not rows: break
                    frame = pd.DataFrame(rows, columns=["normalized_chembl_name", "chembl_id", "source_name", "name_type", "source_table"])
                    frame["chembl_version"] = 36; frame["normalization_version"] = normalization_version
                    writer.write_table(pa.Table.from_pandas(frame, schema=INDEX_SCHEMA, preserve_index=False))
            finally:
                writer.close()
        except Exception:
            if staged_output.exists(): staged_output.unlink()
            raise
        finally:
            staging.close()
    os.replace(staged_output, output)
    return verification | counts | {"index_rows": index_rows, "unique_normalized_names": unique_names,
                                    "unique_chembl_ids": unique_ids, "ambiguous_normalized_names": ambiguous}


def write_summary(metrics: dict, output: Path, summary_path: Path, *, normalization_version: str,
                  index_artifact_version: str, overwrite: bool) -> dict:
    if summary_path.exists() and not overwrite:
        raise FileExistsError(f"summary exists; use --overwrite: {summary_path}")
    commit, dirty = _git_metadata()
    payload = metrics | {"chembl_version": 36, "normalization_version": normalization_version,
        "index_artifact_version": index_artifact_version, "source_tables": list(SOURCE_TABLES),
        "output_path": repository_path(output), "output_sha256": sha256_file(output),
        "created_at": datetime.now(timezone.utc).isoformat(), "code_commit": commit,
        "worktree_dirty": dirty, "all_invariants_passed": True}
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    staged = summary_path.with_name(summary_path.name + ".part")
    try:
        with staged.open("w", encoding="utf-8") as stream: json.dump(payload, stream, indent=2, sort_keys=True); stream.write("\n")
        os.replace(staged, summary_path)
    finally:
        if staged.exists(): staged.unlink()
    return payload


def parse_args():
    parser = argparse.ArgumentParser(description="Build a verified, versioned ChEMBL 36 preferred/synonym name index.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--normalization-version", default="norm_v1")
    parser.add_argument("--index-artifact-version", default="index_v1")
    parser.add_argument("--batch-size", type=int, default=100_000)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        connection = psycopg2.connect(dbname=config.CHEMBL_DB_NAME, user=config.CHEMBL_DB_USER,
            password=config.CHEMBL_DB_PASSWORD, host=config.CHEMBL_DB_HOST, port=config.CHEMBL_DB_PORT)
    except psycopg2.OperationalError as exc:
        raise SystemExit(
            f"Cannot connect to configured database {config.CHEMBL_DB_NAME!r}; "
            "ChEMBL 36 was not verified and no index was produced. "
            f"PostgreSQL error: {exc}"
        ) from exc
    try:
        metrics = build_name_index(connection, args.output, normalization_version=args.normalization_version,
                                   batch_size=args.batch_size, overwrite=args.overwrite)
    finally:
        connection.close()
    summary = write_summary(metrics, args.output, args.summary, normalization_version=args.normalization_version,
                            index_artifact_version=args.index_artifact_version, overwrite=args.overwrite)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__": main()
