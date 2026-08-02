"""Read-only audit of globally selected cases missing resolved REAC rows."""

from __future__ import annotations

import argparse
import json
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

from download_faers import load_manifest, select_sources, sha256_file
from integrate_global_case_assignments import (
    DEFAULT_ASSIGNMENTS,
    DEFAULT_REPORT_DIR,
    load_global_assignments,
    normalize_identifier,
    repository_path,
    selected_cases_for_year,
)
from profile_faers_raw import DEFAULT_MANIFEST, DEFAULT_RAW_DIR, locate_archive, resolve_member


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_VERSIONS = ROOT / "data/interim/case_resolution/global_demo_versions.parquet"
DEFAULT_RESOLVED = ROOT / "data/processed/faers_case_level_annual/2024/resolved_reac_rows.parquet"
VERSION_COLUMNS = [
    "caseid", "primaryid", "caseversion_numeric", "init_fda_dt_raw", "fda_dt_raw",
    "source_year", "source_quarter", "is_selected", "assigned_report_year", "assignment_method",
]


def identify_missing_selected_reac(
    assignments: pd.DataFrame, resolved_reac: pd.DataFrame, year: int,
) -> pd.DataFrame:
    selected = selected_cases_for_year(assignments, year)
    resolved_ids = set(normalize_identifier(resolved_reac["primaryid"]).dropna())
    return selected[~selected["selected_primaryid"].isin(resolved_ids)].copy().reset_index(drop=True)


def validate_selected_provenance(missing: pd.DataFrame, versions: pd.DataFrame) -> pd.DataFrame:
    work = versions.copy()
    work["caseid"] = normalize_identifier(work["caseid"])
    work["primaryid"] = normalize_identifier(work["primaryid"])
    selected_versions = work[work["is_selected"].eq(True)]
    joined = missing.merge(
        selected_versions,
        left_on=["caseid", "selected_primaryid"], right_on=["caseid", "primaryid"],
        how="left", suffixes=("_assignment", "_version"), validate="one_to_one",
    )
    if joined["primaryid"].isna().any():
        raise ValueError("a missing selected assignment has no selected MS4A version record")
    for field in ("source_year", "source_quarter"):
        if not joined[f"selected_{field}"].astype(int).eq(joined[field].astype(int)).all():
            raise ValueError(f"selected {field} disagrees between assignment and version provenance")
    return joined


def scan_reac_sources(sources, raw_dir: Path, target_primaryids: set[str], target_caseids: set[str], chunk_size: int):
    """Scan each manifest REAC member once using production parser semantics."""
    occurrences = defaultdict(lambda: {"rows": 0, "caseids": set(), "pts": set()})
    case_primary_counts = defaultdict(int)
    chunks = rows = 0
    for source in sources:
        archive_path = locate_archive(raw_dir, source.local_filename)
        with zipfile.ZipFile(archive_path) as archive:
            member = resolve_member(archive, source.expected_inner_reac)
            with archive.open(member) as stream:
                iterator = pd.read_csv(
                    stream, sep="$", dtype=str, chunksize=chunk_size, low_memory=False,
                    usecols=lambda name: str(name).strip().lower() in {"primaryid", "caseid", "pt"},
                )
                for chunk in iterator:
                    chunks += 1
                    rows += len(chunk)
                    chunk.columns = chunk.columns.astype(str).str.strip().str.lower()
                    if "primaryid" not in chunk:
                        raise ValueError(f"REAC member {member!r} has no primaryid")
                    chunk["primaryid"] = normalize_identifier(chunk["primaryid"])
                    if "caseid" in chunk:
                        chunk["caseid"] = normalize_identifier(chunk["caseid"])
                    else:
                        chunk["caseid"] = pd.NA
                    selected_rows = chunk[chunk["primaryid"].isin(target_primaryids)]
                    for primaryid, group in selected_rows.groupby("primaryid", sort=True):
                        key = (primaryid, source.year, source.quarter)
                        item = occurrences[key]
                        item["rows"] += len(group)
                        item["caseids"].update(group["caseid"].dropna())
                        if "pt" in group:
                            item["pts"].update(group["pt"].dropna().astype(str).str.strip())
                    case_rows = chunk[chunk["caseid"].isin(target_caseids)]
                    for key, count in case_rows.groupby(["caseid", "primaryid"], dropna=False).size().items():
                        case_primary_counts[key] += int(count)
    return occurrences, case_primary_counts, chunks, rows


def classify_case(selected_found_expected: bool, selected_found_other: bool, other_available: bool,
                  raw_caseid_conflict: bool = False) -> str:
    if raw_caseid_conflict:
        return "raw_caseid_conflict"
    if selected_found_expected:
        return "found_expected_partition_integration_error"
    if selected_found_other:
        return "found_different_partition_provenance_error"
    if other_available:
        return "selected_version_missing_reac_nonselected_version_available"
    return "confirmed_missing_in_raw_reac"


def run_audit(args) -> dict:
    assignments = load_global_assignments(args.global_case_assignments)
    resolved = pd.read_parquet(args.resolved_reac, columns=["primaryid"])
    missing = identify_missing_selected_reac(assignments, resolved, args.year)
    if args.expected_missing is not None and len(missing) != args.expected_missing:
        raise AssertionError(f"expected {args.expected_missing} missing cases, found {len(missing)}")
    if missing["caseid"].duplicated().any() or missing["selected_primaryid"].duplicated().any():
        raise ValueError("missing selected cases are not one-to-one")

    versions = pd.read_parquet(args.global_demo_versions, columns=VERSION_COLUMNS)
    provenance = validate_selected_provenance(missing, versions)
    version_rows = versions[versions["caseid"].isin(set(missing["caseid"]))].copy()
    version_rows["caseid"] = normalize_identifier(version_rows["caseid"])
    version_rows["primaryid"] = normalize_identifier(version_rows["primaryid"])

    sources = select_sources(load_manifest(args.manifest), args.loaded_years, (1, 2, 3, 4))
    source_map = {(s.year, s.quarter): s for s in sources}
    expected_keys = set(zip(missing["selected_source_year"].astype(int), missing["selected_source_quarter"].astype(int)))
    if not expected_keys.issubset(source_map):
        raise ValueError(f"selected source provenance absent from manifest: {sorted(expected_keys - source_map.keys())}")

    raw_paths = [locate_archive(Path(args.raw_dir), s.local_filename) for s in sources]
    hashes_before = {p: sha256_file(p) for p in raw_paths}
    target_pids = set(version_rows["primaryid"].dropna())
    occurrences, case_counts, chunks, rows = scan_reac_sources(
        sources, Path(args.raw_dir), target_pids, set(missing["caseid"]), args.chunk_size,
    )
    hashes_after = {p: sha256_file(p) for p in raw_paths}
    raw_modified = hashes_before != hashes_after
    if raw_modified:
        raise AssertionError("a raw ZIP hash changed during the read-only audit")

    occurrence_rows = []
    for (pid, year, quarter), item in sorted(occurrences.items()):
        source = source_map[(year, quarter)]
        occurrence_rows.append({
            "primaryid": pid, "source_year": year, "source_quarter": quarter,
            "source_archive": repository_path(locate_archive(Path(args.raw_dir), source.local_filename)),
            "reac_member": source.expected_inner_reac.replace("\\", "/"), "row_count": item["rows"],
            "distinct_raw_caseids": "|".join(sorted(item["caseids"])),
            "sample_reaction_values": "|".join(sorted(item["pts"])[:10]),
        })

    case_rows = []
    selected_version_map = {
        (row.caseid, row.primaryid): row
        for row in version_rows[version_rows["is_selected"].eq(True)].itertuples(index=False)
    }
    for row in missing.itertuples(index=False):
        sy, sq = int(row.selected_source_year), int(row.selected_source_quarter)
        source = source_map[(sy, sq)]
        selected_occ = [(key, value) for key, value in occurrences.items() if key[0] == row.selected_primaryid]
        expected_occ = [value for key, value in selected_occ if key[1:] == (sy, sq)]
        other_occ = [value for key, value in selected_occ if key[1:] != (sy, sq)]
        expected_caseids = set().union(*(x["caseids"] for x in expected_occ)) if expected_occ else set()
        conflict = any(cid != row.caseid for cid in expected_caseids)
        alternatives = sorted(
            pid for (cid, pid), count in case_counts.items()
            if cid == row.caseid and pid != row.selected_primaryid and count > 0
        )
        expected_count = sum(x["rows"] for x in expected_occ)
        total_count = sum(x["rows"] for _, x in selected_occ)
        category = classify_case(bool(expected_count), bool(other_occ), bool(alternatives), conflict)
        demo = selected_version_map[(row.caseid, row.selected_primaryid)]
        case_rows.append({
            "caseid": row.caseid, "selected_primaryid": row.selected_primaryid,
            "assigned_report_year": int(row.assigned_report_year),
            "selected_caseversion": int(row.selected_caseversion) if pd.notna(row.selected_caseversion) else "",
            "assignment_method": row.assignment_method,
            "assignment_status": "included", "selected_init_fda_dt": demo.init_fda_dt_raw,
            "selected_fda_dt": demo.fda_dt_raw,
            "selected_source_year": sy, "selected_source_quarter": sq,
            "selected_source_archive": repository_path(locate_archive(Path(args.raw_dir), source.local_filename)),
            "selected_reac_member": source.expected_inner_reac.replace("\\", "/"),
            "found_in_expected_reac_partition": bool(expected_count),
            "expected_partition_match_count": expected_count,
            "expected_partition_raw_caseids": "|".join(sorted(expected_caseids)),
            "expected_partition_reaction_values": "|".join(sorted(set().union(*(x["pts"] for x in expected_occ)))[:10]) if expected_occ else "",
            "found_in_any_reac_partition": bool(total_count), "total_reac_match_count": total_count,
            "found_in_other_partition": bool(other_occ),
            "same_caseid_other_primaryid_has_reac": bool(alternatives),
            "other_primaryids_with_reac": "|".join(alternatives), "classification": category,
            "notes": "Nonselected PRIMARYIDs are diagnostic only and were not substituted.",
        })

    version_report = []
    missing_pid_map = dict(zip(missing["caseid"], missing["selected_primaryid"]))
    for row in version_rows.sort_values(["caseid", "source_year", "source_quarter", "primaryid"]).itertuples(index=False):
        count = case_counts.get((row.caseid, row.primaryid), 0)
        version_report.append({
            "caseid": row.caseid, "primaryid": row.primaryid, "caseversion": row.caseversion_numeric,
            "source_year": int(row.source_year), "source_quarter": int(row.source_quarter),
            "is_selected_primaryid": row.primaryid == missing_pid_map[row.caseid],
            "has_reac_rows": count > 0, "reac_row_count": count,
        })

    out = Path(args.report_dir); out.mkdir(parents=True, exist_ok=True)
    cases_path = out / f"ms4b_missing_reac_cases_{args.year}.csv"
    occurrences_path = out / f"ms4b_missing_reac_partition_occurrences_{args.year}.csv"
    versions_path = out / f"ms4b_missing_reac_case_versions_{args.year}.csv"
    pd.DataFrame(case_rows).to_csv(cases_path, index=False)
    pd.DataFrame(occurrence_rows, columns=["primaryid", "source_year", "source_quarter", "source_archive", "reac_member", "row_count", "distinct_raw_caseids", "sample_reaction_values"]).to_csv(occurrences_path, index=False)
    pd.DataFrame(version_report).to_csv(versions_path, index=False)
    counts = pd.Series([row["classification"] for row in case_rows]).value_counts()
    summary = {
        "analytical_year": args.year, "missing_selected_reac_cases_expected": args.expected_missing,
        "missing_selected_reac_cases_found": len(case_rows),
        "confirmed_missing_in_raw_reac": int(counts.get("confirmed_missing_in_raw_reac", 0)),
        "integration_errors": int(counts.get("found_expected_partition_integration_error", 0)),
        "provenance_errors": int(counts.get("found_different_partition_provenance_error", 0)),
        "selected_version_missing_nonselected_available": int(counts.get("selected_version_missing_reac_nonselected_version_available", 0)),
        "identifier_normalization_errors": int(counts.get("identifier_normalization_error", 0)),
        "unclassified_cases": int(counts.get("unclassified_requires_review", 0)),
        "archives_scanned": len(sources), "raw_rows_scanned_or_chunks_processed": {"rows": rows, "chunks": chunks},
        "absolute_paths_persisted": False, "raw_files_modified": raw_modified,
        "reports": [repository_path(p) for p in (cases_path, occurrences_path, versions_path)],
    }
    json_path = out / f"ms4b_missing_reac_audit_{args.year}.json"
    summary["reports"].append(repository_path(json_path))
    with json_path.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True); stream.write("\n")
    return {"summary": summary, "cases": case_rows, "provenance": provenance}


def parse_args():
    parser = argparse.ArgumentParser(description="Audit selected MS4A cases absent from resolved REAC without modifying data.")
    parser.add_argument("--year", type=int, default=2024)
    parser.add_argument("--global-case-assignments", type=Path, default=DEFAULT_ASSIGNMENTS)
    parser.add_argument("--global-demo-versions", type=Path, default=DEFAULT_VERSIONS)
    parser.add_argument("--resolved-reac", type=Path, default=DEFAULT_RESOLVED)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--loaded-years", type=int, nargs="+", default=[2023, 2024, 2025])
    parser.add_argument("--chunk-size", type=int, default=200_000)
    parser.add_argument("--expected-missing", type=int, default=3)
    return parser.parse_args()


if __name__ == "__main__":
    result = run_audit(parse_args())
    print(json.dumps(result["summary"], indent=2, sort_keys=True))
