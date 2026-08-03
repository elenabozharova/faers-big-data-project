"""Rebuild final annual pair-reaction artifacts from the validated MS4B boundary.

This MS6 orchestrator deliberately does not read raw FAERS tables. It treats
the retained/resolved MS4B Parquet files as immutable inputs, applies the
versioned MS5 mapping, delegates all transformations to the existing
case-level pipeline, validates the result, and atomically publishes one year.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import tempfile
import time
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

import integrate_global_case_assignments as ms4b
import process_faers_case_level as pipeline
from download_faers import sha256_file


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BOUNDARY_ROOT = ROOT / "data/processed/faers_case_level_annual"
DEFAULT_OUTPUT_ROOT = ROOT / "data/processed/faers_case_level_annual_final"
DEFAULT_ASSIGNMENTS = ROOT / "data/interim/case_resolution/global_case_assignments.parquet"
DEFAULT_MAPPING = ROOT / "data/processed/mappings/faers_to_chembl_2023_2025_chembl36_norm_v1_exact_v1.parquet"
DEFAULT_REPORT_DIR = ROOT / "reports/data_quality"
EXPECTED_CASES = {2023: 1_361_395, 2024: 1_308_773, 2025: 1_295_857}
BOUNDARY_FILES = ("retained_cases.parquet", "resolved_drug_rows.parquet", "resolved_reac_rows.parquet")
MISSING_REAC_2024 = {"24346610", "24346611", "24346612"}


def repository_path(path: Path | str) -> str:
    return ms4b.repository_path(path)


def _git_metadata() -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    dirty = bool(subprocess.run(
        ["git", "status", "--porcelain"], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout)
    return commit, dirty


def _observed_memory_mb() -> float | None:
    """Return observed process working-set memory without adding a dependency."""
    if os.name != "nt":
        return None
    try:
        import ctypes
        from ctypes import wintypes

        class ProcessMemoryCounters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        get_process = ctypes.windll.kernel32.GetCurrentProcess
        get_process.restype = wintypes.HANDLE
        process = get_process()
        get_memory = ctypes.windll.psapi.GetProcessMemoryInfo
        get_memory.argtypes = [wintypes.HANDLE, ctypes.POINTER(ProcessMemoryCounters), wintypes.DWORD]
        get_memory.restype = wintypes.BOOL
        ok = get_memory(
            process, ctypes.byref(counters), counters.cb,
        )
        return round(counters.PeakWorkingSetSize / 1024 / 1024, 2) if ok else None
    except Exception:
        return None


def _directory_size(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


class StageRecorder:
    def __init__(self, year: int, output_dir: Path):
        self.year = year
        self.output_dir = output_dir
        self.rows: list[dict] = []

    @contextmanager
    def stage(self, name: str, input_rows: int, output_rows, chunk_count: int = 1):
        started = datetime.now(timezone.utc)
        start_clock = time.perf_counter()
        size_before = _directory_size(self.output_dir)
        record = {}
        try:
            yield record
        finally:
            end = datetime.now(timezone.utc)
            resolved_output_rows = output_rows() if callable(output_rows) else output_rows
            self.rows.append({
                "year": self.year,
                "stage": name,
                "start_time": started.isoformat(),
                "end_time": end.isoformat(),
                "duration_seconds": round(time.perf_counter() - start_clock, 3),
                "input_rows": int(input_rows),
                "output_rows": int(resolved_output_rows),
                "chunk_count": int(record.get("chunk_count", chunk_count)),
                "temporary_storage_bytes": max(0, _directory_size(self.output_dir) - size_before),
                "observed_or_peak_memory_mb": _observed_memory_mb(),
            })


def input_paths(years: list[int], boundary_root: Path, assignments: Path, mapping: Path) -> list[Path]:
    paths = [assignments, mapping]
    for year in years:
        paths.extend(boundary_root / str(year) / name for name in BOUNDARY_FILES)
    raw_root = ROOT / "data/raw/faers"
    for year in years:
        matches = sorted(
            (item for item in raw_root.glob("*.zip") if f"{year}q" in item.name.casefold()),
            key=lambda item: item.name.casefold(),
        )
        if len(matches) != 4:
            raise FileNotFoundError(f"expected four raw ZIP archives for {year}, found {len(matches)}")
        paths.extend(matches)
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"production input missing: {missing[0]}")
    return paths


def hash_inputs(paths: list[Path]) -> dict[str, str]:
    return {repository_path(path): sha256_file(path) for path in paths}


def load_boundary(boundary_root: Path, year: int, assignments: pd.DataFrame):
    source_dir = boundary_root / str(year)
    retained = pd.read_parquet(source_dir / "retained_cases.parquet")
    drug = pd.read_parquet(source_dir / "resolved_drug_rows.parquet")
    reac = pd.read_parquet(source_dir / "resolved_reac_rows.parquet")

    for frame, columns, label in (
        (retained, {"caseid", "retained_primaryid"}, "retained_cases"),
        (drug, {"caseid", "primaryid"}, "resolved_drug_rows"),
        (reac, {"caseid", "primaryid"}, "resolved_reac_rows"),
    ):
        missing = columns - set(frame.columns)
        if missing:
            raise ValueError(f"{year} {label} missing columns: {sorted(missing)}")

    selected = ms4b.selected_cases_for_year(assignments, year)[["caseid", "selected_primaryid"]].copy()
    selected["caseid"] = ms4b.normalize_identifier(selected["caseid"])
    selected["selected_primaryid"] = ms4b.normalize_identifier(selected["selected_primaryid"])
    retained["caseid"] = ms4b.normalize_identifier(retained["caseid"])
    retained["retained_primaryid"] = ms4b.normalize_identifier(retained["retained_primaryid"])
    for frame in (drug, reac):
        frame["caseid"] = ms4b.normalize_identifier(frame["caseid"])
        frame["primaryid"] = ms4b.normalize_identifier(frame["primaryid"])

    if retained["caseid"].duplicated().any():
        raise AssertionError(f"{year}: retained_cases is not unique by CASEID")
    expected = EXPECTED_CASES.get(year, len(selected))
    if len(selected) != expected or len(retained) != expected:
        raise AssertionError(
            f"{year}: expected {expected:,} assignments/retained cases, "
            f"found {len(selected):,}/{len(retained):,}"
        )
    contract = retained.merge(selected, on="caseid", how="outer", indicator=True)
    if not contract["_merge"].eq("both").all():
        raise AssertionError(f"{year}: retained CASEIDs do not exactly match included assignments")
    if not contract["retained_primaryid"].eq(contract["selected_primaryid"]).all():
        raise AssertionError(f"{year}: retained PRIMARYIDs do not match global selections")

    retained_caseids = set(retained["caseid"])
    selected_pairs = set(selected.itertuples(index=False, name=None))
    for frame, label in ((drug, "DRUG"), (reac, "REAC")):
        if not set(frame["caseid"]).issubset(retained_caseids):
            raise AssertionError(f"{year}: resolved {label} contains an unexpected CASEID")
        observed = set(frame[["caseid", "primaryid"]].drop_duplicates().itertuples(index=False, name=None))
        if not observed.issubset(selected_pairs):
            raise AssertionError(f"{year}: resolved {label} contains a nonselected PRIMARYID")
    return retained, drug, reac


def _link_boundary(source_root: Path, year: int, output_dir: Path) -> None:
    for name in BOUNDARY_FILES:
        source = (source_root / str(year) / name).resolve()
        destination = output_dir / name
        os.link(source, destination)


def build_case_usability(year: int, retained: pd.DataFrame, long_drugs: pd.DataFrame,
                         case_drugs: pd.DataFrame) -> dict:
    retained_ids = set(retained["caseid"])
    drug_caseids = set(long_drugs["caseid"])
    mapped_counts = case_drugs.groupby("caseid")["mapped_drug_id"].nunique()
    no_mapped = retained_ids - set(mapped_counts.index)
    statuses = long_drugs.groupby("caseid")["mapping_status"].agg(set)
    combo_status = "mapped_combination_components"
    ambiguous = {"unmapped_ambiguous", "unmapped_combination_ambiguous"}
    usable_name_cases = set(long_drugs.loc[long_drugs["normalized_drug_name"].fillna("").ne(""), "caseid"])
    return {
        "year": year,
        "retained_cases": len(retained_ids),
        "cases_with_drug_rows": len(drug_caseids),
        "cases_with_zero_single_id_mapped_drugs": len(no_mapped),
        "cases_with_one_single_id_mapped_drug": int((mapped_counts == 1).sum()),
        "cases_with_at_least_two_single_id_mapped_drugs": int((mapped_counts >= 2).sum()),
        "cases_with_combination_only_mappings": sum(
            caseid in no_mapped and combo_status in values for caseid, values in statuses.items()
        ),
        "cases_with_ambiguous_only_mappings": sum(
            caseid in no_mapped and bool(values & ambiguous) for caseid, values in statuses.items()
        ),
        "cases_with_no_usable_drug_name": len(drug_caseids - usable_name_cases),
    }


def _atomic_publish(staged_year: Path, final_year: Path, overwrite: bool) -> None:
    if final_year.exists() and not overwrite:
        raise FileExistsError(f"output exists; use --overwrite: {final_year}")
    backup = final_year.with_name(f".{final_year.name}.backup-{uuid.uuid4().hex}")
    try:
        if final_year.exists():
            os.replace(final_year, backup)
        os.replace(staged_year, final_year)
    except Exception:
        if backup.exists() and not final_year.exists():
            os.replace(backup, final_year)
        raise
    else:
        if backup.exists():
            shutil.rmtree(backup)


def rebuild_year(year: int, *, boundary_root: Path, output_root: Path, mapping_path: Path,
                 assignments: pd.DataFrame, pair_role_filter: str,
                 max_rows_per_join_chunk: int, overwrite: bool) -> tuple[dict, list[dict]]:
    final_year = output_root / str(year)
    if final_year.exists() and not overwrite:
        raise FileExistsError(f"output exists; use --overwrite: {final_year}")
    output_root.mkdir(parents=True, exist_ok=True)
    staging_root = Path(tempfile.mkdtemp(prefix=f".{output_root.name}-{year}-", dir=output_root.parent))
    staged_year = staging_root / str(year)
    staged_year.mkdir()
    recorder = StageRecorder(year, staged_year)
    started = time.perf_counter()
    try:
        retained, resolved_drug, resolved_reac = load_boundary(boundary_root, year, assignments)
        _link_boundary(boundary_root, year, staged_year)

        reactions_holder = {}
        with recorder.stage("case_reaction_construction", len(resolved_reac),
                            lambda: len(reactions_holder["data"])):
            reactions, reaction_diag = pipeline.build_case_reactions(resolved_reac)
            reactions_holder["data"] = reactions
            pipeline.run_case_reactions_validation(resolved_reac, reactions)
            pipeline.save_case_reactions_artifacts(staged_year, reactions, reaction_diag)

        drug_holder = {}
        with recorder.stage("case_drug_mapping", len(resolved_drug),
                            lambda: len(drug_holder["long"])):
            long_drugs, case_drugs, mapping_diag, unmatched = pipeline.build_case_drugs(
                resolved_drug, str(mapping_path), include_fuzzy=False,
            )
            drug_holder["long"] = long_drugs
            pipeline.run_case_drugs_validation(resolved_drug, long_drugs, case_drugs)
            pipeline.save_case_drugs_artifacts(staged_year, long_drugs, case_drugs, mapping_diag, unmatched)

        recorder.rows.append({
            "year": year, "stage": "case_drug_deduplication",
            "start_time": recorder.rows[-1]["start_time"], "end_time": recorder.rows[-1]["end_time"],
            "duration_seconds": 0.0, "input_rows": mapping_diag["n_mapped_rows"],
            "output_rows": len(case_drugs), "chunk_count": 1, "temporary_storage_bytes": 0,
            "observed_or_peak_memory_mb": _observed_memory_mb(),
        })
        usability = build_case_usability(year, retained, long_drugs, case_drugs)

        pairs_holder = {}
        with recorder.stage("pair_generation", len(case_drugs), lambda: len(pairs_holder["pairs"])):
            pairs, pair_diag, high_poly = pipeline.build_case_drug_pairs(case_drugs)
            pairs_holder["pairs"] = pairs
            pipeline.run_case_drug_pairs_validation(case_drugs, pairs)
            pipeline.save_case_drug_pairs_artifacts(staged_year, pairs, pair_diag, high_poly)
            role_audit, role_diag = pipeline.build_case_drug_pair_role_audit(pairs, case_drugs)
            pipeline.run_drug_pair_role_audit_validation(pairs, role_audit, role_diag)
            pipeline.save_drug_pair_role_audit_artifacts(staged_year, role_audit, role_diag)

        selected_holder = {}
        with recorder.stage("pair_role_filtering", len(pairs), lambda: len(selected_holder["selected"])):
            selected, filter_diag = pipeline.apply_pair_role_filter(pairs, role_audit, pair_role_filter)
            selected_holder["selected"] = selected
            pipeline.run_pair_role_filter_validation(pairs, role_audit, selected, pair_role_filter)
            pipeline.save_pair_role_filter_artifacts(staged_year, selected, filter_diag)

        join_holder = {"rows": 0}
        pair_reaction_path = staged_year / "case_pair_reactions.parquet"
        with recorder.stage("pair_reaction_join", len(selected) + len(reactions),
                            lambda: join_holder["rows"]) as stage:
            total_rows, counts, case_stats, row_validation = pipeline.build_case_pair_reactions_chunked(
                selected, reactions, str(pair_reaction_path), pair_role_filter,
                max_rows_per_chunk=max_rows_per_join_chunk,
            )
            join_holder["rows"] = total_rows
            stage["chunk_count"] = len(pipeline._plan_caseid_chunks(case_stats, max_rows_per_join_chunk))
            pipeline.run_case_pair_reactions_validation(row_validation, pair_role_filter)

        with recorder.stage("support_aggregation", total_rows, lambda: len(counts)):
            pair_reaction_diag = pipeline.compute_pair_reaction_diagnostics(
                selected, reactions, total_rows, counts, case_stats, pair_role_filter,
            )
            pipeline.run_drug_pair_reaction_counts_validation(counts, str(pair_reaction_path))
            counts = counts.copy()
            counts.insert(0, "year", year)
            pipeline.save_pair_reaction_artifacts(staged_year, counts, pair_reaction_diag)

        retained_ids = set(retained["caseid"])
        downstream_sets = {
            "case_reaction": set(reactions["caseid"]), "case_drug": set(case_drugs["caseid"]),
            "case_pair": set(pairs["caseid"]), "pair_reaction": set(case_stats["caseid"]),
        }
        subset_violations = sum(len(values - retained_ids) for values in downstream_sets.values())
        resolved_drug_cases = set(resolved_drug["caseid"])
        resolved_reac_cases = set(resolved_reac["caseid"])
        selected_pids = set(retained["retained_primaryid"])
        nonselected = len(
            (set(resolved_drug["primaryid"]) | set(resolved_reac["primaryid"])) - selected_pids
        )
        if subset_violations or nonselected:
            raise AssertionError(f"{year}: downstream subset/nonselected violations")

        missing_audit = {}
        if year == 2024:
            missing_audit = {
                "retained": MISSING_REAC_2024 <= retained_ids,
                "in_resolved_drug": MISSING_REAC_2024 <= resolved_drug_cases,
                "absent_resolved_reac": not bool(MISSING_REAC_2024 & resolved_reac_cases),
                "absent_case_reactions": not bool(MISSING_REAC_2024 & downstream_sets["case_reaction"]),
                "absent_pair_reactions": not bool(MISSING_REAC_2024 & downstream_sets["pair_reaction"]),
            }
            if not all(missing_audit.values()):
                raise AssertionError(f"2024 missing-REAC audit failed: {missing_audit}")

        source_hashes = {
            repository_path(boundary_root / str(year) / name):
            sha256_file(boundary_root / str(year) / name) for name in BOUNDARY_FILES
        }
        mapping_hash = sha256_file(mapping_path)
        summary = {
            "year": year,
            "case_resolution": {
                "method": "global_initial_year_frozen",
                "expected_included_assignments": EXPECTED_CASES.get(year, len(retained)),
                "retained_cases": len(retained),
                "nonselected_primaryids_excluded": nonselected,
                "missing_drug_cases": len(retained_ids - resolved_drug_cases),
                "missing_reac_cases": len(retained_ids - resolved_reac_cases),
            },
            "reused_boundary": {
                "source_root": repository_path(boundary_root),
                "files": source_hashes,
                "publication_method": "same-volume hard links",
            },
            "drug_mapping": {
                "artifact": repository_path(mapping_path), "sha256": mapping_hash,
                "chembl_version": 36, "normalization_version": "norm_v1",
                "mapping_scope": "exact_only", "fuzzy_mappings_accepted": False,
            },
            "case_usability": usability,
            "pair_reaction_outputs": {
                "case_drugs": len(case_drugs), "unique_mapped_drugs": case_drugs["mapped_drug_id"].nunique(),
                "case_drug_pairs": len(pairs), "case_drug_pairs_selected": len(selected),
                "case_pair_reactions": total_rows, "unique_pair_reaction_patterns": len(counts),
                "min3_patterns": int((counts["n_cases"] >= 3).sum()),
                "min5_patterns": int((counts["n_cases"] >= 5).sum()),
                "min10_patterns": int((counts["n_cases"] >= 10).sum()),
            },
            "validation": {
                "resolved_drug_caseids": len(resolved_drug_cases),
                "resolved_reac_caseids": len(resolved_reac_cases),
                "case_reaction_caseids": len(downstream_sets["case_reaction"]),
                "case_drug_caseids": len(downstream_sets["case_drug"]),
                "case_pair_caseids": len(downstream_sets["case_pair"]),
                "pair_reaction_caseids": len(downstream_sets["pair_reaction"]),
                "unexpected_caseids": 0, "unexpected_primaryids": 0,
                "nonselected_primaryids": nonselected,
                "downstream_caseid_subset_violations": subset_violations,
                "pair_order_violations": int((pairs["drug_a"] >= pairs["drug_b"]).sum()),
                "self_pair_violations": int((pairs["drug_a"] == pairs["drug_b"]).sum()),
                "pair_reaction_duplicate_rows": row_validation["n_duplicate_keys"],
                "support_count_violations": 0,
                "missing_reac_audit": missing_audit,
            },
            "performance": recorder.rows,
            "runtime_seconds": round(time.perf_counter() - started, 3),
            "output_path": repository_path(final_year),
        }
        with (staged_year / "annual_summary.json").open("w", encoding="utf-8") as stream:
            json.dump(summary, stream, indent=2, sort_keys=True)
            stream.write("\n")
        _atomic_publish(staged_year, final_year, overwrite)
        return summary, recorder.rows
    finally:
        if staging_root.exists():
            shutil.rmtree(staging_root)


def write_reports(summaries: list[dict], performance: list[dict], *, report_dir: Path,
                  assignments_path: Path, mapping_path: Path, input_hashes_before: dict,
                  input_hashes_after: dict, output_root: Path, pair_role_filter: str,
                  max_rows_per_join_chunk: int) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    years = [item["year"] for item in summaries]
    cross_year_overlap = 0
    ms4_rows = []
    ms6_rows = []
    usability_rows = []
    for item in summaries:
        year = item["year"]
        case = item["case_resolution"]
        valid = item["validation"]
        outputs = item["pair_reaction_outputs"]
        usability = item["case_usability"]
        usability_rows.append(usability)
        ms4_rows.append({
            "year": year,
            "expected_included_assignments": case["expected_included_assignments"],
            "retained_cases": case["retained_cases"],
            "resolved_drug_caseids": valid["resolved_drug_caseids"],
            "resolved_reac_caseids": valid["resolved_reac_caseids"],
            "case_reaction_caseids": valid["case_reaction_caseids"],
            "case_drug_caseids": valid["case_drug_caseids"],
            "case_pair_caseids": valid["case_pair_caseids"],
            "pair_reaction_caseids": valid["pair_reaction_caseids"],
            "missing_drug_cases": case["missing_drug_cases"],
            "missing_reac_cases": case["missing_reac_cases"],
            "unexpected_caseids": valid["unexpected_caseids"],
            "unexpected_primaryids": valid["unexpected_primaryids"],
            "nonselected_primaryids": valid["nonselected_primaryids"],
            "cross_year_caseid_overlap": cross_year_overlap,
            "downstream_caseid_subset_violations": valid["downstream_caseid_subset_violations"],
            "validation_status": "passed",
        })
        ms6_rows.append({
            "year": year, "retained_cases": case["retained_cases"],
            "cases_with_at_least_two_mapped_drugs": usability["cases_with_at_least_two_single_id_mapped_drugs"],
            "case_drug_rows": outputs["case_drugs"], "unique_mapped_drugs": outputs["unique_mapped_drugs"],
            "case_drug_pair_rows": outputs["case_drug_pairs"],
            "case_drug_pair_selected_rows": outputs["case_drug_pairs_selected"],
            "case_pair_reaction_rows": outputs["case_pair_reactions"],
            "unique_pair_reaction_patterns": outputs["unique_pair_reaction_patterns"],
            "min3_patterns": outputs["min3_patterns"], "min5_patterns": outputs["min5_patterns"],
            "min10_patterns": outputs["min10_patterns"],
            "pair_order_violations": valid["pair_order_violations"],
            "self_pair_violations": valid["self_pair_violations"],
            "pair_reaction_duplicate_rows": valid["pair_reaction_duplicate_rows"],
            "support_count_violations": valid["support_count_violations"],
            "validation_status": "passed",
        })

    pd.DataFrame(ms4_rows).to_csv(report_dir / "ms4b_downstream_validation_by_year.csv", index=False)
    pd.DataFrame(ms6_rows).to_csv(report_dir / "ms6_annual_pair_reaction_validation.csv", index=False)
    pd.DataFrame(usability_rows).to_csv(report_dir / "ms6_case_mapping_usability_by_year.csv", index=False)
    pd.DataFrame(performance).to_csv(report_dir / "ms6_stage_performance_by_year.csv", index=False)

    commit, dirty = _git_metadata()
    common = {"created_at": datetime.now(timezone.utc).isoformat(), "code_commit": commit,
              "worktree_dirty": dirty, "output_root": repository_path(output_root)}
    ms4_summary = common | {
        "years": years, "case_resolution_method": "global_initial_year_frozen",
        "global_assignments_input": repository_path(assignments_path),
        "mapping_artifact": repository_path(mapping_path),
        "full_downstream_rebuild_completed": True, "all_ms4b_invariants_passed": True,
        "raw_inputs_modified": False, "resolved_inputs_modified": False,
        "temporary_files_remaining": 0, "absolute_paths_persisted": False,
    }
    ms6_summary = common | {
        "years": years, "mapping_artifact": repository_path(mapping_path),
        "mapping_artifact_sha256": sha256_file(mapping_path),
        "pair_role_filter": pair_role_filter,
        "max_rows_per_join_chunk": max_rows_per_join_chunk,
        "observed_peak_memory_mb": max(
            (row["observed_or_peak_memory_mb"] for row in performance
             if row["observed_or_peak_memory_mb"] is not None),
            default=None,
        ),
        "all_ms6_invariants_passed": True, "production_rebuild_completed": True,
        "temporary_files_remaining": 0, "absolute_paths_persisted": False,
    }
    integrity = common | {
        "input_hashes_before": input_hashes_before, "input_hashes_after": input_hashes_after,
        "all_input_hashes_unchanged": input_hashes_before == input_hashes_after,
        "raw_inputs_modified": False, "resolved_inputs_modified": False,
        "temporary_files_remaining": 0, "absolute_paths_persisted": False,
    }
    for name, payload in (
        ("ms4b_downstream_validation_summary.json", ms4_summary),
        ("ms6_annual_pair_reaction_summary.json", ms6_summary),
        ("ms6_production_input_integrity.json", integrity),
    ):
        with (report_dir / name).open("w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")

    legacy = []
    previous = {2023: 1_544_991, 2024: 1_484_350, 2025: 1_469_305}
    for item in summaries:
        year = item["year"]
        outputs = item["pair_reaction_outputs"]
        legacy.append({
            "year": year, "legacy_retained_cases": previous.get(year, ""),
            "global_retained_cases": item["case_resolution"]["retained_cases"],
            "removed_annual_appearances": previous.get(year, 0) - item["case_resolution"]["retained_cases"],
            "legacy_case_reactions": "unavailable", "global_case_reactions": item["validation"]["case_reaction_caseids"],
            "legacy_case_drugs": "unavailable", "global_case_drugs": outputs["case_drugs"],
            "legacy_case_drug_pairs": "unavailable", "global_case_drug_pairs": outputs["case_drug_pairs"],
            "legacy_case_pair_reactions": "unavailable", "global_case_pair_reactions": outputs["case_pair_reactions"],
            "legacy_unique_pair_reaction_patterns": "unavailable",
            "global_unique_pair_reaction_patterns": outputs["unique_pair_reaction_patterns"],
            "legacy_min3_patterns": "unavailable", "global_min3_patterns": outputs["min3_patterns"],
            "legacy_min5_patterns": "unavailable", "global_min5_patterns": outputs["min5_patterns"],
            "legacy_min10_patterns": "unavailable", "global_min10_patterns": outputs["min10_patterns"],
        })
    pd.DataFrame(legacy).to_csv(report_dir / "ms4b_ms6_legacy_global_downstream_comparison.csv", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reuse the immutable MS4B boundary and atomically build final MS6 annual outputs."
    )
    parser.add_argument("--years", nargs="+", type=int, required=True)
    parser.add_argument("--boundary-root", type=Path, default=DEFAULT_BOUNDARY_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--global-case-assignments", type=Path, default=DEFAULT_ASSIGNMENTS)
    parser.add_argument("--mapping-file", type=Path, default=DEFAULT_MAPPING)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--pair-role-filter", choices=pipeline.PAIR_ROLE_FILTERS,
                        default=pipeline.DEFAULT_PAIR_ROLE_FILTER)
    parser.add_argument("--max-rows-per-join-chunk", type=int,
                        default=pipeline.DEFAULT_MAX_ROWS_PER_JOIN_CHUNK)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    years = sorted(set(args.years))
    if args.out_dir.resolve() == args.boundary_root.resolve():
        raise ValueError("output root must differ from the immutable boundary root")
    paths = input_paths(years, args.boundary_root, args.global_case_assignments, args.mapping_file)
    report_outputs = [
        args.report_dir / name for name in (
            "ms4b_downstream_validation_by_year.csv",
            "ms4b_downstream_validation_summary.json",
            "ms6_case_mapping_usability_by_year.csv",
            "ms6_annual_pair_reaction_validation.csv",
            "ms6_annual_pair_reaction_summary.json",
            "ms6_stage_performance_by_year.csv",
            "ms6_production_input_integrity.json",
            "ms4b_ms6_legacy_global_downstream_comparison.csv",
        )
    ]
    combined_summary = args.out_dir / "annual_summary_all_years.json"
    if not args.overwrite:
        existing = [path for path in [combined_summary, *report_outputs] if path.exists()]
        if existing:
            raise FileExistsError(f"output exists; use --overwrite: {existing[0]}")
    before = hash_inputs(paths)
    assignments = ms4b.load_global_assignments(args.global_case_assignments)
    included = assignments[
        assignments["assignment_status"].eq("included")
        & assignments["assigned_report_year"].isin(years)
    ]
    if included["caseid"].duplicated().any():
        raise AssertionError("included CASEID appears in multiple requested analytical years")

    summaries = []
    performance = []
    for year in years:
        summary, stages = rebuild_year(
            year, boundary_root=args.boundary_root, output_root=args.out_dir,
            mapping_path=args.mapping_file, assignments=assignments,
            pair_role_filter=args.pair_role_filter,
            max_rows_per_join_chunk=args.max_rows_per_join_chunk,
            overwrite=args.overwrite,
        )
        summaries.append(summary)
        performance.extend(stages)

    after = hash_inputs(paths)
    if before != after:
        raise AssertionError("one or more production input hashes changed during rebuild")
    write_reports(
        summaries, performance, report_dir=args.report_dir,
        assignments_path=args.global_case_assignments, mapping_path=args.mapping_file,
        input_hashes_before=before, input_hashes_after=after, output_root=args.out_dir,
        pair_role_filter=args.pair_role_filter,
        max_rows_per_join_chunk=args.max_rows_per_join_chunk,
    )
    with (args.out_dir / "annual_summary_all_years.json").open("w", encoding="utf-8") as stream:
        json.dump(summaries, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"years": years, "output_root": repository_path(args.out_dir),
                      "all_invariants_passed": True}, indent=2))


if __name__ == "__main__":
    main()
