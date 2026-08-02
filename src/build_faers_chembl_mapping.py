"""Build the versioned exact-only FAERS vocabulary to ChEMBL mapping."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sqlite3
import subprocess
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from build_multi_year_drug_vocabulary import (
    AnnualDrugInput, DEFAULT_INPUT_ROOT, INVALID_NORMALIZED_SENTINELS,
    REQUIRED_COLUMNS, parse_years, resolve_annual_drug_inputs,
)
from download_faers import sha256_file
from integrate_global_case_assignments import normalize_identifier, repository_path
from normalize_drug_names import normalize_drug_string, split_combo
from process_faers_case_level import _choose_raw_drug_name


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_VOCABULARY = ROOT / "data/interim/drug_mapping/multi_year_drug_vocabulary_2023_2025_norm_v1.parquet"
DEFAULT_INDEX = ROOT / "data/interim/drug_mapping/chembl36_name_index_norm_v1_index_v1.parquet"
DEFAULT_OUTPUT = ROOT / "data/processed/mappings/faers_to_chembl_2023_2025_chembl36_norm_v1_exact_v1.parquet"
DEFAULT_REPORT_DIR = ROOT / "reports/data_quality"


def _git_metadata() -> tuple[str, bool]:
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True, check=True).stdout)
    return commit, dirty


def load_name_index(path: Path, *, chembl_version: int, normalization_version: str) -> tuple[dict, dict]:
    required = {"normalized_chembl_name", "chembl_id", "name_type", "chembl_version", "normalization_version"}
    frame = pd.read_parquet(path, columns=list(required))
    missing = required - set(frame.columns)
    if missing: raise ValueError(f"name index missing columns: {sorted(missing)}")
    if set(frame["chembl_version"].dropna().astype(int)) != {chembl_version}:
        raise ValueError("name index ChEMBL version mismatch")
    if set(frame["normalization_version"].dropna().astype(str)) != {normalization_version}:
        raise ValueError("name index normalization version mismatch")
    indexes = {"preferred": defaultdict(set), "synonym": defaultdict(set)}
    for row in frame.itertuples(index=False):
        if row.name_type not in indexes: raise ValueError(f"invalid index name_type: {row.name_type}")
        indexes[row.name_type][row.normalized_chembl_name].add(str(row.chembl_id))
    return ({key: sorted(values) for key, values in indexes["preferred"].items()},
            {key: sorted(values) for key, values in indexes["synonym"].items()})


def exact_component(name: str, preferred: dict, synonyms: dict) -> dict:
    candidates = preferred.get(name)
    method = "exact_pref"
    if candidates is None:
        candidates = synonyms.get(name)
        method = "exact_synonym"
    candidates = sorted(set(candidates or []))
    if not candidates:
        return {"component": name, "status": "unmapped_no_match", "method": "unmatched",
                "chembl_id": None, "candidates": []}
    if len(candidates) > 1:
        return {"component": name, "status": "unmapped_ambiguous", "method": method,
                "chembl_id": None, "candidates": candidates}
    return {"component": name, "status": "mapped", "method": method,
            "chembl_id": candidates[0], "candidates": candidates}


def map_normalized_name(name: str, preferred: dict, synonyms: dict) -> dict:
    if "|" not in name:
        result = exact_component(name, preferred, synonyms)
        return {"chembl_id": result["chembl_id"], "mapping_status": result["status"],
                "mapping_method": result["method"],
                "mapping_score": 1.0 if result["method"] != "unmatched" else None,
                "mapping_candidates": result["candidates"], "candidate_count": len(result["candidates"]),
                "is_ambiguous": result["status"] == "unmapped_ambiguous",
                "component_chembl_ids": [], "component_mapping_details": "[]"}
    details = [exact_component(component, preferred, synonyms) for component in split_combo(name)]
    has_ambiguous = any(item["status"] == "unmapped_ambiguous" for item in details)
    has_unmatched = any(item["status"] == "unmapped_no_match" for item in details)
    if has_ambiguous:
        status = "unmapped_combination_ambiguous"
    elif has_unmatched:
        status = "unmapped_combination_partial"
    else:
        status = "mapped_combination_components"
    ids = [item["chembl_id"] for item in details if item["chembl_id"]]
    return {"chembl_id": None, "mapping_status": status,
            "mapping_method": "combo_all_exact" if status == "mapped_combination_components" else "combo_exact_components",
            "mapping_score": 1.0 if status == "mapped_combination_components" else None,
            "mapping_candidates": [], "candidate_count": 0, "is_ambiguous": has_ambiguous,
            "component_chembl_ids": ids,
            "component_mapping_details": json.dumps(details, sort_keys=True, separators=(",", ":"))}


def build_exact_mapping(vocabulary: pd.DataFrame, preferred: dict, synonyms: dict, *,
                        chembl_version: int, normalization_version: str,
                        vocabulary_artifact_version: str, mapping_artifact_version: str,
                        name_index_artifact_version: str) -> pd.DataFrame:
    required = {"normalized_drug_name", "first_seen_year", "last_seen_year", "years_seen", "year_count",
                "raw_row_count", "distinct_case_count", "distinct_primaryid_count", "distinct_raw_name_count", "sample_raw_names"}
    missing = required - set(vocabulary.columns)
    if missing: raise ValueError(f"vocabulary missing columns: {sorted(missing)}")
    if not vocabulary["normalized_drug_name"].is_unique: raise ValueError("vocabulary names are not unique")
    outcomes = pd.DataFrame([map_normalized_name(name, preferred, synonyms)
                             for name in vocabulary["normalized_drug_name"]])
    result = pd.concat([vocabulary.reset_index(drop=True), outcomes], axis=1)
    result["chembl_version"] = chembl_version; result["normalization_version"] = normalization_version
    result["vocabulary_artifact_version"] = vocabulary_artifact_version
    result["mapping_artifact_version"] = mapping_artifact_version
    result["name_index_artifact_version"] = name_index_artifact_version
    leading = ["normalized_drug_name", "chembl_id", "mapping_status", "mapping_method", "mapping_score",
               "mapping_candidates", "candidate_count", "is_ambiguous", "component_chembl_ids",
               "component_mapping_details"]
    result = result[leading + [column for column in result.columns if column not in leading]]
    result = result.sort_values("normalized_drug_name", kind="stable").reset_index(drop=True)
    if len(result) != len(vocabulary) or not result["normalized_drug_name"].is_unique:
        raise AssertionError("mapping is not one row per vocabulary key")
    if (result["mapping_status"].eq("mapped") & result["chembl_id"].isna()).any():
        raise AssertionError("mapped rows require chembl_id")
    if (result["mapping_status"].ne("mapped") & result["chembl_id"].notna()).any():
        raise AssertionError("non-single-ID rows must not expose chembl_id")
    return result


def build_coverage(mapping: pd.DataFrame, annual_inputs: list[AnnualDrugInput], *, batch_size: int = 250_000):
    lookup = mapping.set_index("normalized_drug_name")[["mapping_status", "mapping_method", "chembl_id"]].to_dict("index")
    annual_metrics = []
    status_rows = []
    requested_years = {item.analytical_year for item in annual_inputs}
    scoped_mapping = mapping[mapping["years_seen"].map(lambda values: bool(set(values) & requested_years))]
    with tempfile.TemporaryDirectory(prefix="ms5_mapping_coverage_") as directory:
        database = sqlite3.connect(Path(directory) / "status_cases.sqlite")
        database.execute("CREATE TABLE status_cases (year INTEGER,status TEXT,method TEXT,caseid TEXT,"
                         "PRIMARY KEY(year,status,method,caseid)) WITHOUT ROWID")
        for item in sorted(annual_inputs, key=lambda value: value.analytical_year):
            rows_by = Counter(); names_by = defaultdict(set); eligible_names = set(); eligible_rows = invalid_rows = 0
            parquet = pq.ParquetFile(item.path)
            for batch in parquet.iter_batches(batch_size=batch_size, columns=list(REQUIRED_COLUMNS)):
                frame = batch.to_pandas(); frame["caseid"] = normalize_identifier(frame["caseid"])
                chosen = _choose_raw_drug_name(frame); normalized = chosen.map(normalize_drug_string)
                valid = normalized.ne("") & ~normalized.str.casefold().isin(INVALID_NORMALIZED_SENTINELS)
                invalid_rows += int((~valid).sum()); frame = frame.loc[valid, ["caseid"]].copy()
                frame["normalized"] = normalized[valid]; eligible_rows += len(frame); eligible_names.update(frame["normalized"])
                frame["status"] = frame["normalized"].map(lambda value: lookup[value]["mapping_status"])
                frame["method"] = frame["normalized"].map(lambda value: lookup[value]["mapping_method"])
                for key, count in frame.groupby(["status", "method"]).size().items(): rows_by[key] += int(count)
                for key, values in frame.groupby(["status", "method"])["normalized"]: names_by[key].update(values)
                database.executemany("INSERT OR IGNORE INTO status_cases VALUES (?,?,?,?)",
                    ((item.analytical_year, status, method, caseid) for status, method, caseid
                     in frame[["status", "method", "caseid"]].drop_duplicates().itertuples(index=False, name=None)))
                database.commit()
            year_mapping = mapping[mapping["years_seen"].map(lambda values: item.analytical_year in values)]
            if eligible_names != set(year_mapping["normalized_drug_name"]):
                raise AssertionError(f"year {item.analytical_year}: annual normalized names do not match vocabulary")
            annual_metrics.append(_coverage_row(str(item.analytical_year), year_mapping, eligible_rows, invalid_rows, rows_by))
            for (status, method), count in sorted(rows_by.items()):
                cases = database.execute("SELECT COUNT(*) FROM status_cases WHERE year=? AND status=? AND method=?",
                                         (item.analytical_year, status, method)).fetchone()[0]
                status_rows.append({"year": item.analytical_year, "mapping_status": status, "mapping_method": method,
                    "unique_normalized_names": len(names_by[(status, method)]), "raw_row_count": count,
                    "distinct_case_count": cases})
        database.close()
    all_rows = Counter()
    for row in status_rows: all_rows[(row["mapping_status"], row["mapping_method"])] += row["raw_row_count"]
    annual_metrics.append(_coverage_row("ALL", scoped_mapping, sum(row["eligible_drug_rows"] for row in annual_metrics),
                                        sum(row["invalid_or_ineligible_rows"] for row in annual_metrics), all_rows))
    for (status, method), count in sorted(all_rows.items()):
        subset = scoped_mapping[(scoped_mapping["mapping_status"] == status) & (scoped_mapping["mapping_method"] == method)]
        all_cases = sum(row["distinct_case_count"] for row in status_rows
                        if row["mapping_status"] == status and row["mapping_method"] == method)
        status_rows.append({"year": "ALL", "mapping_status": status, "mapping_method": method,
            "unique_normalized_names": len(subset), "raw_row_count": count,
            "distinct_case_count": all_cases})
    return pd.DataFrame(annual_metrics), pd.DataFrame(status_rows)


def _coverage_row(year: str, names: pd.DataFrame, eligible_rows: int, invalid_rows: int, rows_by: Counter) -> dict:
    single = names[names["mapping_status"] == "mapped"]
    single_rows = sum(count for (status, _), count in rows_by.items() if status == "mapped")
    combo = names[names["mapping_status"] == "mapped_combination_components"]
    combo_rows = sum(count for (status, _), count in rows_by.items() if status == "mapped_combination_components")
    ambiguous = names[names["mapping_status"] == "unmapped_ambiguous"]
    unmatched = names[names["mapping_status"] == "unmapped_no_match"]
    return {"year": year, "eligible_drug_rows": eligible_rows, "invalid_or_ineligible_rows": invalid_rows,
        "eligible_unique_normalized_names": len(names), "single_id_mapped_rows": single_rows,
        "single_id_unmapped_rows": eligible_rows-single_rows, "single_id_mapped_unique_names": len(single),
        "single_id_unmapped_unique_names": len(names)-len(single),
        "single_id_row_mapping_percentage": round(100*single_rows/eligible_rows, 4) if eligible_rows else 0,
        "single_id_name_mapping_percentage": round(100*len(single)/len(names), 4) if len(names) else 0,
        "exact_pref_unique_names": int((single["mapping_method"] == "exact_pref").sum()),
        "exact_synonym_unique_names": int((single["mapping_method"] == "exact_synonym").sum()),
        "mapped_combination_unique_names": len(combo), "mapped_combination_rows": combo_rows,
        "ambiguous_unique_names": len(ambiguous),
        "ambiguous_rows": sum(count for (status,_),count in rows_by.items() if status == "unmapped_ambiguous"),
        "unmatched_unique_names": len(unmatched),
        "unmatched_rows": sum(count for (status,_),count in rows_by.items() if status == "unmapped_no_match"),
        "partial_combination_unique_names": int((names["mapping_status"] == "unmapped_combination_partial").sum()),
        "ambiguous_combination_unique_names": int((names["mapping_status"] == "unmapped_combination_ambiguous").sum())}


def write_outputs(mapping: pd.DataFrame, coverage: pd.DataFrame, status: pd.DataFrame, *, output: Path,
                  report_dir: Path, vocabulary_path: Path, index_path: Path, requested_years: list[int],
                  normalization_version: str, mapping_artifact_version: str, overwrite: bool) -> dict:
    reports = {"coverage": report_dir/"ms5_chembl_mapping_coverage_by_year.csv",
        "status": report_dir/"ms5_chembl_mapping_status_by_year.csv",
        "unmatched": report_dir/"ms5_chembl_unmatched_names_top.csv",
        "ambiguous": report_dir/"ms5_chembl_ambiguous_mappings.csv",
        "combination": report_dir/"ms5_chembl_combination_mapping.csv",
        "summary": report_dir/"ms5_chembl_mapping_summary.json"}
    paths = [output, *reports.values()]
    if not overwrite:
        existing = [path for path in paths if path.exists()]
        if existing: raise FileExistsError(f"output exists; use --overwrite: {existing[0]}")
    output.parent.mkdir(parents=True, exist_ok=True); report_dir.mkdir(parents=True, exist_ok=True)
    unmatched = mapping[mapping["mapping_status"] == "unmapped_no_match"][["normalized_drug_name","raw_row_count",
        "distinct_case_count","first_seen_year","last_seen_year","years_seen","sample_raw_names"]].sort_values(
            ["raw_row_count","normalized_drug_name"], ascending=[False,True])
    ambiguous = mapping[mapping["mapping_status"] == "unmapped_ambiguous"][[
        "normalized_drug_name","mapping_method","candidate_count","mapping_candidates","raw_row_count",
        "distinct_case_count","sample_raw_names"]]
    combinations = mapping[mapping["normalized_drug_name"].str.contains("|",regex=False)][[
        "normalized_drug_name","mapping_status","component_chembl_ids","component_mapping_details",
        "raw_row_count","distinct_case_count","sample_raw_names"]].copy()
    combinations.insert(2,"component_count",combinations["normalized_drug_name"].str.count(r"\|")+1)
    all_row = coverage[coverage["year"].astype(str)=="ALL"].iloc[0]
    commit, dirty = _git_metadata()
    summary = {"requested_years": requested_years, "vocabulary_input": repository_path(vocabulary_path),
        "vocabulary_sha256": sha256_file(vocabulary_path), "vocabulary_row_count": len(mapping),
        "chembl_name_index": repository_path(index_path), "chembl_name_index_sha256": sha256_file(index_path),
        "chembl_version": 36, "normalization_version": normalization_version,
        "mapping_artifact_version": mapping_artifact_version,
        "single_id_mapped_unique_names": int(all_row.single_id_mapped_unique_names),
        "single_id_mapped_rows": int(all_row.single_id_mapped_rows),
        "single_id_name_mapping_percentage": all_row.single_id_name_mapping_percentage,
        "single_id_row_mapping_percentage": all_row.single_id_row_mapping_percentage,
        "exact_pref_count": int(((mapping["mapping_status"] == "mapped") &
                                 (mapping["mapping_method"] == "exact_pref")).sum()),
        "exact_synonym_count": int(((mapping["mapping_status"] == "mapped") &
                                    (mapping["mapping_method"] == "exact_synonym")).sum()),
        "ambiguous_count": int((mapping["mapping_status"] == "unmapped_ambiguous").sum()),
        "unmatched_count": int((mapping["mapping_status"]=="unmapped_no_match").sum()),
        "mapped_combination_count": int((mapping["mapping_status"]=="mapped_combination_components").sum()),
        "partial_combination_count": int((mapping["mapping_status"]=="unmapped_combination_partial").sum()),
        "ambiguous_combination_count": int((mapping["mapping_status"]=="unmapped_combination_ambiguous").sum()),
        "all_invariants_passed": True, "input_files_modified": False, "temporary_files_remaining": 0,
        "absolute_paths_persisted": False, "created_at": datetime.now(timezone.utc).isoformat(),
        "code_commit": commit, "worktree_dirty": dirty}
    staging = {path:path.with_name(path.name+".part") for path in paths}
    try:
        mapping.to_parquet(staging[output],index=False); coverage.to_csv(staging[reports["coverage"]],index=False)
        status.to_csv(staging[reports["status"]],index=False)
        for frame,path in ((unmatched,reports["unmatched"]),(ambiguous,reports["ambiguous"]),(combinations,reports["combination"])):
            serial=frame.copy()
            for column in set(serial.columns)&{"years_seen","sample_raw_names","mapping_candidates","component_chembl_ids"}:
                serial[column]=serial[column].map(lambda values: json.dumps(list(values),default=lambda value:value.item()))
            serial.to_csv(staging[path],index=False)
        with staging[reports["summary"]].open("w",encoding="utf-8") as stream: json.dump(summary,stream,indent=2,sort_keys=True); stream.write("\n")
        for path in paths: os.replace(staging[path],path)
    finally:
        for path in staging.values():
            if path.exists(): path.unlink()
    return summary


def parse_args():
    parser=argparse.ArgumentParser(description="Build exact-only versioned FAERS-to-ChEMBL mapping and coverage reports.")
    parser.add_argument("--years",nargs="+",required=True); parser.add_argument("--vocabulary",type=Path,default=DEFAULT_VOCABULARY)
    parser.add_argument("--name-index",type=Path,default=DEFAULT_INDEX); parser.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument("--annual-input-root",type=Path,default=DEFAULT_INPUT_ROOT); parser.add_argument("--input-filename",default="resolved_drug_rows.parquet")
    parser.add_argument("--report-dir",type=Path,default=DEFAULT_REPORT_DIR); parser.add_argument("--normalization-version",default="norm_v1")
    parser.add_argument("--vocabulary-artifact-version",default="vocab_v1"); parser.add_argument("--mapping-artifact-version",default="exact_v1")
    parser.add_argument("--name-index-artifact-version",default="index_v1"); parser.add_argument("--batch-size",type=int,default=250_000)
    parser.add_argument("--overwrite",action="store_true"); return parser.parse_args()


def main():
    args=parse_args(); years=parse_years(args.years); vocabulary=pd.read_parquet(args.vocabulary)
    preferred,synonyms=load_name_index(args.name_index,chembl_version=36,normalization_version=args.normalization_version)
    mapping=build_exact_mapping(vocabulary,preferred,synonyms,chembl_version=36,normalization_version=args.normalization_version,
        vocabulary_artifact_version=args.vocabulary_artifact_version,mapping_artifact_version=args.mapping_artifact_version,
        name_index_artifact_version=args.name_index_artifact_version)
    annual=resolve_annual_drug_inputs(args.annual_input_root,years,args.input_filename)
    coverage,status=build_coverage(mapping,annual,batch_size=args.batch_size)
    summary=write_outputs(mapping,coverage,status,output=args.output,report_dir=args.report_dir,vocabulary_path=args.vocabulary,
        index_path=args.name_index,requested_years=years,normalization_version=args.normalization_version,
        mapping_artifact_version=args.mapping_artifact_version,overwrite=args.overwrite)
    print(json.dumps(summary,indent=2,sort_keys=True))


if __name__=="__main__": main()
