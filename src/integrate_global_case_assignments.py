"""MS4B boundary: apply MS4A selections to raw FAERS DRUG and REAC rows."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ASSIGNMENTS = ROOT / "data" / "interim" / "case_resolution" / "global_case_assignments.parquet"
DEFAULT_REPORT_DIR = ROOT / "reports" / "data_quality"
REQUIRED_COLUMNS = {
    "caseid", "selected_primaryid", "assigned_report_year", "assignment_status",
    "selected_source_year", "selected_source_quarter",
}


def normalize_identifier(values: pd.Series) -> pd.Series:
    """Trim identifiers as strings without numeric coercion (leading zeros survive)."""
    result = values.astype("string").str.strip()
    return result.mask(result == "", pd.NA)


def repository_path(path: Path | str) -> str:
    """Serialize repository paths as relative POSIX paths; external paths as names only."""
    candidate = Path(path).resolve()
    try:
        return candidate.relative_to(ROOT).as_posix()
    except ValueError:
        return candidate.name


def load_global_assignments(path: Path | str) -> pd.DataFrame:
    assignments = pd.read_parquet(path)
    missing = sorted(REQUIRED_COLUMNS - set(assignments.columns))
    if missing:
        raise ValueError(f"global assignments missing required columns: {missing}")
    result = assignments.copy()
    for column in ("caseid", "selected_primaryid"):
        result[column] = normalize_identifier(result[column])
    included = result["assignment_status"].eq("included")
    if result.loc[included, ["caseid", "selected_primaryid"]].isna().any().any():
        raise ValueError("included assignments contain blank caseid or selected_primaryid")
    selected = result.loc[included]
    if selected["caseid"].duplicated().any():
        raise ValueError("an included CASEID is assigned more than once")
    if selected.groupby("selected_primaryid")["caseid"].nunique().gt(1).any():
        raise ValueError("a selected_primaryid maps to multiple CASEIDs")
    if selected.duplicated(["caseid", "selected_primaryid"]).any():
        raise ValueError("duplicate included caseid/selected_primaryid pair")
    if selected.groupby("caseid")["assigned_report_year"].nunique().gt(1).any():
        raise ValueError("an included CASEID maps to multiple analytical years")
    return result


def selected_cases_for_year(assignments: pd.DataFrame, year: int) -> pd.DataFrame:
    """Select only the four approved inclusion fields; audit flags are ignored."""
    mask = assignments["assignment_status"].eq("included") & assignments["assigned_report_year"].eq(year)
    optional = [c for c in ("selected_caseversion", "assignment_method", "later_followup_primaryids") if c in assignments]
    columns = [
        "caseid", "selected_primaryid", "assigned_report_year",
        "selected_source_year", "selected_source_quarter", *optional,
    ]
    selected = assignments.loc[mask, columns].copy().sort_values("caseid", kind="stable").reset_index(drop=True)
    if selected[["caseid", "selected_primaryid", "selected_source_year", "selected_source_quarter"]].isna().any().any():
        raise ValueError(f"year {year}: included assignment lacks identifier or source provenance")
    return selected


def adapt_retained_cases(selected: pd.DataFrame) -> pd.DataFrame:
    rename = {"selected_primaryid": "retained_primaryid", "selected_caseversion": "caseversion"}
    return selected.rename(columns=rename).copy()


def filter_selected_rows(raw: pd.DataFrame, selected: pd.DataFrame, table_name: str) -> tuple[pd.DataFrame, dict]:
    """Many-to-one PRIMARYID join with assignment CASEID authoritative and mismatch-fatal."""
    if "primaryid" not in raw:
        raise ValueError(f"{table_name} missing primaryid")
    working = raw.copy()
    working["primaryid"] = normalize_identifier(working["primaryid"])
    mapping = selected[["selected_primaryid", "caseid"]].rename(
        columns={"selected_primaryid": "primaryid", "caseid": "assignment_caseid"}
    )
    filtered = working.merge(mapping, on="primaryid", how="inner", validate="many_to_one")
    blank = mismatch = 0
    if "caseid" in filtered:
        raw_caseid = normalize_identifier(filtered["caseid"])
        blank_mask = raw_caseid.isna()
        mismatch_mask = ~blank_mask & raw_caseid.ne(filtered["assignment_caseid"])
        blank = int(blank_mask.sum())
        mismatch = int(mismatch_mask.sum())
        if mismatch:
            raise ValueError(f"{table_name}: {mismatch:,} nonblank raw CASEIDs conflict with global assignments")
        filtered = filtered.drop(columns="caseid")
    filtered = filtered.rename(columns={"assignment_caseid": "caseid"})
    return filtered, {
        "raw_rows_scanned": int(len(raw)), "raw_rows_retained": int(len(filtered)),
        "raw_caseid_blank_count": blank, "raw_caseid_mismatch_count": mismatch,
        "matched_primaryids": set(filtered["primaryid"].dropna().unique()),
    }


def integrate_year(
    assignments: pd.DataFrame,
    year: int,
    partition_loader: Callable[[int, int, str], pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """Route selected IDs to exact physical partitions and return the MS4B boundary."""
    selected = selected_cases_for_year(assignments, year)
    drug_parts, reac_parts = [], []
    totals = {"drug_scanned": 0, "reac_scanned": 0, "blank": 0, "mismatch": 0}
    matched = {"drug": set(), "reac": set()}
    pid_partitions: dict[str, set[tuple[int, int]]] = {}
    for (source_year, source_quarter), group in selected.groupby(
        ["selected_source_year", "selected_source_quarter"], sort=True
    ):
        partition = (int(source_year), int(source_quarter))
        for pid in group["selected_primaryid"]:
            pid_partitions.setdefault(pid, set()).add(partition)
        for table_name, destination in (("drug", drug_parts), ("reac", reac_parts)):
            raw = partition_loader(*partition, table_name)
            filtered, stats = filter_selected_rows(raw, group, table_name.upper())
            destination.append(filtered)
            totals[f"{table_name}_scanned"] += stats["raw_rows_scanned"]
            totals["blank"] += stats["raw_caseid_blank_count"]
            totals["mismatch"] += stats["raw_caseid_mismatch_count"]
            matched[table_name].update(stats["matched_primaryids"])

    def concatenate(parts: list[pd.DataFrame]) -> pd.DataFrame:
        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["primaryid", "caseid"])

    drug, reac = concatenate(drug_parts), concatenate(reac_parts)
    selected_pids = set(selected["selected_primaryid"])
    retained_caseids = set(selected["caseid"])
    for name, frame in (("DRUG", drug), ("REAC", reac)):
        if not set(frame["primaryid"]).issubset(selected_pids):
            raise AssertionError(f"{name} contains an unexpected PRIMARYID")
        if not set(frame["caseid"]).issubset(retained_caseids):
            raise AssertionError(f"{name} contains an unexpected CASEID")
    drug_cases, reac_cases = set(drug["caseid"]), set(reac["caseid"])
    later = set()
    if "later_followup_primaryids" in selected:
        for values in selected["later_followup_primaryids"].dropna():
            # PyArrow list columns are commonly materialized as numpy arrays,
            # while synthetic/unit data often uses Python lists.
            if not isinstance(values, (str, bytes)):
                later.update(
                    normalized for value in values
                    if (normalized := str(value).strip())
                )
    leaked_later = later & (set(drug["primaryid"]) | set(reac["primaryid"]))
    if leaked_later:
        raise AssertionError(f"later follow-up PRIMARYIDs reached downstream: {sorted(leaked_later)[:5]}")
    report = {
        "analytical_year": int(year),
        "included_assignments_expected": int(len(selected)),
        "retained_cases_written": int(len(selected)),
        "unique_selected_primaryids": int(selected["selected_primaryid"].nunique()),
        "selected_primaryids_matched_in_drug": len(matched["drug"]),
        "selected_primaryids_matched_in_reac": len(matched["reac"]),
        "selected_cases_with_drug": len(drug_cases), "selected_cases_with_reac": len(reac_cases),
        "selected_cases_with_both": len(drug_cases & reac_cases),
        "selected_cases_missing_drug": len(retained_caseids - drug_cases),
        "selected_cases_missing_reac": len(retained_caseids - reac_cases),
        "selected_cases_missing_both": len(retained_caseids - (drug_cases | reac_cases)),
        "raw_drug_rows_retained": int(len(drug)), "raw_reac_rows_retained": int(len(reac)),
        "raw_drug_rows_scanned": totals["drug_scanned"], "raw_reac_rows_scanned": totals["reac_scanned"],
        "raw_caseid_blank_count": totals["blank"], "raw_caseid_mismatch_count": totals["mismatch"],
        "selected_primaryid_multiple_source_match_count": sum(len(v) > 1 for v in pid_partitions.values()),
        "unexpected_primaryid_count": 0,
        "later_followup_primaryids_excluded": len(later),
    }
    return adapt_retained_cases(selected), drug, reac, report


def write_integration_report(report: dict, report_dir: Path | str = DEFAULT_REPORT_DIR) -> Path:
    directory = Path(report_dir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"ms4b_global_assignment_integration_{report['analytical_year']}.json"
    payload = dict(report)
    payload["report_path"] = repository_path(path)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return path
