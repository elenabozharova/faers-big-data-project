"""Chunked, read-only schema and quality profiling for raw FAERS ZIP tables."""
from __future__ import annotations

import argparse, csv, hashlib, json, os, sys, time, zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import pandas as pd
from download_faers import ManifestError, Source, load_manifest, select_sources

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "config/faers_sources.csv"
DEFAULT_RAW_DIR = ROOT / "data/raw/faers"
DEFAULT_OUT_DIR = ROOT / "reports/data_quality"
TABLES = ("DEMO", "DRUG", "REAC")
EXPECTED = {
    "DEMO": ("primaryid", "caseid", "caseversion", "fda_dt", "init_fda_dt"),
    "DRUG": ("primaryid", "caseid", "drugname", "prod_ai", "role_cod", "drug_seq"),
    "REAC": ("primaryid", "caseid", "pt"),
}
COMMON = ("year", "quarter", "table", "zip_filename", "zip_path", "inner_filename",
          "zip_size_bytes", "inner_file_size_bytes", "row_count", "is_sampled",
          "max_rows_per_table", "column_count", "raw_column_names", "normalized_column_names",
          "schema_signature", "missing_expected_columns", "unexpected_columns", "missing_primaryid",
          "missing_caseid", "unique_primaryids", "unique_caseids", "profile_status", "warnings", "error")
METRICS = ("missing_caseversion", "invalid_caseversion", "missing_fda_dt", "invalid_fda_dt",
           "missing_init_fda_dt", "invalid_init_fda_dt", "duplicate_primaryid_rows",
           "duplicate_caseid_caseversion_rows", "missing_drugname", "missing_prod_ai",
           "missing_both_drugname_and_prod_ai", "missing_role_cod", "unique_raw_drug_names",
           "duplicate_primaryid_drug_seq_rows", "missing_pt", "unique_reaction_pts",
           "duplicate_primaryid_pt_rows")
PROFILE_FIELDS = COMMON + METRICS
LIST_FIELDS = ("raw_column_names", "normalized_column_names", "missing_expected_columns", "unexpected_columns", "warnings")

def normalized(name: Any) -> str: return str(name).strip().lower()
def signature(columns: list[str]) -> str:
    return hashlib.sha256(json.dumps(columns, separators=(",", ":")).encode()).hexdigest()
def missing(series: pd.Series) -> pd.Series:
    return series.isna() | series.fillna("").astype(str).str.strip().eq("")
def clean(series: pd.Series) -> pd.Series:
    return series[~missing(series)].astype(str).str.strip()

def repository_path(path: str | Path) -> str:
    """Return a portable, non-personal path for persisted reports.

    Repository paths are relative POSIX paths. For an external path, where a
    repository-relative value is impossible, only its filename is retained;
    runtime access continues to use the original Path object.
    """
    candidate = Path(path)
    try:
        return candidate.resolve().relative_to(ROOT.resolve()).as_posix()
    except (OSError, ValueError):
        return candidate.name.replace("\\", "/")

def locate_archive(raw_dir: Path, filename: str) -> Path:
    path = raw_dir / filename
    if path.exists(): return path
    if not raw_dir.exists(): raise FileNotFoundError(f"raw directory does not exist: {raw_dir}")
    found = [p for p in raw_dir.iterdir() if p.is_file() and p.name.casefold() == filename.casefold()]
    if len(found) != 1: raise FileNotFoundError(f"{'missing' if not found else 'ambiguous'} archive {filename!r}")
    return found[0]

def resolve_member(zf: zipfile.ZipFile, wanted: str) -> str:
    wanted = wanted.replace("\\", "/")
    exact = [n for n in zf.namelist() if n.replace("\\", "/") == wanted]
    found = exact or [n for n in zf.namelist() if n.replace("\\", "/").casefold() == wanted.casefold()]
    if len(found) != 1: raise KeyError(f"{'missing' if not found else 'ambiguous'} ZIP member {wanted!r}")
    return found[0]

def blank(source: Source, table: str, path: Path, limit: int | None) -> dict[str, Any]:
    row = dict.fromkeys(PROFILE_FIELDS)
    row.update(year=source.year, quarter=source.quarter, table=table, zip_filename=source.local_filename,
               zip_path=repository_path(path), row_count=0, is_sampled=limit is not None,
               max_rows_per_table=limit, profile_status="failed", warnings=[], error="")
    return row

def add_duplicate(state: dict, metric: str, *series: pd.Series | None) -> None:
    if any(s is None for s in series): return
    frame = pd.concat(series, axis=1)
    valid = ~frame.apply(missing).any(axis=1)
    for key in frame.loc[valid].astype(str).apply(lambda s: s.str.strip()).itertuples(index=False, name=None):
        if key in state[metric + "_seen"]: state[metric] += 1
        else: state[metric + "_seen"].add(key)

def profile_table(source: Source, table: str, path: Path, chunk_size: int, limit: int | None) -> dict[str, Any]:
    row = blank(source, table, path, limit)
    counts = {m: 0 for m in METRICS}
    sets = {"primaryids": set(), "caseids": set(), "drugnames": set(), "pts": set(),
            "duplicate_primaryid_rows_seen": set(), "duplicate_caseid_caseversion_rows_seen": set(),
            "duplicate_primaryid_drug_seq_rows_seen": set(), "duplicate_primaryid_pt_rows_seen": set()}
    state = counts | sets
    try:
        with zipfile.ZipFile(path) as zf:
            wanted = getattr(source, f"expected_inner_{table.lower()}")
            member = resolve_member(zf, wanted); info = zf.getinfo(member)
            with zf.open(member) as stream: raw = [str(c) for c in pd.read_csv(stream, sep="$", dtype=str, nrows=0).columns]
            cols = [normalized(c) for c in raw]
            if len(cols) != len(set(cols)): raise ValueError("columns collide after trim/lower normalization")
            remaining = limit
            with zf.open(member) as stream:
                for chunk in pd.read_csv(stream, sep="$", dtype=str, chunksize=chunk_size, low_memory=False):
                    chunk.columns = cols
                    if remaining is not None:
                        chunk = chunk.iloc[:remaining]; remaining -= len(chunk)
                    if chunk.empty: break
                    state["row_count"] = state.get("row_count", 0) + len(chunk)
                    get = lambda c: chunk[c] if c in chunk else None
                    pid, cid = get("primaryid"), get("caseid")
                    if pid is not None:
                        state["missing_primaryid"] = state.get("missing_primaryid", 0) + int(missing(pid).sum())
                        state["primaryids"].update(clean(pid))
                    if cid is not None:
                        state["missing_caseid"] = state.get("missing_caseid", 0) + int(missing(cid).sum())
                        state["caseids"].update(clean(cid))
                    if table == "DEMO":
                        ver = get("caseversion")
                        if pid is not None:
                            for value in clean(pid):
                                if value in state["duplicate_primaryid_rows_seen"]: state["duplicate_primaryid_rows"] += 1
                                else: state["duplicate_primaryid_rows_seen"].add(value)
                        if ver is not None:
                            miss = missing(ver); state["missing_caseversion"] += int(miss.sum())
                            valid = ver.fillna("").astype(str).str.strip().str.fullmatch(r"[+-]?\d+(?:\.0+)?")
                            state["invalid_caseversion"] += int((~miss & ~valid).sum())
                        add_duplicate(state, "duplicate_caseid_caseversion_rows", cid, ver)
                        for col in ("fda_dt", "init_fda_dt"):
                            values = get(col)
                            if values is not None:
                                miss = missing(values); state[f"missing_{col}"] += int(miss.sum())
                                parsed = pd.to_datetime(values.where(~miss), format="%Y%m%d", errors="coerce")
                                state[f"invalid_{col}"] += int((~miss & parsed.isna()).sum())
                    elif table == "DRUG":
                        name, ai, role = get("drugname"), get("prod_ai"), get("role_cod")
                        mn, ma = (missing(name) if name is not None else None), (missing(ai) if ai is not None else None)
                        if name is not None: state["missing_drugname"] += int(mn.sum()); state["drugnames"].update(clean(name))
                        if ai is not None: state["missing_prod_ai"] += int(ma.sum())
                        if mn is not None and ma is not None: state["missing_both_drugname_and_prod_ai"] += int((mn & ma).sum())
                        if role is not None: state["missing_role_cod"] += int(missing(role).sum())
                        add_duplicate(state, "duplicate_primaryid_drug_seq_rows", pid, get("drug_seq"))
                    else:
                        pt = get("pt")
                        if pt is not None: state["missing_pt"] += int(missing(pt).sum()); state["pts"].update(clean(pt))
                        add_duplicate(state, "duplicate_primaryid_pt_rows", pid, pt)
                    if remaining == 0: break
        absent = [c for c in EXPECTED[table] if c not in cols]
        warnings = [f"missing expected column: {c}" for c in absent]
        row.update(inner_filename=member, zip_size_bytes=path.stat().st_size, inner_file_size_bytes=info.file_size,
                   row_count=state.get("row_count", 0), column_count=len(raw), raw_column_names=raw,
                   normalized_column_names=cols, schema_signature=signature(cols), missing_expected_columns=absent,
                   unexpected_columns=[c for c in cols if c not in EXPECTED[table]],
                   missing_primaryid=state.get("missing_primaryid") if "primaryid" in cols else None,
                   missing_caseid=state.get("missing_caseid") if "caseid" in cols else None,
                   unique_primaryids=len(state["primaryids"]) if "primaryid" in cols else None,
                   unique_caseids=len(state["caseids"]) if "caseid" in cols else None,
                   profile_status="profiled", warnings=warnings)
        deps = {"missing_caseversion":("caseversion",), "invalid_caseversion":("caseversion",),
                "missing_fda_dt":("fda_dt",), "invalid_fda_dt":("fda_dt",),
                "missing_init_fda_dt":("init_fda_dt",), "invalid_init_fda_dt":("init_fda_dt",),
                "duplicate_primaryid_rows":("primaryid",), "duplicate_caseid_caseversion_rows":("caseid","caseversion"),
                "missing_drugname":("drugname",), "missing_prod_ai":("prod_ai",),
                "missing_both_drugname_and_prod_ai":("drugname","prod_ai"), "missing_role_cod":("role_cod",),
                "duplicate_primaryid_drug_seq_rows":("primaryid","drug_seq"), "missing_pt":("pt",),
                "duplicate_primaryid_pt_rows":("primaryid","pt")}
        applicable = {"DEMO": METRICS[:8], "DRUG": METRICS[8:14], "REAC": METRICS[14:]}[table]
        for metric in applicable:
            if metric == "unique_raw_drug_names": row[metric] = len(state["drugnames"]) if "drugname" in cols else None
            elif metric == "unique_reaction_pts": row[metric] = len(state["pts"]) if "pt" in cols else None
            elif all(c in cols for c in deps[metric]): row[metric] = state[metric]
        warning_metrics = ("missing_primaryid", "missing_caseid") + tuple(
            metric for metric in METRICS if metric.startswith(("missing_", "invalid_", "duplicate_"))
        )
        for metric in warning_metrics:
            if isinstance(row.get(metric), int) and row[metric] > 0: warnings.append(f"{metric}={row[metric]}")
        return row
    except (OSError, KeyError, ValueError, zipfile.BadZipFile, pd.errors.ParserError) as exc:
        row["error"] = str(exc); row["warnings"] = ["profile failed"]; return row

def schema_analysis(rows: list[dict]) -> tuple[dict, list, dict, dict, list]:
    variants, drift, first, forms, absent = {}, [], {}, {}, []
    good = [r for r in rows if r["profile_status"] == "profiled"]
    for table in TABLES:
        prior = None; seen = set(); grouped = {}; raw_forms = {}
        table_rows = sorted((r for r in good if r["table"] == table), key=lambda r:(r["year"],r["quarter"]))
        universe = sorted(set().union(*(set(r["normalized_column_names"]) for r in table_rows))) if table_rows else []
        for r in table_rows:
            sig=r["schema_signature"]; grouped.setdefault(sig,{"normalized_column_names":r["normalized_column_names"],"partitions":[]})["partitions"].append(f"{r['year']}Q{r['quarter']}")
            for raw,col in zip(r["raw_column_names"],r["normalized_column_names"]):
                raw_forms.setdefault(col,set()).add(raw)
                if col not in seen: first[f"{table}.{col}"]=f"{r['year']}Q{r['quarter']}"; seen.add(col)
            missing_here=sorted(set(universe)-set(r["normalized_column_names"]))
            if missing_here: absent.append({"year":r["year"],"quarter":r["quarter"],"table":table,"columns":missing_here})
            if prior and prior["schema_signature"] != sig:
                old,new=set(prior["normalized_column_names"]),set(r["normalized_column_names"])
                drift.append({"table":table,"previous_year_quarter":f"{prior['year']}Q{prior['quarter']}","current_year_quarter":f"{r['year']}Q{r['quarter']}","columns_added":sorted(new-old),"columns_removed":sorted(old-new),"only_casing_or_spacing_changed":old==new})
            prior=r
        variants[table]=grouped; forms[table]={k:sorted(v) for k,v in sorted(raw_forms.items()) if len(v)>1}
    return variants,drift,first,forms,absent

def build_summary(rows: list[dict], args, csv_path: Path, runtime: float) -> dict:
    good=[r for r in rows if r["profile_status"]=="profiled"]; bad=[r for r in rows if r["profile_status"]!="profiled"]
    variants,drift,first,forms,absent=schema_analysis(rows); complete=[r for r in good if not r["is_sampled"]]
    return {"generated_at_utc":datetime.now(timezone.utc).isoformat().replace("+00:00","Z"),"runtime_seconds":round(runtime,3),
      "requested_years":sorted(set(args.years)),"requested_quarters":sorted(set(args.quarters)),"requested_tables":list(args.tables),
      "is_sampled":args.max_rows_per_table is not None,"max_rows_per_table":args.max_rows_per_table,"chunk_size":args.chunk_size,
      "archives_requested":len({(r['year'],r['quarter']) for r in rows}),"archives_profiled":len({(r['year'],r['quarter']) for r in good}),
      "tables_requested":len(rows),"tables_profiled":len(good),"successful_tables":len(good),"failed_tables":len(bad),
      "total_rows_by_table":{t:sum(r['row_count'] for r in complete if r['table']==t) for t in args.tables},
      "total_rows_by_year":{str(y):sum(r['row_count'] for r in complete if r['year']==y) for y in sorted(set(args.years))},
      "total_rows_by_year_and_table":{str(y):{t:sum(r['row_count'] for r in complete if r['year']==y and r['table']==t) for t in args.tables} for y in sorted(set(args.years))},
      "schema_variant_count_by_table":{t:len(variants.get(t,{})) for t in args.tables},"schema_variants":variants,
      "schema_drift_events":drift,"columns_first_seen":first,"columns_absent_by_partition":absent,"normalized_column_raw_forms":forms,
      "missing_expected_columns_by_partition":[{"year":r['year'],"quarter":r['quarter'],"table":r['table'],"columns":r['missing_expected_columns']} for r in good if r['missing_expected_columns']],
      "quality_warnings":[{"year":r['year'],"quarter":r['quarter'],"table":r['table'],"warnings":r['warnings'],"error":r['error']} for r in rows if r['warnings']],
      "largest_tables":[{"year":r['year'],"quarter":r['quarter'],"table":r['table'],"row_count":r['row_count']} for r in sorted(good,key=lambda r:(-r['row_count'],r['year'],r['quarter'],r['table']))[:10]],
      "smallest_tables":[{"year":r['year'],"quarter":r['quarter'],"table":r['table'],"row_count":r['row_count']} for r in sorted(good,key=lambda r:(r['row_count'],r['year'],r['quarter'],r['table']))[:10]],
      "output_profile_csv":repository_path(csv_path),"source_manifest":repository_path(args.manifest),"raw_directory":repository_path(args.raw_dir)}

def write_reports(rows: list[dict], summary: dict, out: Path) -> None:
    out.mkdir(parents=True,exist_ok=True); cp=out/"raw_schema_profile.csv"; jp=out/"raw_schema_summary.json"; ct=out/"raw_schema_profile.csv.tmp"; jt=out/"raw_schema_summary.json.tmp"
    try:
        with ct.open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=PROFILE_FIELDS); w.writeheader()
            for row in rows:
                item=dict(row)
                for field in LIST_FIELDS: item[field]="" if item[field] is None else json.dumps(item[field],ensure_ascii=False,separators=(",",":"))
                w.writerow(item)
        with jt.open("w",encoding="utf-8") as f: json.dump(summary,f,indent=2,sort_keys=True,ensure_ascii=False); f.write("\n")
        os.replace(ct,cp); os.replace(jt,jp)
    finally:
        for p in (ct,jt):
            if p.exists(): p.unlink()

def parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--years",type=int,nargs="+",required=True)
    p.add_argument("--quarters",type=int,choices=(1,2,3,4),nargs="+",default=(1,2,3,4)); p.add_argument("--manifest",type=Path,default=DEFAULT_MANIFEST)
    p.add_argument("--raw-dir",type=Path,default=DEFAULT_RAW_DIR); p.add_argument("--out-dir",type=Path,default=DEFAULT_OUT_DIR)
    p.add_argument("--tables",type=str.upper,choices=TABLES,nargs="+",default=TABLES); p.add_argument("--chunk-size",type=int,default=250000)
    p.add_argument("--max-rows-per-table",type=int); p.add_argument("--fail-fast",action="store_true"); return p

def main(argv=None) -> int:
    p=parser(); a=p.parse_args(argv)
    if any(y<=0 for y in a.years): p.error("years must be positive")
    if a.chunk_size<=0: p.error("--chunk-size must be positive")
    if a.max_rows_per_table is not None and a.max_rows_per_table<=0: p.error("--max-rows-per-table must be positive")
    if a.max_rows_per_table is not None and a.out_dir.resolve()==DEFAULT_OUT_DIR.resolve(): a.out_dir=DEFAULT_OUT_DIR/"sample"
    a.tables=tuple(dict.fromkeys(a.tables))
    try: sources=select_sources(load_manifest(a.manifest),a.years,a.quarters)
    except ManifestError as exc: print(f"ERROR: {exc}",file=sys.stderr); return 2
    started=time.perf_counter(); rows=[]
    for source in sources:
        try: path=locate_archive(a.raw_dir,source.local_filename)
        except OSError as exc:
            for table in a.tables:
                row=blank(source,table,a.raw_dir/source.local_filename,a.max_rows_per_table); row["error"]=str(exc); row["warnings"]=["profile failed"]; rows.append(row); print(f"{source.year} Q{source.quarter} {table}  FAILED {exc}")
                if a.fail_fast:return 1
            continue
        for table in a.tables:
            row=profile_table(source,table,path,a.chunk_size,a.max_rows_per_table); rows.append(row)
            print(f"{source.year} Q{source.quarter} {table}  {'PROFILED rows='+str(row['row_count']) if row['profile_status']=='profiled' else 'FAILED '+row['error']}")
            if a.fail_fast and row["profile_status"]!="profiled": return 1
    rows.sort(key=lambda r:(r["year"],r["quarter"],TABLES.index(r["table"])))
    summary=build_summary(rows,a,a.out_dir/"raw_schema_profile.csv",time.perf_counter()-started); write_reports(rows,summary,a.out_dir)
    print(f"requested={len(rows)} profiled={summary['successful_tables']} failed={summary['failed_tables']} runtime_seconds={summary['runtime_seconds']}")
    return 1 if summary["failed_tables"] else 0
if __name__=="__main__": raise SystemExit(main())
