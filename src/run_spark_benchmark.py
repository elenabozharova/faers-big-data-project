"""MS8C: run the Spark scaling benchmark unattended.

Phases (each config runs in a fresh Spark session via spark_pair_reactions.py):
  A  scaling grid: all years at every core count       -> reports/spark/benchmark_scaling.csv
  B  repeats: one year, N repeats at every core count  -> reports/spark/benchmark_repeats.csv
  C  compute only (--noop), same configs as B          -> reports/spark/benchmark_noop.csv

Repeats (phase A: --a-repeats) alternate the core-count order (1..24, then 24..1) so slow drift such
as thermal throttling does not bias one end of the curve. A failed config is
recorded and the run continues. The machine is kept awake while it runs.
Benchmark Parquet output goes to %TEMP%/faers_spark_bench (not indexed by
Windows Search). Per-config logs go to data/interim/spark/bench_logs/ (gitignored); one manifest
row per config (exit code, wall time, free memory before start) goes to
reports/spark/benchmark_runs.csv.
"""

from __future__ import annotations

import argparse
import ctypes
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

from spark_pair_reactions import append_csv

ROOT = Path(__file__).resolve().parents[1]
JOB = ROOT / "src" / "spark_pair_reactions.py"
REPORTS = ROOT / "reports" / "spark"
LOGS = ROOT / "data" / "interim" / "spark" / "bench_logs"
# Outside the repo: Windows Search does not index %TEMP%, so indexing new Parquet files cannot interfere.
BENCH_OUT = Path(tempfile.gettempdir()) / "faers_spark_bench"


def keep_awake(on: bool) -> None:
    """Windows: block sleep while this process runs (no system setting is changed)."""
    if sys.platform == "win32":
        es_continuous, es_system_required = 0x80000000, 0x00000001
        ctypes.windll.kernel32.SetThreadExecutionState(es_continuous | (es_system_required if on else 0))


def free_memory_gb() -> float | None:
    if sys.platform != "win32":
        return None

    class MemoryStatus(ctypes.Structure):
        _fields_ = [("length", ctypes.c_ulong), ("load", ctypes.c_ulong), ("total_phys", ctypes.c_ulonglong),
                    ("avail_phys", ctypes.c_ulonglong), ("total_page", ctypes.c_ulonglong),
                    ("avail_page", ctypes.c_ulonglong), ("total_virtual", ctypes.c_ulonglong),
                    ("avail_virtual", ctypes.c_ulonglong), ("avail_extended", ctypes.c_ulonglong)]

    status = MemoryStatus(length=ctypes.sizeof(MemoryStatus))
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
    return round(status.avail_phys / 2**30, 1)


def plan(phases: list[str], cores: list[int], years: str, repeat_year: int, repeats: int,
         a_repeats: int = 1) -> list[dict]:
    configs = []
    if "A" in phases:
        for repeat in range(a_repeats):
            order = cores if repeat % 2 == 0 else cores[::-1]
            configs += [{"phase": "A", "cores": c, "years": years, "noop": False, "csv": "benchmark_scaling.csv",
                         "repeat": repeat + 1} for c in order]
    for phase, noop, csv_name in (("B", False, "benchmark_repeats.csv"), ("C", True, "benchmark_noop.csv")):
        if phase in phases:
            for repeat in range(repeats):
                order = cores if repeat % 2 == 0 else cores[::-1]
                configs += [{"phase": phase, "cores": c, "years": str(repeat_year), "noop": noop,
                             "csv": csv_name, "repeat": repeat + 1} for c in order]
    return configs


def run(config: dict, index: int, total: int, report_dir: Path) -> dict:
    name = f"{config['phase']}_c{config['cores']:02d}_r{config.get('repeat', 1)}_{config['years']}"
    command = [sys.executable, str(JOB), "--years", config["years"], "--cores", str(config["cores"]),
               "--out-dir", str(BENCH_OUT), "--benchmark-csv", str(report_dir / config["csv"])]
    if config["noop"]:
        command.append("--noop")
    row = {"config": name, "phase": config["phase"], "cores": config["cores"], "years": config["years"],
           "repeat": config.get("repeat", 1), "free_memory_gb_before": free_memory_gb(),
           "started_at": datetime.now(timezone.utc).isoformat()}
    print(f"[{index}/{total}] {name} (free memory {row['free_memory_gb_before']} GB)", flush=True)
    started = time.perf_counter()
    with (LOGS / f"{name}.log").open("w", encoding="utf-8") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)
    row.update(exit_code=result.returncode, wall_seconds=round(time.perf_counter() - started, 1))
    print(f"    exit {result.returncode} in {row['wall_seconds']} s", flush=True)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--phases", nargs="+", default=["A", "B", "C"], choices=["A", "B", "C"])
    parser.add_argument("--cores", nargs="+", type=int, default=[1, 2, 4, 8, 16, 24])
    parser.add_argument("--years", default="2013-2025", help="phase A year range")
    parser.add_argument("--repeat-year", type=int, default=2023, help="phases B and C")
    parser.add_argument("--repeats", type=int, default=3, help="phases B and C")
    parser.add_argument("--a-repeats", type=int, default=1, help="phase A; repeats alternate core order")
    parser.add_argument("--report-dir", type=Path, default=REPORTS, help="override for smoke tests")
    parser.add_argument("--dry-run", action="store_true", help="print the plan and exit")
    args = parser.parse_args()

    configs = plan(args.phases, args.cores, args.years, args.repeat_year, args.repeats, args.a_repeats)
    if args.dry_run:
        for config in configs:
            print(config)
        return
    if shutil.which("java") is None:
        raise SystemExit("java not found on PATH: install Java 17+ and open a new terminal")

    args.report_dir.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    keep_awake(True)
    failed = []
    try:
        for index, config in enumerate(configs, 1):
            row = run(config, index, len(configs), args.report_dir)
            append_csv(args.report_dir / "benchmark_runs.csv", row)
            if row["exit_code"]:
                failed.append(row["config"])
    finally:
        keep_awake(False)
    print("done; failed configs: " + (", ".join(failed) if failed else "none"))


if __name__ == "__main__":
    main()
