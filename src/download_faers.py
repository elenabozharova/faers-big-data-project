"""Reproducibly acquire and validate quarterly FAERS ASCII archives."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys
import urllib.error
import urllib.request
import zipfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = PROJECT_ROOT / "config" / "faers_sources.csv"
DEFAULT_OUT_DIR = PROJECT_ROOT / "data" / "raw" / "faers"
DEFAULT_LOG_PATH = PROJECT_ROOT / "data" / "interim" / "acquisition" / "faers_download_log.csv"

MANIFEST_COLUMNS = (
    "year", "quarter", "source_system", "format", "url", "source_filename",
    "local_filename", "expected_inner_demo", "expected_inner_drug",
    "expected_inner_reac", "notes",
)
LOG_COLUMNS = (
    "attempted_at_utc", "year", "quarter", "source_url", "source_filename",
    "local_path", "status", "file_size_bytes", "sha256", "previous_sha256",
    "changed_from_previous", "http_etag", "http_last_modified", "error",
)
SUCCESS_STATUSES = {"downloaded", "overwritten", "verified_existing", "skipped_valid"}


class ManifestError(ValueError):
    """The source manifest is invalid."""


class ArchiveError(ValueError):
    def __init__(self, status: str, message: str):
        super().__init__(message)
        self.status = status


@dataclass(frozen=True)
class Source:
    year: int
    quarter: int
    source_system: str
    format: str
    url: str
    source_filename: str
    local_filename: str
    expected_inner_demo: str
    expected_inner_drug: str
    expected_inner_reac: str
    notes: str = ""


@dataclass
class Event:
    attempted_at_utc: str
    year: int
    quarter: int
    source_url: str
    source_filename: str
    local_path: str
    status: str
    file_size_bytes: str = ""
    sha256: str = ""
    previous_sha256: str = ""
    changed_from_previous: str = ""
    http_etag: str = ""
    http_last_modified: str = ""
    error: str = ""


def _safe_basename(value: str, field: str) -> str:
    value = value.strip()
    path = PurePosixPath(value.replace("\\", "/"))
    if not value or path.is_absolute() or len(path.parts) != 1 or value in {".", ".."}:
        raise ManifestError(f"{field} must be a safe basename: {value!r}")
    return value


def load_manifest(path: Path | str) -> list[Source]:
    path = Path(path)
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            missing = [name for name in MANIFEST_COLUMNS if name not in (reader.fieldnames or [])]
            if missing:
                raise ManifestError(f"manifest missing required column(s): {', '.join(missing)}")
            rows = list(reader)
    except OSError as exc:
        raise ManifestError(f"cannot read manifest {path}: {exc}") from exc

    sources: list[Source] = []
    keys: set[tuple[int, int]] = set()
    filenames: set[str] = set()
    for line, row in enumerate(rows, start=2):
        try:
            year = int(row["year"].strip())
            quarter = int(row["quarter"].strip())
        except (AttributeError, ValueError) as exc:
            raise ManifestError(f"line {line}: year and quarter must be integers") from exc
        if not 1900 <= year <= 2100:
            raise ManifestError(f"line {line}: year is not sensible: {year}")
        if quarter not in {1, 2, 3, 4}:
            raise ManifestError(f"line {line}: quarter must be one of 1, 2, 3, 4")
        key = (year, quarter)
        if key in keys:
            raise ManifestError(f"line {line}: duplicate year-quarter {year} Q{quarter}")
        keys.add(key)
        local = _safe_basename(row["local_filename"], "local_filename")
        source_name = _safe_basename(row["source_filename"], "source_filename")
        if local.casefold() in filenames:
            raise ManifestError(f"line {line}: duplicate local_filename {local!r}")
        filenames.add(local.casefold())
        url = row["url"].strip()
        if not url.lower().startswith("https://"):
            raise ManifestError(f"line {line}: url must be a non-empty HTTPS URL")
        members = [row[name].strip() for name in
                   ("expected_inner_demo", "expected_inner_drug", "expected_inner_reac")]
        if not all(members):
            raise ManifestError(f"line {line}: expected DEMO, DRUG, and REAC members are required")
        sources.append(Source(
            year, quarter, row["source_system"].strip(), row["format"].strip(), url,
            source_name, local, *members, row.get("notes", "").strip(),
        ))
    return sorted(sources, key=lambda item: (item.year, item.quarter))


def select_sources(sources: Iterable[Source], years: Iterable[int], quarters: Iterable[int]) -> list[Source]:
    requested = {(year, quarter) for year in years for quarter in quarters}
    selected = sorted((s for s in sources if (s.year, s.quarter) in requested),
                      key=lambda item: (item.year, item.quarter))
    found = {(s.year, s.quarter) for s in selected}
    missing = sorted(requested - found)
    if missing:
        labels = ", ".join(f"{year} Q{quarter}" for year, quarter in missing)
        raise ManifestError(f"requested quarter(s) absent from manifest: {labels}")
    return selected


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_archive(path: Path, source: Source) -> dict[str, str]:
    try:
        with zipfile.ZipFile(path) as archive:
            bad_member = archive.testzip()
            if bad_member is not None:
                raise ArchiveError("failed_invalid_zip", f"ZIP CRC check failed for {bad_member!r}")
            names = archive.namelist()
    except (zipfile.BadZipFile, OSError) as exc:
        raise ArchiveError("failed_invalid_zip", f"invalid ZIP: {exc}") from exc

    actual: dict[str, str] = {}
    for role, expected in (
        ("demo", source.expected_inner_demo), ("drug", source.expected_inner_drug),
        ("reac", source.expected_inner_reac),
    ):
        normalized = expected.replace("\\", "/")
        exact = [name for name in names if name.replace("\\", "/") == normalized]
        matches = exact or [name for name in names if name.replace("\\", "/").casefold() == normalized.casefold()]
        if len(matches) != 1:
            detail = "missing" if not matches else "ambiguous case-insensitive matches"
            raise ArchiveError("failed_missing_members", f"{detail} for {expected!r}")
        actual[role] = matches[0]
    return actual


def _latest_checksums(log_path: Path) -> dict[tuple[int, int], str]:
    checksums: dict[tuple[int, int], str] = {}
    if not log_path.exists():
        return checksums
    try:
        with log_path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("status") in SUCCESS_STATUSES and row.get("sha256"):
                    checksums[(int(row["year"]), int(row["quarter"]))] = row["sha256"]
    except (OSError, ValueError, csv.Error) as exc:
        raise ManifestError(f"cannot read runtime log {log_path}: {exc}") from exc
    return checksums


def _append_event(log_path: Path, event: Event) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not log_path.exists() or log_path.stat().st_size == 0
    with log_path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOG_COLUMNS)
        if new_file:
            writer.writeheader()
        writer.writerow(asdict(event))


def _event(source: Source, final: Path, status: str, **values: object) -> Event:
    return Event(
        attempted_at_utc=datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        year=source.year, quarter=source.quarter, source_url=source.url,
        source_filename=source.source_filename, local_path=str(final), status=status,
        **values,
    )


def _download(source: Source, part: Path) -> tuple[str, str, str, int]:
    digest = hashlib.sha256()
    request = urllib.request.Request(source.url, headers={"User-Agent": "faers-big-data-project/1.0"})
    with urllib.request.urlopen(request) as response, part.open("wb") as output:
        etag = response.headers.get("ETag", "")
        modified = response.headers.get("Last-Modified", "")
        size = 0
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            output.write(chunk)
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), etag, modified, size


def process_source(source: Source, out_dir: Path, log_path: Path, previous: str,
                   verify_only: bool, overwrite: bool) -> Event:
    final = out_dir / source.local_filename
    part = final.with_name(final.name + ".part")
    replacing_existing = final.exists() and overwrite
    base = {"previous_sha256": previous}
    if final.exists() and not overwrite:
        try:
            validate_archive(final, source)
            checksum = sha256_file(final)
            changed = bool(previous and checksum != previous)
            if changed:
                return _event(source, final, "failed_checksum_changed", file_size_bytes=str(final.stat().st_size),
                              sha256=checksum, changed_from_previous="true", **base,
                              error="existing archive checksum differs from latest successful log entry")
            status = "verified_existing" if verify_only or not previous else "skipped_valid"
            return _event(source, final, status, file_size_bytes=str(final.stat().st_size),
                          sha256=checksum, changed_from_previous="false" if previous else "", **base)
        except ArchiveError as exc:
            return _event(source, final, exc.status, file_size_bytes=str(final.stat().st_size),
                          **base, error=str(exc))
    if verify_only:
        return _event(source, final, "failed_missing", **base, error="archive does not exist")

    out_dir.mkdir(parents=True, exist_ok=True)
    if part.exists():
        part.unlink()
    etag = modified = checksum = ""
    size: int | str = ""
    try:
        checksum, etag, modified, size = _download(source, part)
        validate_archive(part, source)
        os.replace(part, final)
        status = "overwritten" if replacing_existing else "downloaded"
        return _event(source, final, status, file_size_bytes=str(size), sha256=checksum,
                      changed_from_previous=str(bool(previous and checksum != previous)).lower() if previous else "",
                      http_etag=etag, http_last_modified=modified, **base)
    except ArchiveError as exc:
        return _event(source, final, exc.status, file_size_bytes=str(size), sha256=checksum,
                      http_etag=etag, http_last_modified=modified, **base, error=str(exc))
    except (OSError, urllib.error.URLError, urllib.error.HTTPError) as exc:
        return _event(source, final, "failed_download", file_size_bytes=str(size), sha256=checksum,
                      http_etag=etag, http_last_modified=modified, **base, error=str(exc))
    finally:
        if part.exists():
            part.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--years", type=int, nargs="+", required=True)
    parser.add_argument("--quarters", type=int, choices=(1, 2, 3, 4), nargs="+", default=(1, 2, 3, 4))
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--log-path", type=Path, default=DEFAULT_LOG_PATH)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.verify_only and args.overwrite:
        parser.error("--verify-only and --overwrite cannot be combined")
    try:
        selected = select_sources(load_manifest(args.manifest), args.years, args.quarters)
        previous = {} if args.dry_run else _latest_checksums(args.log_path)
    except ManifestError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if args.dry_run:
        for source in selected:
            print(f"{source.year} Q{source.quarter}  DRY_RUN  {args.out_dir / source.local_filename}")
        print(f"requested={len(selected)} downloaded=0 overwritten=0 skipped_valid=0 verified_existing=0 failed=0")
        return 0

    events = []
    for source in selected:
        event = process_source(source, args.out_dir, args.log_path,
                               previous.get((source.year, source.quarter), ""),
                               args.verify_only, args.overwrite)
        _append_event(args.log_path, event)
        events.append(event)
        print(f"{source.year} Q{source.quarter}  {event.status.upper()}")

    counts = {name: sum(event.status == name for event in events) for name in
              ("downloaded", "overwritten", "skipped_valid", "verified_existing")}
    failed = sum(event.status not in SUCCESS_STATUSES for event in events)
    print(f"requested={len(events)} downloaded={counts['downloaded']} overwritten={counts['overwritten']} "
          f"skipped_valid={counts['skipped_valid']} verified_existing={counts['verified_existing']} failed={failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
