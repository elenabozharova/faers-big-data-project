"""Deterministic tests for FAERS archive acquisition."""

import csv
import io
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import download_faers as dl  # noqa: E402


HEADERS = list(dl.MANIFEST_COLUMNS)


def zip_bytes(missing=""):
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as archive:
        for role in ("DEMO", "DRUG", "REAC"):
            if role != missing:
                archive.writestr(f"ASCII/{role}24Q1.txt", "header\nvalue\n")
    return output.getvalue()


class Response(io.BytesIO):
    def __init__(self, body):
        super().__init__(body)
        self.headers = {"ETag": '"test"', "Last-Modified": "today"}

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.close()


class DownloadFaersTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.manifest = self.root / "manifest.csv"
        self.out = self.root / "raw"
        self.log = self.root / "interim" / "log.csv"

    def tearDown(self):
        self.temp.cleanup()

    def row(self, year=2024, quarter=1, filename="faers_ascii_2024q1.zip"):
        yy = str(year)[2:]
        return {
            "year": year, "quarter": quarter, "source_system": "FAERS", "format": "ASCII",
            "url": f"https://example.test/{filename}", "source_filename": filename,
            "local_filename": filename, "expected_inner_demo": f"ASCII/DEMO{yy}Q{quarter}.txt",
            "expected_inner_drug": f"ASCII/DRUG{yy}Q{quarter}.txt",
            "expected_inner_reac": f"ASCII/REAC{yy}Q{quarter}.txt", "notes": "test",
        }

    def write_manifest(self, rows, headers=HEADERS):
        with self.manifest.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers)
            writer.writeheader()
            writer.writerows(rows)

    def args(self, *extra):
        return ["--years", "2024", "--quarters", "1", "--manifest", str(self.manifest),
                "--out-dir", str(self.out), "--log-path", str(self.log), *extra]

    def test_valid_manifest_and_extra_column(self):
        row = self.row()
        row["future_metadata"] = "ok"
        self.write_manifest([row], HEADERS + ["future_metadata"])
        source = dl.load_manifest(self.manifest)[0]
        self.assertEqual((source.year, source.quarter), (2024, 1))

    def test_missing_required_column(self):
        headers = [h for h in HEADERS if h != "url"]
        self.write_manifest([{k: v for k, v in self.row().items() if k in headers}], headers)
        with self.assertRaisesRegex(dl.ManifestError, "missing required"):
            dl.load_manifest(self.manifest)

    def test_duplicate_year_quarter_and_filename_rejected(self):
        row = self.row()
        self.write_manifest([row, row])
        with self.assertRaisesRegex(dl.ManifestError, "duplicate year-quarter"):
            dl.load_manifest(self.manifest)
        self.write_manifest([row, self.row(2024, 2)])
        with self.assertRaisesRegex(dl.ManifestError, "duplicate local_filename"):
            dl.load_manifest(self.manifest)

    def test_unsafe_filename_rejected(self):
        for unsafe in ("../x.zip", "/x.zip", "folder/x.zip", "folder\\x.zip"):
            row = self.row(filename=unsafe)
            self.write_manifest([row])
            with self.assertRaisesRegex(dl.ManifestError, "safe basename"):
                dl.load_manifest(self.manifest)

    def test_selection_is_ordered_and_missing_fails(self):
        sources = [dl.Source(**self.row(2024, q, f"q{q}.zip")) for q in (3, 1, 2)]
        self.assertEqual([s.quarter for s in dl.select_sources(sources, [2024], [3, 1])], [1, 3])
        with self.assertRaisesRegex(dl.ManifestError, "absent"):
            dl.select_sources(sources, [2024], [4])

    def test_existing_zip_baseline_then_skip(self):
        self.write_manifest([self.row()])
        self.out.mkdir()
        archive = self.out / self.row()["local_filename"]
        archive.write_bytes(zip_bytes())
        with patch("download_faers.urllib.request.urlopen") as request:
            self.assertEqual(dl.main(self.args()), 0)
            self.assertEqual(dl.main(self.args()), 0)
            request.assert_not_called()
        with self.log.open(encoding="utf-8") as handle:
            statuses = [r["status"] for r in csv.DictReader(handle)]
        self.assertEqual(statuses, ["verified_existing", "skipped_valid"])

    def test_checksum_change_fails_without_overwrite(self):
        self.write_manifest([self.row()])
        self.out.mkdir()
        archive = self.out / self.row()["local_filename"]
        archive.write_bytes(zip_bytes())
        self.assertEqual(dl.main(self.args()), 0)
        with zipfile.ZipFile(archive, "a") as target:
            target.writestr("extra.txt", "changed")
        self.assertEqual(dl.main(self.args()), 1)
        with self.log.open(encoding="utf-8") as handle:
            self.assertEqual(list(csv.DictReader(handle))[-1]["status"], "failed_checksum_changed")

    def test_invalid_and_missing_member_failures(self):
        self.write_manifest([self.row()])
        self.out.mkdir()
        archive = self.out / self.row()["local_filename"]
        archive.write_bytes(b"not zip")
        self.assertEqual(dl.main(self.args("--verify-only")), 1)
        archive.write_bytes(zip_bytes("REAC"))
        self.assertEqual(dl.main(self.args("--verify-only")), 1)
        with self.log.open(encoding="utf-8") as handle:
            statuses = [r["status"] for r in csv.DictReader(handle)]
        self.assertEqual(statuses, ["failed_invalid_zip", "failed_missing_members"])

    def test_crc_failure_is_rejected(self):
        self.write_manifest([self.row()])
        self.out.mkdir()
        archive = self.out / self.row()["local_filename"]
        raw = io.BytesIO()
        with zipfile.ZipFile(raw, "w", zipfile.ZIP_STORED) as target:
            target.writestr("ASCII/DEMO24Q1.txt", "unique-demo-payload")
            target.writestr("ASCII/DRUG24Q1.txt", "drug")
            target.writestr("ASCII/REAC24Q1.txt", "reac")
        damaged = bytearray(raw.getvalue())
        offset = damaged.index(b"unique-demo-payload")
        damaged[offset] ^= 0x01
        archive.write_bytes(damaged)
        self.assertEqual(dl.main(self.args("--verify-only")), 1)
        with self.log.open(encoding="utf-8") as handle:
            self.assertEqual(list(csv.DictReader(handle))[-1]["status"], "failed_invalid_zip")

    def test_successful_download_is_atomic_and_repeat_skips(self):
        self.write_manifest([self.row()])
        with patch("download_faers.urllib.request.urlopen", return_value=Response(zip_bytes())) as request:
            self.assertEqual(dl.main(self.args()), 0)
            self.assertTrue((self.out / self.row()["local_filename"]).exists())
            self.assertFalse((self.out / (self.row()["local_filename"] + ".part")).exists())
            self.assertEqual(dl.main(self.args()), 0)
            self.assertEqual(request.call_count, 1)

    def test_failed_download_and_overwrite_preservation(self):
        self.write_manifest([self.row()])
        with patch("download_faers.urllib.request.urlopen", side_effect=OSError("offline")):
            self.assertEqual(dl.main(self.args()), 1)
        self.assertFalse((self.out / self.row()["local_filename"]).exists())
        self.out.mkdir(exist_ok=True)
        final = self.out / self.row()["local_filename"]
        original = zip_bytes()
        final.write_bytes(original)
        with patch("download_faers.urllib.request.urlopen", return_value=Response(b"bad")):
            self.assertEqual(dl.main(self.args("--overwrite")), 1)
        self.assertEqual(final.read_bytes(), original)

    def test_verify_only_never_calls_http_and_missing_fails(self):
        self.write_manifest([self.row()])
        with patch("download_faers.urllib.request.urlopen") as request:
            self.assertEqual(dl.main(self.args("--verify-only")), 1)
            request.assert_not_called()

    def test_dry_run_creates_nothing(self):
        self.write_manifest([self.row()])
        with patch("download_faers.urllib.request.urlopen") as request:
            self.assertEqual(dl.main(self.args("--dry-run")), 0)
            request.assert_not_called()
        self.assertFalse(self.out.exists())
        self.assertFalse(self.log.exists())

    def test_one_failure_makes_overall_exit_nonzero(self):
        self.write_manifest([self.row(2024, 1, "q1.zip"), self.row(2024, 2, "q2.zip")])
        self.out.mkdir()
        (self.out / "q1.zip").write_bytes(zip_bytes())
        args = self.args("--verify-only")
        args[args.index("--quarters") + 1:args.index("--manifest")] = ["1", "2"]
        self.assertEqual(dl.main(args), 1)
        with self.log.open(encoding="utf-8") as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 2)


if __name__ == "__main__":
    unittest.main()
