import io
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import faers_utils  # noqa: E402


def make_zip(directory: Path, inner_name: str, body_bytes: bytes) -> Path:
    zip_path = directory / "sample.zip"
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr(f"ascii/{inner_name}", body_bytes)
    return zip_path


class FaersUtilsEncodingTests(unittest.TestCase):
    def test_clean_utf8_loads_normally(self):
        body = "primaryid$caseid$drugname\n1$1$aspirin\n".encode("utf-8")
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            zip_path = make_zip(directory, "DRUG13Q1.txt", body)
            df = faers_utils.load_faers_table(str(zip_path), "ascii/DRUG13Q1.txt")
            self.assertEqual(df["drugname"].tolist(), ["aspirin"])

    def test_windows_1252_soft_hyphen_falls_back_and_decodes(self):
        # 0xAD (soft hyphen) is invalid as a standalone UTF-8 byte but valid cp1252.
        body = "primaryid$caseid$drugname\n1$1$co\xadamoxiclav\n".encode("cp1252")
        with self.assertRaises(UnicodeDecodeError):
            body.decode("utf-8")
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            zip_path = make_zip(directory, "DRUG19Q3.txt", body)
            df = faers_utils.load_faers_table(str(zip_path), "ascii/DRUG19Q3.txt")
            self.assertEqual(df["drugname"].tolist(), ["co\xadamoxiclav"])


if __name__ == "__main__":
    unittest.main()
