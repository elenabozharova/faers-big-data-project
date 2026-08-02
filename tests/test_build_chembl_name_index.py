import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import build_chembl_name_index as indexer  # noqa: E402


class FakeCursor:
    def __init__(self, database="chembl_36", version="ChEMBL 36"):
        self.database = database; self.version = version; self.description = []
    def __enter__(self): return self
    def __exit__(self, *_): return False
    def execute(self, query, params=None):
        if "current_database" in query: self.rows = [(self.database, "public")]
        elif "information_schema" in query: self.rows = [(name,) for name in indexer.SOURCE_TABLES]
        else:
            self.rows = [(self.version,)]; self.description = [type("D", (), {"name": "chembl_release"})()]
    def fetchone(self): return self.rows[0]
    def fetchall(self): return self.rows


class FakeConnection:
    def __init__(self, database="chembl_36", version="ChEMBL 36"):
        self.database = database; self.version = version
    def cursor(self): return FakeCursor(self.database, self.version)


class ChemblIndexTests(unittest.TestCase):
    def test_release_and_required_tables_are_verified(self):
        result = indexer.verify_chembl36(FakeConnection())
        self.assertEqual(result["database_name"], "chembl_36")
        self.assertIn("36", result["version_evidence"])

    def test_wrong_release_fails(self):
        with self.assertRaisesRegex(RuntimeError, "cannot verify"):
            indexer.verify_chembl36(FakeConnection("chembl_35", "ChEMBL 35"))


if __name__ == "__main__": unittest.main()
