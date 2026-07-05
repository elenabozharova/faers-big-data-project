"""
SQLite-backed HTTP response cache.

Every URL fetched by the pipeline is stored as JSON so re-runs never
re-hit the API.  The key is the full URL (SHA-256 hashed for compactness).
Thread-safe for single-process sequential access.
"""
import hashlib
import json
import os
import sqlite3
from typing import Optional


class ResponseCache:
    def __init__(self, db_path: str) -> None:
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        self._conn = sqlite3.connect(db_path, check_same_thread=False)
        self._conn.execute("""
            CREATE TABLE IF NOT EXISTS cache (
                key       TEXT PRIMARY KEY,
                value     TEXT NOT NULL,
                url       TEXT NOT NULL,
                hit_count INTEGER NOT NULL DEFAULT 0,
                created   TEXT DEFAULT (datetime('now'))
            )
        """)
        self._conn.commit()
        self._hits = 0
        self._misses = 0

    # ------------------------------------------------------------------
    def _key(self, url: str) -> str:
        return hashlib.sha256(url.encode()).hexdigest()

    def get(self, url: str) -> Optional[dict]:
        row = self._conn.execute(
            "SELECT value FROM cache WHERE key = ?", (self._key(url),)
        ).fetchone()
        if row is None:
            self._misses += 1
            return None
        self._hits += 1
        self._conn.execute(
            "UPDATE cache SET hit_count = hit_count + 1 WHERE key = ?",
            (self._key(url),),
        )
        return json.loads(row[0])

    def set(self, url: str, data: dict) -> None:
        self._conn.execute(
            "INSERT OR REPLACE INTO cache (key, url, value) VALUES (?, ?, ?)",
            (self._key(url), url, json.dumps(data)),
        )
        self._conn.commit()

    # ------------------------------------------------------------------
    def stats(self) -> dict:
        total = self._hits + self._misses
        n_cached = self._conn.execute("SELECT COUNT(*) FROM cache").fetchone()[0]
        return {
            "cache_hits":   self._hits,
            "cache_misses": self._misses,
            "hit_rate":     self._hits / total if total else 0.0,
            "entries":      n_cached,
        }

    def close(self) -> None:
        self._conn.close()
