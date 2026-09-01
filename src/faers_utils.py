"""Shared FAERS ingestion and reaction-filtering utilities."""

import zipfile

import pandas as pd


# Administrative, medication-error, and effectiveness terms that are not
# treated as adverse reactions in the case-level reporting-pattern pipeline.
REACTION_PT_EXCLUDE: frozenset[str] = frozenset({
    "drug ineffective",
    "off label use",
    "product use in unapproved indication",
    "intentional overdose",
    "accidental overdose",
    "drug interaction",
    "drug dose omission",
    "wrong drug administered",
    "drug dispensing error",
    "expired product administered",
    "inappropriate schedule of product administration",
    "product administered to patient of inappropriate age",
    "medication error",
    "no adverse event",
    "therapeutic response unexpected",
    "incorrect dose administered",
    "drug administration error",
    "product quality issue",
    "underdose",
    "overdose",
    "drug level increased",
    "drug level decreased",
    "drug level above therapeutic",
    "drug level below therapeutic",
    "drug resistance",
    "lack of efficacy",
})


def resolve_member(zf: zipfile.ZipFile, wanted: str) -> str:
    """Resolve a ZIP member case-insensitively (older FAERS archives use
    lowercase 'ascii/' directories; the manifest declares 'ASCII/')."""
    wanted = wanted.replace("\\", "/")
    exact = [n for n in zf.namelist() if n.replace("\\", "/") == wanted]
    found = exact or [n for n in zf.namelist() if n.replace("\\", "/").casefold() == wanted.casefold()]
    if len(found) != 1:
        raise KeyError(f"{'missing' if not found else 'ambiguous'} ZIP member {wanted!r}")
    return found[0]


def load_faers_table(zip_path: str, inner_path: str, max_rows: int | None = None) -> pd.DataFrame:
    """Load a dollar-delimited FAERS ASCII table from a quarterly ZIP."""
    with zipfile.ZipFile(zip_path) as archive:
        resolved_path = resolve_member(archive, inner_path)
        try:
            with archive.open(resolved_path) as source:
                df = pd.read_csv(
                    source, sep="$", dtype=str, low_memory=False, nrows=max_rows, encoding="utf-8",
                )
        except UnicodeDecodeError:
            # A handful of pre-2020 FAERS quarters carry Windows-1252 bytes (e.g. a
            # soft hyphen, 0xAD) in free-text fields instead of valid UTF-8. Only
            # retry on a confirmed UTF-8 failure -- never silently reinterpret bytes
            # that already decoded correctly as UTF-8.
            print(f"  [WARN] {inner_path}: not valid UTF-8, retrying as cp1252")
            with archive.open(resolved_path) as source:
                df = pd.read_csv(
                    source, sep="$", dtype=str, low_memory=False, nrows=max_rows, encoding="cp1252",
                )
    df.columns = df.columns.str.lower().str.strip()
    print(f"  {inner_path}: {len(df):,} rows, columns={list(df.columns)}")
    return df
