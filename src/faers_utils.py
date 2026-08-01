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


def load_faers_table(zip_path: str, inner_path: str, max_rows: int | None = None) -> pd.DataFrame:
    """Load a dollar-delimited FAERS ASCII table from a quarterly ZIP."""
    with zipfile.ZipFile(zip_path) as archive:
        with archive.open(inner_path) as source:
            df = pd.read_csv(
                source,
                sep="$",
                dtype=str,
                low_memory=False,
                nrows=max_rows,
            )
    df.columns = df.columns.str.lower().str.strip()
    print(f"  {inner_path}: {len(df):,} rows, columns={list(df.columns)}")
    return df
