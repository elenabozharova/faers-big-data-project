import zipfile
import pandas as pd

ZIP_PATH = "data/raw/faers/faers_ascii_2025Q4.zip"

# ---------- Helpers ----------
def read_faers_table(inner_path: str, nrows: int = 100) -> pd.DataFrame:
    """Read a FAERS ASCII table from inside the quarterly ZIP."""
    with zipfile.ZipFile(ZIP_PATH) as z:
        with z.open(inner_path) as f:
            # Your files use '$' delimiter (confirmed from DRUG header)
            return pd.read_csv(f, sep="$", dtype=str, nrows=nrows, low_memory=False)

def preview(df: pd.DataFrame, title: str, cols: list[str] | None = None) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)
    if cols:
        missing = [c for c in cols if c not in df.columns]
        if missing:
            print("Missing columns in file:", missing)
            print("Available columns:", list(df.columns))
        use = [c for c in cols if c in df.columns]
        df = df[use]
    print(df.head(100).to_string(index=False))

    # Save a CSV snapshot so VS Code opens it as a nice grid
    safe_name = title.lower().replace(" ", "_").replace("/", "_").replace("-", "_")
    out_path = f"data/processed/{safe_name}_top100.csv"

    print(out_path)
    df.head(100).to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")

# ---------- Table configs ----------
TABLES = {
    "FAERS DRUG (relevant fields)": {
        "inner": "ASCII/DRUG25Q4.txt",
        "cols": ["primaryid", "caseid", "drug_seq", "role_cod", "drugname", "prod_ai"]
    },
    "FAERS REAC (relevant fields)": {
        "inner": "ASCII/REAC25Q4.txt",
        # We'll attempt typical columns; script will warn if names differ
        "cols": ["primaryid", "caseid", "pt"]  # pt = MedDRA Preferred Term (often)
    },
    "FAERS DEMO (relevant fields)": {
        "inner": "ASCII/DEMO25Q4.txt",
        "cols": ["primaryid", "caseid", "age", "age_cod", "sex", "wt", "wt_cod", "reporter_country"]
    },
    # Optional: indications, useful for source-data inspection
    "FAERS INDI (optional)": {
        "inner": "ASCII/INDI25Q4.txt",
        "cols": ["primaryid", "caseid", "indi_pt"]  # often indi_pt or indi; script will warn if mismatch
    }
}

if __name__ == "__main__":
    # Make sure output folder exists
    import os
    os.makedirs("data/processed", exist_ok=True)

    for title, cfg in TABLES.items():
        df = read_faers_table(cfg["inner"], nrows=100)
        preview(df, title, cols=cfg["cols"])
