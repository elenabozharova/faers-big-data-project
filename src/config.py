"""
Central configuration for FAERS temporal reporting-pattern preprocessing.

Every value that differs between machines (file paths, DB credentials) is
read from an environment variable with a sensible default that matches the
conventional data/ layout under the project root.

Setting variables
-----------------
Shell (Linux / macOS):
    export PGPASSWORD=secret
    export FAERS_ZIP_PATH=/data/faers/faers_ascii_2025Q4.zip

PowerShell (Windows):
    $env:PGPASSWORD = "secret"
    $env:FAERS_ZIP_PATH = "C:/data/faers/faers_ascii_2025Q4.zip"

Or copy .env.example -> .env, fill in values, and load with python-dotenv
before running any script:
    from dotenv import load_dotenv; load_dotenv()

See .env.example for the full list of variables and their defaults.
"""

import os

_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(_SRC_DIR)
DATA_RAW = os.path.join(PROJECT_ROOT, "data", "raw")
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")

# ---------------------------------------------------------------------------
# ChEMBL PostgreSQL connection
# ---------------------------------------------------------------------------
CHEMBL_DB_NAME     = os.environ.get("CHEMBL_DB_NAME",     "chembl_36")
CHEMBL_DB_USER     = os.environ.get("CHEMBL_DB_USER",     "postgres")
CHEMBL_DB_PASSWORD = os.environ.get("PGPASSWORD",         "")
CHEMBL_DB_HOST     = os.environ.get("CHEMBL_DB_HOST",     "localhost")
CHEMBL_DB_PORT     = int(os.environ.get("CHEMBL_DB_PORT", "5432"))

# ---------------------------------------------------------------------------
# FAERS
# ---------------------------------------------------------------------------
# The ZIP filename and inner file names change with each quarterly release.
# Update FAERS_ZIP_PATH, FAERS_DRUG_FILE, and FAERS_REAC_FILE when upgrading.
#
# Default ZIP:  data/raw/faers/faers_ascii_2025Q4.zip
# Default drug file inside ZIP:  ASCII/DRUG25Q4.txt
# Default reac file inside ZIP:  ASCII/REAC25Q4.txt
FAERS_ZIP_PATH   = os.environ.get(
    "FAERS_ZIP_PATH",
    os.path.join(DATA_RAW, "faers", "faers_ascii_2025Q4.zip"),
)
FAERS_DRUG_INNER = os.environ.get("FAERS_DRUG_FILE", "ASCII/DRUG25Q4.txt")
FAERS_REAC_INNER = os.environ.get("FAERS_REAC_FILE", "ASCII/REAC25Q4.txt")
# DEMO carries caseid/caseversion, needed by process_faers_case_level.py to
# resolve one row per real-world case (see Milestone 2 of that pipeline).
FAERS_DEMO_INNER = os.environ.get("FAERS_DEMO_FILE", "ASCII/DEMO25Q4.txt")
