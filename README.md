# decagon-preprocessing

Data engineering pipeline that builds the graph inputs required by the
[Decagon](https://snap.stanford.edu/decagon/) polypharmacy side-effect model
(Zitnik et al., 2018).

It converts four public data sources into a set of sparse adjacency matrices
and index files that a Decagon training container can consume directly.

---

## What the pipeline produces

| Output | Description |
|--------|-------------|
| `ppi_adj.npz` | Symmetric protein–protein interaction graph (UniProt-indexed) |
| `drug_target_adj.npz` | Bipartite drug–protein graph (ChEMBL drugs × UniProt proteins) |
| `combo_adj_coo.npz` | Polypharmacy side-effect triples (drug_i, drug_j, se_idx) |
| `mono_adj.npz` | Single-drug side-effect associations |
| `drug_to_idx.json` | ChEMBL molecule ID → integer index |
| `uniprot_to_idx.json` | UniProt accession → integer index |
| `se_names.json` | Side-effect index → MedDRA Preferred Term string |

All outputs land in `data/processed/decagon_inputs/` (strict drug set) or
`data/processed/decagon_inputs_all/` (broad drug set).

---

## Prerequisites

### 1. Python environment

```bash
python -m venv .venv
# Linux / macOS
source .venv/bin/activate
# Windows PowerShell
.venv\Scripts\Activate.ps1

pip install -r requirements.txt
```

### 2. Raw data files

Place the files below under `data/raw/` (or override the paths via environment
variables — see [Configuration](#configuration)).

| Variable | Default path | Where to get it |
|----------|-------------|-----------------|
| `STRING_PPI_PATH` | `data/raw/string/9606.protein.physical.links.v12.0.txt.gz` | [STRING DB](https://string-db.org) → Download → Homo sapiens → physical links |
| `UNIPROT_IDMAPPING_PATH` | `data/raw/uniprot/HUMAN_9606_idmapping.dat.gz` | Downloaded automatically on first run |
| `FAERS_ZIP_PATH` | `data/raw/faers/faers_ascii_2025Q4.zip` | [FDA FAERS](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html) → ASCII quarterly extract |
| — | ChEMBL PostgreSQL (local) | [ChEMBL downloads](https://chembl.gitbook.io/chembl-interface-documentation/downloads) → PostgreSQL dump |

### 3. ChEMBL database

Restore the ChEMBL PostgreSQL dump locally, then set `PGPASSWORD` (and
optionally `CHEMBL_DB_*` variables) so the scripts can connect.

```bash
# Example restore (adjust paths as needed)
createdb chembl_36
pg_restore -d chembl_36 chembl_36_postgresql.dmp
```

---

## Configuration

Copy `.env.example` to `.env` and fill in your values:

```bash
cp .env.example .env
```

The key variables (all optional if you follow the default `data/raw/` layout):

```
PGPASSWORD=your_db_password
CHEMBL_DB_NAME=chembl_36
CHEMBL_DB_USER=postgres
CHEMBL_DB_HOST=localhost
CHEMBL_DB_PORT=5432

STRING_PPI_PATH=data/raw/string/9606.protein.physical.links.v12.0.txt.gz
FAERS_ZIP_PATH=data/raw/faers/faers_ascii_2025Q4.zip
FAERS_DRUG_FILE=ASCII/DRUG25Q4.txt
FAERS_REAC_FILE=ASCII/REAC25Q4.txt
```

When upgrading to a new FAERS quarterly release, update `FAERS_ZIP_PATH`,
`FAERS_DRUG_FILE`, and `FAERS_REAC_FILE` to match the new filenames.

---

## Running the pipeline

Run the scripts in order from the project root. Each step depends on the
outputs of all previous steps.

### Step 1 — STRING PPI

Builds a UniProt-indexed protein–protein interaction graph from the STRING
physical-links file. Downloads the UniProt ID-mapping file automatically on
first run.

```bash
python src/process_string.py
```

Output: `data/processed/string/`

**Key configuration inside the script** (`__main__` block):

| Variable | Default | Description |
|----------|---------|-------------|
| `SCORE_MIN` | `700` | Minimum STRING combined score. 700 = high confidence (~10k proteins, ~80k edges). Use 400 for medium confidence or 150 for all physical evidence (~18k proteins, ~700k edges). |

### Step 2 — ChEMBL drug–target interactions

Extracts drug–protein pairs from a local ChEMBL PostgreSQL database.
Produces a strict set (mechanism-based) and a broad set (activity-based).

```bash
PGPASSWORD=secret python src/process_chembl.py
```

Output: `data/processed/chembl/` (strict) and `data/processed/chembl_all/` (broad)

### Step 3 — Unified Decagon graph

Unions the STRING protein space with ChEMBL targets so both graphs share one
integer protein index.

```bash
# Strict drug set (3,471 drugs, mechanism-based targets)
python src/build_decagon_graph.py strict

# Broad drug set (672k drugs, activity-based targets)
python src/build_decagon_graph.py all
```

Output: `data/processed/decagon/` and `data/processed/decagon_all/`

### Step 4 — FAERS drug name mapping

Maps FAERS drug names to ChEMBL molecule IDs using a tiered exact + fuzzy
matching strategy.

```bash
PGPASSWORD=secret python src/map_faers_to_chembl.py
```

Output: `data/processed/faers/faers_drug_to_chembl.parquet`

### Step 5 — FAERS side-effect signals

Extracts polypharmacy (combo) and single-drug (mono) side-effect co-occurrence
signals from FAERS adverse-event reports.

```bash
# Strict drug universe
python src/process_faers.py --mode strict

# Broad drug universe
python src/process_faers.py --mode all
```

Output: `data/processed/faers/combo_counts.parquet`, `mono_counts.parquet`

### Step 6 — Assemble Decagon inputs

Combines all processed data into the sparse tensor format expected by Decagon.

```bash
# Strict — top 200 side effects, high-confidence drug set
python src/build_decagon_inputs.py --mode strict --min-se-pairs 10 --top-se 200

# Broad
python src/build_decagon_inputs.py --mode all --min-se-pairs 10 --top-se 964
```

Output: `data/processed/decagon_inputs/` or `data/processed/decagon_inputs_all/`

---

## Loading outputs in a training script

```python
import scipy.sparse as sp
import numpy as np
from src.build_decagon_inputs import load_combo_adjacencies

inputs_dir = "data/processed/decagon_inputs"

ppi_adj        = sp.load_npz(f"{inputs_dir}/ppi_adj.npz")
drug_target    = sp.load_npz(f"{inputs_dir}/drug_target_adj.npz")
mono_adj       = sp.load_npz(f"{inputs_dir}/mono_adj.npz")
combo_adj_list = load_combo_adjacencies(inputs_dir)
# combo_adj_list[i] -> scipy csr_matrix (n_drugs x n_drugs) for side effect i

adj_mats = {
    (0, 0): [ppi_adj],
    (0, 1): [drug_target.T],
    (1, 0): [drug_target],
    (1, 1): combo_adj_list,
}
```

---

## Repository structure

```
decagon-preprocessing/
├── src/
│   ├── config.py                # Environment variable config (all scripts import this)
│   ├── process_string.py        # Step 1: STRING PPI -> UniProt adjacency
│   ├── process_chembl.py        # Step 2: ChEMBL drug-target extraction
│   ├── build_decagon_graph.py   # Step 3: Unified protein index + adjacency matrices
│   ├── map_faers_to_chembl.py   # Step 4: FAERS drug name -> ChEMBL ID mapping
│   ├── normalize_drug_names.py  # Drug name normalization utility
│   ├── process_faers.py         # Step 5: FAERS side-effect signal extraction
│   └── build_decagon_inputs.py  # Step 6: Assemble final Decagon tensors
├── .env.example                 # Environment variable template
├── .gitignore
├── requirements.txt
└── README.md
```

`data/` is gitignored. Raw data files must be obtained from the sources listed
above and placed under `data/raw/` (or pointed to via environment variables).

---

## Reference

Zitnik, M., Agrawal, M., & Leskovec, J. (2018). Modeling polypharmacy side
effects with graph convolutional networks. *Bioinformatics*, 34(13), i457–i466.
https://doi.org/10.1093/bioinformatics/bty294
