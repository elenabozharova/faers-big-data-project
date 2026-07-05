"""
Visualize ChEMBL drug–target network (drug ↔ protein).

Assumes ChEMBL is restored into PostgreSQL (chembl_36).
Draws a sampled bipartite subgraph for readability.
"""

import random
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from sqlalchemy import create_engine, text

# ---- CHANGE THESE ----
USER = "postgres"
PASSWORD = "1886" # --  postgres password
HOST = "localhost"
PORT = 5432
DB = "chembl_36"

# ---- Controls ----
MAX_ROWS = 200_000     # how many rows to pull from DB for exploration
N_DRUGS = 200          # how many drugs to sample for visualization
SEED = 42

engine = create_engine(f"postgresql+psycopg2://{USER}:{PASSWORD}@{HOST}:{PORT}/{DB}")

QUERY = text("""
SELECT
    md.chembl_id AS drug,
    cs.accession          AS protein_uniprot
FROM drug_mechanism dm
JOIN molecule_dictionary md ON dm.molregno = md.molregno
JOIN target_dictionary td ON dm.tid = td.tid
JOIN target_components tc ON td.tid = tc.tid
JOIN component_sequences cs ON tc.component_id = cs.component_id
WHERE cs.tax_id = 9606
  AND cs.accession IS NOT NULL
  AND md.chembl_id IS NOT NULL
LIMIT :max_rows;
""")

def build_bipartite_graph(df: pd.DataFrame) -> nx.Graph:
    G = nx.Graph()

    # tag nodes by type (helps later)
    for d in df["drug"].unique():
        G.add_node(f"D:{d}", kind="drug")
    for p in df["protein_uniprot"].unique():
        G.add_node(f"P:{p}", kind="protein")

    # edges drug-protein
    for row in df.itertuples(index=False):
        G.add_edge(f"D:{row.drug}", f"P:{row.protein_uniprot}")
    return G


def sample_drug_induced_subgraph(G: nx.Graph, n_drugs: int, seed: int) -> nx.Graph:
    random.seed(seed)
    drug_nodes = [n for n, data in G.nodes(data=True) if data.get("kind") == "drug"]
    chosen_drugs = set(random.sample(drug_nodes, min(n_drugs, len(drug_nodes))))

    # include all proteins connected to those drugs
    chosen = set(chosen_drugs)
    for d in chosen_drugs:
        chosen.update(G.neighbors(d))

    return G.subgraph(chosen).copy()


def draw_bipartite(H: nx.Graph, title: str) -> None:
    drugs = [n for n, data in H.nodes(data=True) if data.get("kind") == "drug"]
    prots = [n for n, data in H.nodes(data=True) if data.get("kind") == "protein"]

    # bipartite layout (no SciPy required)
    pos = {}
    for i, n in enumerate(drugs):
        pos[n] = (0, i)
    for i, n in enumerate(prots):
        pos[n] = (1, i)

    plt.figure(figsize=(12, 10))
    nx.draw_networkx_nodes(H, pos, nodelist=drugs, node_size=30)
    nx.draw_networkx_nodes(H, pos, nodelist=prots, node_size=30)
    nx.draw_networkx_edges(H, pos, alpha=0.2)
    plt.title(title)
    plt.axis("off")
    plt.show()


if __name__ == "__main__":
    df = pd.read_sql(QUERY, engine, params={"max_rows": MAX_ROWS})

    print("Pulled rows:", len(df))
    print("Unique drugs:", df["drug"].nunique())
    print("Unique proteins:", df["protein_uniprot"].nunique())
    print(df.head(5))

    G = build_bipartite_graph(df)
    print(f"\nChEMBL bipartite graph (from sample): {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    H = sample_drug_induced_subgraph(G, n_drugs=N_DRUGS, seed=SEED)
    print(f"Visualization subgraph: {H.number_of_nodes():,} nodes, {H.number_of_edges():,} edges")

    draw_bipartite(H, f"ChEMBL Drug–Target (sample) | drugs={min(N_DRUGS, df['drug'].nunique())}, rows={MAX_ROWS:,}")
