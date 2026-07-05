"""
Visualize a small sample of the STRING physical PPI network (human-only).
Reads STRING .gz without loading the full file into memory.
"""

import gzip
import networkx as nx
from common import sample_subgraph, draw_graph

STRING_PATH = "data/raw/string/9606.protein.physical.links.v12.0.txt"  # <-- adjust if needed
SCORE_MIN = 700
MAX_EDGES = 500_000  # limit for fast exploration


def load_string_human_ppi(path: str) -> nx.Graph:
    G = nx.Graph()
    edges = 0

    with gzip.open(path, "rt") as f:
        header = f.readline().strip()
        print("HEADER:", header)

        for line in f:
            p1, p2, score = line.split()
            score = int(score)

            # keep human only
            if not (p1.startswith("9606.") and p2.startswith("9606.")):
                continue

            # confidence threshold
            if score < SCORE_MIN:
                continue

            G.add_edge(p1, p2, weight=score)
            edges += 1

            if edges >= MAX_EDGES:
                break

    return G


if __name__ == "__main__":
    G = load_string_human_ppi(STRING_PATH)
    print(f"STRING graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    H = sample_subgraph(G, n_nodes=500)
    draw_graph(H, f"STRING human physical PPI (sample) | score>={SCORE_MIN}")
