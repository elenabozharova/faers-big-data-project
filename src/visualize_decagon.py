"""
Visualize a small sample of the Decagon PPI network.
"""

import tarfile
import pandas as pd
import networkx as nx
from common import sample_subgraph, draw_graph

DECAGON_TAR = "data/raw/decagon/bio-decagon-ppi.tar.gz"


def load_decagon_ppi(tar_path: str) -> nx.Graph:
    with tarfile.open(tar_path, "r:gz") as tar:
        member = next(
            m for m in tar.getmembers()
            if m.isfile() and m.name.endswith("bio-decagon-ppi.csv") and "._" not in m.name
        )
        f = tar.extractfile(member)
        assert f is not None
        df = pd.read_csv(f)  # columns: Gene 1, Gene 2

    return nx.from_pandas_edgelist(df, source="Gene 1", target="Gene 2")


if __name__ == "__main__":
    G = load_decagon_ppi(DECAGON_TAR)
    print(f"Decagon PPI graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    H = sample_subgraph(G, n_nodes=500)
    draw_graph(H, "Decagon PPI (sample)")
