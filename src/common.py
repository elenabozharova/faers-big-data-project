import random
import networkx as nx
import matplotlib.pyplot as plt


def sample_subgraph(G: nx.Graph, n_nodes: int = 500, seed: int = 42) -> nx.Graph:
    """
    Sample an induced subgraph with up to n_nodes for visualization.

    Large graphs (STRING / Decagon) cannot be drawn fully, so we
    randomly sample nodes and visualize only that subgraph.
    """
    nodes = list(G.nodes())
    if not nodes:
        raise ValueError("Graph has no nodes to sample.")

    random.seed(seed)
    k = min(n_nodes, len(nodes))
    chosen = set(random.sample(nodes, k))
    return G.subgraph(chosen).copy()


def draw_graph(G: nx.Graph, title: str) -> None:
    """
    Draw a graph using a spring layout.
    Intended for small sampled graphs only.
    """
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42)

    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos, alpha=0.2)

    plt.title(title)
    plt.axis("off")
    plt.show()
