
import networkx as nx

class PathwayGraph:
    """Provides a minimal ETC-like graph for demos.
    Nodes are 0..8; 4 acts as 'sink' node (Complex IV) in the pipeline.
    """
    def build_etc_mini(self, n_nodes=9):
        G = nx.Graph()
        for i in range(n_nodes):
            G.add_node(i)
        # A simple chain with some shortcuts
        edges = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),
                 (0,2),(2,4),(4,6),(6,8)]
        G.add_edges_from(edges)
        return G
