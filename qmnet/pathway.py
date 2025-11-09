
import networkx as nx
class PathwayGraph:
    def build_etc_mini(self, n_nodes=9):
        G = nx.Graph()
        for i in range(n_nodes):
            G.add_node(i)
        edges = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),
                 (0,2),(2,4),(4,6),(6,8)]
        G.add_edges_from(edges)
        return G
