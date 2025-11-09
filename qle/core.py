
import numpy as np

class HamiltonianBuilder:
    def __init__(self, J_0=0.02, J_max=0.1, epsilon_0=0.0):
        self.J_0 = J_0
        self.J_max = J_max
        self.epsilon_0 = epsilon_0

    def build_from_graph(self, graph, site_energies=None, couplings=None):
        nodes = sorted(graph.nodes())
        n = len(nodes)
        H = np.zeros((n, n), dtype=complex)
        if site_energies is None:
            for i in range(n):
                H[i, i] = self.epsilon_0
        else:
            for i in range(n):
                H[i, i] = site_energies.get(i, self.epsilon_0)
        if couplings is None:
            for i, j in graph.edges():
                H[i, j] = H[j, i] = self.J_0
        else:
            for (i, j), val in couplings.items():
                H[i, j] = H[j, i] = np.clip(val, -self.J_max, self.J_max)
        return H
