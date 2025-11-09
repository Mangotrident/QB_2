
import numpy as np

class HamiltonianBuilder:
    """Builds tight-binding Hamiltonians from a graph and parameters.
    Exposed API matches what's used in the pipeline.
    """
    def __init__(self, J_0=0.02, J_max=0.1, epsilon_0=0.0):
        self.J_0 = J_0
        self.J_max = J_max
        self.epsilon_0 = epsilon_0

    def build_from_graph(self, graph, site_energies=None, couplings=None):
        # Graph is expected to have nodes numbered 0..N-1
        nodes = sorted(graph.nodes())
        n = len(nodes)
        H = np.zeros((n, n), dtype=complex)

        # Diagonal (site energies)
        if site_energies is None:
            for i in range(n):
                H[i, i] = self.epsilon_0
        else:
            for i in range(n):
                H[i, i] = site_energies.get(i, self.epsilon_0)

        # Off-diagonal (couplings)
        if couplings is None:
            # default: J_0 for existing edges
            for i, j in graph.edges():
                H[i, j] = H[j, i] = self.J_0
        else:
            for (i, j), val in couplings.items():
                H[i, j] = H[j, i] = np.clip(val, -self.J_max, self.J_max)
        return H
