
import numpy as np

class OmicsMapper:
    """Maps gene expression to site energies and couplings.
    This is a lightweight, local version that matches the pipeline API.
    """
    def __init__(self, alpha=0.01, J_0=0.02, J_max=0.1):
        self.alpha = alpha
        self.J_0 = J_0
        self.J_max = J_max
        self._gene_to_node = {}

    def set_gene_mapping(self, gene_to_node):
        self._gene_to_node = dict(gene_to_node)

    def map_sample(self, expression, graph, nodes):
        # Site energies: average expression mapped per node (z-scored expr expected)
        site_accum = {n: [] for n in nodes}
        for gene, val in expression.items():
            if gene in self._gene_to_node:
                site_accum[self._gene_to_node[gene]].append(float(val))

        site_energies = {}
        for n in nodes:
            vals = site_accum[n]
            mean_val = float(np.mean(vals)) if len(vals) else 0.0
            site_energies[n] = self.alpha * mean_val

        # Couplings: base J_0 for each edge
        couplings = {}
        for (i, j) in graph.edges():
            couplings[(i, j)] = self.J_0
        return {'site_energies': site_energies, 'couplings': couplings}
