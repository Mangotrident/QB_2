
import numpy as np
class OmicsMapper:
    def __init__(self, alpha=0.01, J_0=0.02, J_max=0.1):
        self.alpha = alpha
        self.J_0 = J_0
        self.J_max = J_max
        self._gene_to_node = {}

    def set_gene_mapping(self, gene_to_node):
        self._gene_to_node = dict(gene_to_node)

    def map_sample(self, expression, graph, nodes):
        site_accum = {n: [] for n in nodes}
        for gene, val in expression.items():
            if gene in self._gene_to_node:
                site_accum[self._gene_to_node[gene]].append(float(val))
        site_energies = {n: self.alpha*(float(np.mean(v)) if v else 0.0) for n, v in site_accum.items()}
        couplings = {(i, j): self.J_0 for (i, j) in graph.edges()}
        return {'site_energies': site_energies, 'couplings': couplings}
