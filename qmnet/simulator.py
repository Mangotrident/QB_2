
import numpy as np
from .pathway import PathwayGraph
from qle import HamiltonianBuilder, NoiseModel
from QB_1_fixed_enaqt_physics import FixedQuantumIntegrator, FixedMetricsCalculator

class BatchSimulator:
    """Simplified batch simulator that runs a gamma sweep and reports ETE_peak."""
    def __init__(self, config, mapper):
        self.config = config
        self.mapper = mapper
        self.hbuilder = HamiltonianBuilder(
            J_0=config['hamiltonian']['J_0'],
            J_max=config['hamiltonian']['J_max'],
            epsilon_0=config['hamiltonian']['epsilon_0']
        )
        self.noise = NoiseModel(sigma=config['noise']['sigma'])
        self.integrator = FixedQuantumIntegrator(
            dt=config['integration']['dt'],
            t_max=config['integration']['t_max']
        )
        self.metrics = FixedMetricsCalculator()

    def simulate_single_sample(self, expression, graph, gamma_sweep=True):
        nodes = sorted(graph.nodes())
        params = self.mapper.map_sample(expression, graph, nodes)
        H = self.hbuilder.build_from_graph(graph, params['site_energies'], params['couplings'])
        H = self.noise.add_static_disorder(H)

        n = H.shape[0]
        rho0 = np.zeros((n, n), dtype=complex)
        rho0[0, 0] = 1.0

        sink_node = 4
        k_sink = self.config['sinks']['k_sink']
        k_loss = self.config['sinks']['k_loss']

        if gamma_sweep:
            gamma_array = self.noise.create_gamma_sweep(
                self.config['noise']['gamma_range'][0],
                self.config['noise']['gamma_range'][1],
                self.config['noise']['gamma_step']
            )
        else:
            gamma_array = [self.config['noise']['gamma_default']]

        ete_vals = []
        tau_vals = []
        for gamma in gamma_array:
            t, rho_t = self.integrator.evolve(H, rho0, gamma=gamma, k_sink=k_sink, k_loss=k_loss, sink_node=sink_node, noise_model=self.noise)
            ete = self.metrics.compute_ete_from_flux(rho_t, t, sink_node, k_sink)
            tau = self.metrics.compute_tau_c(rho_t, t)
            ete_vals.append(ete)
            tau_vals.append(tau)

        ete_vals = np.array(ete_vals)
        peak_idx = int(np.argmax(ete_vals))
        return {
            'ETE_peak': float(ete_vals[peak_idx]),
            'gamma_star': float(gamma_array[peak_idx]),
            'tau_c': float(tau_vals[peak_idx]),
        }
