"""
QB_5: Fixed Edge Sensitivity
Proper edge sensitivity with ranking and Complex I edge detection.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt
import warnings


class FixedEdgeSensitivity:
    """
    Fixed edge sensitivity analysis.
    """
    
    def __init__(self):
        self.complexI_nodes = [0]  # Complex I is typically node 0 in etc_mini
        self.complexIII_nodes = [2]  # Complex III is typically node 2
    
    def compute_edge_sensitivity(self, graph, expression: pd.Series,
                                baseline_ete: float, mapper, hamiltonian_builder,
                                integrator, metrics_calc, noise_model,
                                sink_node: int = 4, perturbation: float = 0.05) -> Dict:
        """
        Compute edge sensitivity by perturbing J[i,j].
        
        For each edge (i,j): set J'[i,j] = 1.05 * J[i,j], recompute ETE_peak.
        Score: ΔETE_peak = ETE'_peak - ETE_peak
        """
        edges = list(graph.edges())
        sensitivities = {}
        
        # Get baseline parameters
        nodes = sorted(graph.nodes())
        base_params = mapper.map_sample(expression, graph, nodes)
        base_H = hamiltonian_builder.build_from_graph(
            graph, site_energies=base_params['site_energies'],
            couplings=base_params['couplings']
        )
        
        # Get parameters
        k_sink = 0.10
        k_loss = 0.01
        gamma = 0.02  # Use moderate gamma
        
        # Initial state
        n = base_H.shape[0]
        rho0 = np.zeros((n, n), dtype=complex)
        rho0[0, 0] = 1.0
        
        for edge in edges:
            # Create perturbed couplings
            perturbed_couplings = base_params['couplings'].copy()
            
            # Perturb edge coupling
            if edge in perturbed_couplings:
                perturbed_couplings[edge] = perturbed_couplings[edge] * (1.0 + perturbation)
            elif (edge[1], edge[0]) in perturbed_couplings:
                reverse_edge = (edge[1], edge[0])
                perturbed_couplings[reverse_edge] = perturbed_couplings[reverse_edge] * (1.0 + perturbation)
            else:
                sensitivities[edge] = 0.0
                continue
            
            # Build perturbed Hamiltonian
            perturbed_H = hamiltonian_builder.build_from_graph(
                graph, site_energies=base_params['site_energies'],
                couplings=perturbed_couplings
            )
            
            # Run simulation
            try:
                t_array, rho_t = integrator.evolve(
                    perturbed_H, rho0, gamma=gamma,
                    k_sink=k_sink, k_loss=k_loss,
                    sink_node=sink_node, noise_model=noise_model
                )
                
                ete_perturbed = metrics_calc.compute_ete_from_flux(
                    rho_t, t_array, sink_node, k_sink
                )
                
                # Sensitivity = ΔETE
                delta_ete = ete_perturbed - baseline_ete
                sensitivities[edge] = float(delta_ete)
                
            except Exception as e:
                warnings.warn(f"Edge {edge} sensitivity failed: {e}")
                sensitivities[edge] = 0.0
        
        return sensitivities
    
    def rank_edges(self, sensitivities: Dict, graph) -> pd.DataFrame:
        """
        Rank edges by sensitivity.
        """
        edges = list(sensitivities.keys())
        sens_values = [sensitivities[e] for e in edges]
        
        # Sort by absolute sensitivity
        sorted_indices = np.argsort(np.abs(sens_values))[::-1]
        
        ranking_data = []
        for idx in sorted_indices:
            edge = edges[idx]
            sensitivity = sens_values[idx]
            
            # Check if Complex I/III edge
            is_complexI = edge[0] in self.complexI_nodes or edge[1] in self.complexI_nodes
            is_complexIII = edge[0] in self.complexIII_nodes or edge[1] in self.complexIII_nodes
            
            # Get label
            if hasattr(graph, 'edge_labels') and edge in graph.edge_labels:
                label = graph.edge_labels[edge]
            else:
                label = f"{edge[0]}→{edge[1]}"
            
            ranking_data.append({
                'edge': str(edge),
                'label': label,
                'sensitivity': float(sensitivity),
                'abs_sensitivity': float(abs(sensitivity)),
                'is_complexI': is_complexI,
                'is_complexIII': is_complexIII,
                'rank': len(ranking_data) + 1
            })
        
        return pd.DataFrame(ranking_data)
    
    def check_complexI_in_top25(self, ranking_df: pd.DataFrame) -> bool:
        """
        Check if Complex I edges appear in top 25%.
        """
        top_25_percent = max(1, len(ranking_df) // 4)
        top_edges = ranking_df.head(top_25_percent)
        
        complexI_in_top = top_edges['is_complexI'].any() or top_edges['is_complexIII'].any()
        
        return complexI_in_top
    
    def plot_edge_sensitivity(self, ranking_df: pd.DataFrame, output_path: Path):
        """
        Plot edge sensitivity ranking.
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        top_n = min(10, len(ranking_df))
        top_edges = ranking_df.head(top_n)
        
        y_pos = np.arange(len(top_edges))
        colors = ['red' if (row['is_complexI'] or row['is_complexIII']) else 'steelblue'
                 for _, row in top_edges.iterrows()]
        
        ax.barh(y_pos, top_edges['abs_sensitivity'], color=colors, alpha=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_edges['label'], fontsize=10)
        ax.set_xlabel('|ΔETE_peak|', fontsize=12, fontweight='bold')
        ax.set_title('Edge Sensitivity Ranking', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='x')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', alpha=0.7, label='Complex I/III'),
            Patch(facecolor='steelblue', alpha=0.7, label='Other')
        ]
        ax.legend(handles=legend_elements, loc='lower right')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

