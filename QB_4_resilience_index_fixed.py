"""
QB_4: Fixed Resilience Index
Proper resilience calculation with parameter jittering.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional
import copy
import warnings


class FixedResilienceIndex:
    """
    Fixed resilience index with proper parameter perturbation.
    """
    
    def __init__(self, perturbation: float = 0.1, n_trials: int = 5):
        """
        Initialize resilience calculator.
        
        Args:
            perturbation: Fractional perturbation (±10% = 0.1)
            n_trials: Number of perturbation trials per sample
        """
        self.perturbation = perturbation
        self.n_trials = n_trials
    
    def compute_resilience(self, baseline_ete: float, perturbed_etes: List[float],
                          cohort_median: float) -> float:
        """
        Compute resilience: fraction of trials that keep sign of (ETE - median).
        
        Args:
            baseline_ete: Baseline ETE value
            perturbed_etes: List of ETE values from perturbations
            cohort_median: Cohort median ETE
            
        Returns:
            Resilience index (0-1)
        """
        baseline_sign = np.sign(baseline_ete - cohort_median)
        
        if len(perturbed_etes) == 0:
            return 0.0
        
        # Count trials that preserve sign
        preserved = 0
        for ete in perturbed_etes:
            perturbed_sign = np.sign(ete - cohort_median)
            if perturbed_sign == baseline_sign:
                preserved += 1
        
        resilience = preserved / len(perturbed_etes)
        return resilience
    
    def perturb_parameters(self, config: Dict, parameters: List[str]) -> List[Dict]:
        """
        Create perturbed configurations.
        
        Args:
            config: Base configuration
            parameters: List of parameter paths (e.g., ['hamiltonian.J_0'])
            
        Returns:
            List of perturbed configurations
        """
        perturbed_configs = []
        
        for _ in range(self.n_trials):
            perturbed_config = copy.deepcopy(config)
            
            for param_path in parameters:
                parts = param_path.split('.')
                param_dict = perturbed_config
                
                # Navigate to parameter
                for part in parts[:-1]:
                    if part not in param_dict:
                        param_dict[part] = {}
                    param_dict = param_dict[part]
                
                param_name = parts[-1]
                if param_name in param_dict:
                    base_value = param_dict[param_name]
                    # Random perturbation
                    perturbation_factor = 1.0 + np.random.uniform(
                        -self.perturbation, self.perturbation
                    )
                    param_dict[param_name] = base_value * perturbation_factor
            
            perturbed_configs.append(perturbed_config)
        
        return perturbed_configs
    
    def compute_sample_resilience(self, sample_ete: float, cohort_median: float,
                                  simulator, expression: pd.Series, graph,
                                  base_config: Dict, parameters: List[str]) -> float:
        """
        Compute resilience for a single sample.
        
        Args:
            sample_ete: Sample's baseline ETE
            cohort_median: Cohort median ETE
            simulator: BatchSimulator instance
            expression: Gene expression
            graph: Network graph
            base_config: Base configuration
            parameters: Parameters to perturb
            
        Returns:
            Resilience index
        """
        # Create perturbed configs
        perturbed_configs = self.perturb_parameters(base_config, parameters)
        
        perturbed_etes = []
        
        for perturbed_config in perturbed_configs:
            try:
                # Create new simulator with perturbed config
                from qmnet import BatchSimulator, OmicsMapper
                
                perturbed_mapper = OmicsMapper(
                    alpha=perturbed_config['hamiltonian'].get('alpha', 0.01),
                    J_0=perturbed_config['hamiltonian'].get('J_0', 0.05),
                    J_max=perturbed_config['hamiltonian'].get('J_max', 0.1)
                )
                
                # Set gene mapping
                nodes = sorted(graph.nodes())
                gene_to_node = {gene: i % len(nodes) 
                              for i, gene in enumerate(expression.index)}
                perturbed_mapper.set_gene_mapping(gene_to_node)
                
                perturbed_simulator = BatchSimulator(config=perturbed_config,
                                                   mapper=perturbed_mapper)
                
                # Run simulation
                result = perturbed_simulator.simulate_single_sample(
                    expression, graph, gamma_sweep=True
                )
                
                ete = result.get('ETE_peak', result.get('ETE', 0.0))
                perturbed_etes.append(ete)
                
            except Exception as e:
                warnings.warn(f"Perturbation trial failed: {e}")
                continue
        
        # Compute resilience
        resilience = self.compute_resilience(sample_ete, perturbed_etes, cohort_median)
        
        return resilience

