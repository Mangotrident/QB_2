"""
QB_1: Fixed ENAQT Physics
Fixes the quantum integrator to produce proper bell curve with noise-assisted peak.
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import Optional, Tuple, Dict
import warnings


class FixedQuantumIntegrator:
    """
    Fixed quantum integrator that properly implements ENAQT.
    
    Key fixes:
    - Proper sink flux integration (not final population)
    - Correct dephasing implementation
    - Proper time integration
    """
    
    def __init__(self, dt: float = 0.2, t_max: float = 40.0, rtol: float = 1e-6, atol: float = 1e-8):
        self.dt = dt
        self.t_max = t_max
        self.rtol = rtol
        self.atol = atol
    
    def _vectorize_density_matrix(self, rho: np.ndarray) -> np.ndarray:
        """Convert density matrix to vector form."""
        return rho.flatten()
    
    def _unvectorize_density_matrix(self, rho_vec: np.ndarray, n: int) -> np.ndarray:
        """Convert vector back to density matrix."""
        return rho_vec.reshape(n, n)
    
    def _lindblad_rhs(self, t: float, rho_vec: np.ndarray, H: np.ndarray, gamma: float,
                     k_sink: float, k_loss: float, sink_node: Optional[int],
                     noise_model) -> np.ndarray:
        """Right-hand side of Lindblad equation - FIXED."""
        n = H.shape[0]
        rho = self._unvectorize_density_matrix(rho_vec, n)
        
        # Hamiltonian evolution: -i[H, ρ]
        H_comm = -1j * (H @ rho - rho @ H)
        
        # Dephasing: -γ/2 * Σ_i [|i⟩⟨i|, [|i⟩⟨i|, ρ]]
        drho_dephase = np.zeros_like(rho, dtype=complex)
        for i in range(n):
            P_i = np.zeros((n, n), dtype=complex)
            P_i[i, i] = 1.0
            comm1 = P_i @ rho - rho @ P_i
            comm2 = P_i @ comm1 - comm1 @ P_i
            drho_dephase -= (gamma / 2.0) * comm2
        
        # Sink channel (energy capture at specific node)
        drho_sink = np.zeros_like(rho, dtype=complex)
        if sink_node is not None and sink_node < n:
            # Lindblad operator: L_sink = sqrt(k_sink) |sink⟩⟨sink|
            L_sink = np.zeros((n, n), dtype=complex)
            L_sink[sink_node, sink_node] = np.sqrt(k_sink)
            L_dag_L = L_sink.T.conj() @ L_sink
            drho_sink = (L_sink @ rho @ L_sink.T.conj() - 
                        0.5 * (L_dag_L @ rho + rho @ L_dag_L))
        
        # Loss channel (uniform dissipation)
        drho_loss = -k_loss * rho if k_loss > 0 else np.zeros_like(rho, dtype=complex)
        
        # Total derivative
        drho = H_comm + drho_dephase + drho_sink + drho_loss
        
        return self._vectorize_density_matrix(drho)
    
    def evolve(self, H: np.ndarray, rho0: np.ndarray, gamma: float = 0.01,
              k_sink: float = 0.10, k_loss: float = 0.01, sink_node: Optional[int] = None,
              noise_model=None, seed: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Evolve quantum system - FIXED for proper ENAQT.
        """
        if seed is not None:
            np.random.seed(seed)
        
        if noise_model is None:
            from qle.noise import NoiseModel
            noise_model = NoiseModel()
        
        n = H.shape[0]
        
        # Validate initial state
        trace0 = np.trace(rho0).real
        if abs(trace0 - 1.0) > 1e-4:
            raise ValueError(f"Invalid initial density matrix: trace = {trace0}")
        
        # Time span
        t_span = (0.0, self.t_max)
        t_eval = np.arange(0, self.t_max + self.dt, self.dt)
        
        # Initial condition (vectorized)
        rho0_vec = self._vectorize_density_matrix(rho0)
        
        # ODE function
        def rhs(t, y):
            return self._lindblad_rhs(t, y, H, gamma, k_sink, k_loss, sink_node, noise_model)
        
        # Solve ODE
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solve_ivp(
                rhs, t_span, rho0_vec, t_eval=t_eval,
                method='RK45', rtol=self.rtol, atol=self.atol, dense_output=False
            )
        
        if not sol.success:
            raise RuntimeError(f"ODE integration failed: {sol.message}")
        
        # Reshape solution
        n_times = len(sol.t)
        rho_t = np.zeros((n_times, n, n), dtype=complex)
        
        for i in range(n_times):
            rho_t[i] = self._unvectorize_density_matrix(sol.y[:, i], n)
        
        return sol.t, rho_t


class FixedMetricsCalculator:
    """
    Fixed metrics calculator with proper ETE, τc, and QLS computation.
    """
    
    def compute_ete_from_flux(self, rho_t: np.ndarray, t_array: np.ndarray,
                              sink_node: int, k_sink: float = 0.10) -> float:
        """
        Compute ETE from sink flux (FIXED - not final population).
        
        ETE = ∫₀^T k_sink * ρ_sink(t) dt
        """
        n_times = len(t_array)
        sink_flux = np.zeros(n_times)
        
        for i in range(n_times):
            rho = rho_t[i]
            sink_flux[i] = k_sink * rho[sink_node, sink_node].real
        
        # Integrate flux over time
        ete = np.trapz(sink_flux, t_array)
        
        return ete
    
    def compute_coherence_curve(self, rho_t: np.ndarray) -> np.ndarray:
        """
        Compute coherence curve: C(t) = sum_{i≠j} |ρ_ij(t)|^2.
        """
        n_times = rho_t.shape[0]
        coherence = np.zeros(n_times)
        
        for i in range(n_times):
            rho = rho_t[i]
            # Sum of off-diagonal elements squared
            off_diag = rho - np.diag(np.diag(rho))
            coherence[i] = np.sum(np.abs(off_diag)**2)
        
        return coherence
    
    def compute_tau_c(self, rho_t: np.ndarray, t_array: np.ndarray) -> float:
        """
        Compute coherence lifetime: τc = ∫ C(t) dt.
        """
        coherence = self.compute_coherence_curve(rho_t)
        tau_c = np.trapz(coherence, t_array)
        return tau_c
    
    def compute_qls(self, ete_peak: float, tau_c: float, resilience: float = 0.5,
                   ete_min: float = 0.0, ete_max: float = 1.0,
                   tau_c_min: float = 0.0, tau_c_max: float = 100.0) -> float:
        """
        Compute QLS: QLS = 0.5 * Ê + 0.3 * T̂ + 0.2 * R
        
        where Ê and T̂ are min-max normalized ETE and τc.
        """
        # Normalize ETE
        if ete_max > ete_min:
            ete_norm = (ete_peak - ete_min) / (ete_max - ete_min)
        else:
            ete_norm = 0.0
        ete_norm = np.clip(ete_norm, 0.0, 1.0)
        
        # Normalize τc
        if tau_c_max > tau_c_min:
            tau_c_norm = (tau_c - tau_c_min) / (tau_c_max - tau_c_min)
        else:
            tau_c_norm = 0.0
        tau_c_norm = np.clip(tau_c_norm, 0.0, 1.0)
        
        # QLS
        qls = 0.5 * ete_norm + 0.3 * tau_c_norm + 0.2 * resilience
        
        return qls
    
    def check_integrity(self, rho_t: np.ndarray) -> Dict[str, float]:
        """
        Check quantum integrity: trace and PSD.
        """
        n_times = rho_t.shape[0]
        max_trace_error = 0.0
        min_eigenvalue = 0.0
        
        for i in range(n_times):
            rho = rho_t[i]
            
            # Trace check
            trace = np.trace(rho).real
            trace_error = abs(trace - 1.0)
            max_trace_error = max(max_trace_error, trace_error)
            
            # PSD check
            rho_herm = (rho + rho.T.conj()) / 2.0
            eigenvals = np.linalg.eigvalsh(rho_herm)
            min_eigenvalue = min(min_eigenvalue, np.min(eigenvals))
        
        return {
            'max_trace_error': float(max_trace_error),
            'min_eigenvalue': float(min_eigenvalue),
            'pass_qc': max_trace_error < 1e-6 and min_eigenvalue > -1e-8
        }

