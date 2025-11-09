
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import yaml

from qle import HamiltonianBuilder, NoiseModel
from qmnet import PathwayGraph, OmicsMapper, BatchSimulator

from QB_1_fixed_enaqt_physics import FixedQuantumIntegrator, FixedMetricsCalculator

st.set_page_config(page_title="Quantum Bioenergetics Mapping", layout="centered")

st.title("ENAQT Demo — Energy‑Transfer Efficiency vs Dephasing")
st.write("Physics + sink flux; demonstrates noise‑assisted transport and the non‑zero optimum γ*.")

# Load config
cfg_path = Path(__file__).parent / "config" / "default_params.yaml"
config = yaml.safe_load(open(cfg_path, "r"))

# Controls
J0 = st.slider("Coupling J0", 0.005, 0.08, float(config['hamiltonian']['J_0']), 0.001)
sigma = st.slider("Static disorder σ", 0.0, 0.05, float(config['noise']['sigma']), 0.001)
k_sink = st.slider("Sink rate κ", 0.01, 0.5, float(config['sinks']['k_sink']), 0.01)
k_loss = st.slider("Loss rate Γ", 0.0, 0.05, float(config['sinks']['k_loss']), 0.001)
t_max = st.slider("Integration horizon T", 10.0, 80.0, float(config['integration']['t_max']), 1.0)

# Build graph and Hamiltonian
pathway = PathwayGraph()
graph = pathway.build_etc_mini(n_nodes=9)
hbuilder = HamiltonianBuilder(J_0=J0, J_max=config['hamiltonian']['J_max'], epsilon_0=config['hamiltonian']['epsilon_0'])
H = hbuilder.build_from_graph(graph)

noise = NoiseModel(sigma=sigma)
H = noise.add_static_disorder(H, seed=42)

n = H.shape[0]
rho0 = np.zeros((n, n), dtype=complex)
rho0[0,0] = 1.0

integrator = FixedQuantumIntegrator(dt=config['integration']['dt'], t_max=t_max)
metrics = FixedMetricsCalculator()

gamma_array = noise.create_gamma_sweep(config['noise']['gamma_range'][0], config['noise']['gamma_range'][1], config['noise']['gamma_step'])
sink_node = 4

ete_vals = []
for g in gamma_array:
    t, rho_t = integrator.evolve(H, rho0, gamma=g, k_sink=k_sink, k_loss=k_loss, sink_node=sink_node, noise_model=noise)
    ete = metrics.compute_ete_from_flux(rho_t, t, sink_node, k_sink)
    ete_vals.append(ete)

ete_vals = np.array(ete_vals)
peak_idx = int(np.argmax(ete_vals))

fig = plt.figure(figsize=(7,4.5))
plt.plot(gamma_array, ete_vals, linewidth=2)
plt.axvline(gamma_array[peak_idx], linestyle="--")
plt.axhline(ete_vals[peak_idx], linestyle="--")
plt.xlabel("Dephasing rate γ")
plt.ylabel("ETE (sink flux integral)")
plt.title("ENAQT Bell Curve")
plt.grid(True, alpha=0.3)
st.pyplot(fig)

st.write(f"**γ\* = {gamma_array[peak_idx]:.4f}**, **ETE_peak = {ete_vals[peak_idx]:.4f}**")
