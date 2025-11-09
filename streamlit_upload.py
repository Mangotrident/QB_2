
import streamlit as st
import pandas as pd, numpy as np, json, yaml
from pathlib import Path
from qmnet import OmicsMapper, BatchSimulator
from qle import HamiltonianBuilder, NoiseModel
from QB_1_fixed_enaqt_physics import FixedQuantumIntegrator, FixedMetricsCalculator

st.set_page_config(page_title="Upload Omics → Biomarkers", layout="centered")
st.title("Upload Omics → Physics Biomarkers (ETE, γ*, τc)")

cfg = yaml.safe_load(open(Path(__file__).parent / "config" / "default_params.yaml"))
G_spec = json.load(open(Path(__file__).parent / "data" / "templates" / "graph.json"))
import networkx as nx
G = nx.Graph()
for n in G_spec["nodes"]:
    G.add_node(int(n))
for i,j in G_spec["edges"]:
    G.add_edge(int(i), int(j))

st.markdown("**Upload files**")
up_expr = st.file_uploader("Sample omics CSV (columns: gene,value)", type=["csv"])
up_map = st.file_uploader("Mapping CSV (gene,node)", type=["csv"])

if st.button("Run Simulation", disabled=(up_expr is None or up_map is None)):
    df_expr = pd.read_csv(up_expr)
    df_map = pd.read_csv(up_map)
    gene2node = {str(r["gene"]).strip(): int(r["node"]) for _, r in df_map.iterrows()}
    expr = {str(r["gene"]).strip(): float(r["value"]) for _, r in df_expr.iterrows()}

    vals = np.array(list(expr.values()), dtype=float)
    sd = vals.std() if vals.std()!=0 else 1.0
    mu = vals.mean()
    expr = {k:(v-mu)/sd for k,v in expr.items()}

    mapper = OmicsMapper(alpha=cfg['hamiltonian']['alpha'], J_0=cfg['hamiltonian']['J_0'], J_max=cfg['hamiltonian']['J_max'])
    mapper.set_gene_mapping(gene2node)
    sim = BatchSimulator(cfg, mapper)
    res = sim.simulate_single_sample(expr, G, gamma_sweep=True)
    st.success(f"γ* = {res['gamma_star']:.4f}, ETE_peak = {res['ETE_peak']:.4f}, τc = {res['tau_c']:.4f}")
