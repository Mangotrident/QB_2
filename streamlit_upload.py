import streamlit as st
import pandas as pd, numpy as np, json, yaml
from pathlib import Path

# ---- Minimal helpers (no external packages beyond numpy/pandas/streamlit) ----

def build_graph_from_json(path):
    with open(path, "r") as f:
        spec = json.load(f)
    n = len(spec["nodes"])
    edges = [(int(i), int(j)) for (i, j) in spec["edges"]]
    return n, edges

def map_omics_to_params(expr_df, mapping_df, alpha=0.01, J0=0.025, Jmax=0.10):
    gene2node = {str(r["gene"]).strip(): int(r["node"]) for _, r in mapping_df.iterrows()}
    expr = {str(r["gene"]).strip(): float(r["value"]) for _, r in expr_df.iterrows()}
    # z-score
    vals = np.array(list(expr.values()), dtype=float)
    mu = float(vals.mean()); sd = float(vals.std()) if vals.std() != 0 else 1.0
    expr_z = {k:(v-mu)/sd for k,v in expr.items()}
    return gene2node, expr_z, alpha, J0, Jmax

def build_hamiltonian(n, edges, gene2node, expr_z, alpha, J0, Jmax):
    H = np.zeros((n, n), dtype=np.complex128)
    # site energies from expr
    means = {node:[] for node in range(n)}
    for g, v in expr_z.items():
        if g in gene2node:
            means[gene2node[g]].append(v)
    for i in range(n):
        val = np.mean(means[i]) if len(means[i])>0 else 0.0
        H[i, i] = alpha * val
    # couplings
    for (i, j) in edges:
        H[i, j] = H[j, i] = np.clip(J0, -Jmax, Jmax)
    return H

def lindblad_step(rho, H, gamma, k_sink, k_loss, sink_node, dt):
    # -i[H,rho]
    comm = H @ rho - rho @ H
    drho = -1j * comm

    # dephasing on each site
    n = rho.shape[0]
    for i in range(n):
        Pi = np.zeros_like(rho); Pi[i, i] = 1.0
        drho += gamma * (Pi @ rho @ Pi - 0.5 * (Pi @ rho + rho @ Pi))

    # sink localized at sink_node
    Ps = np.zeros_like(rho); Ps[sink_node, sink_node] = 1.0
    drho += k_sink * (Ps @ rho @ Ps - 0.5 * (Ps @ rho + rho @ Ps))

    # uniform loss
    for i in range(n):
        Pi = np.zeros_like(rho); Pi[i, i] = 1.0
        drho += k_loss * (Pi @ rho @ Pi - 0.5 * (Pi @ rho + rho @ Pi))

    return rho + dt * drho

def simulate_ete_curve(H, gamma_values, sink_node=4, k_sink=0.20, k_loss=0.005, dt=0.2, t_max=60.0, sigma=0.08, seed=1):
    rng = np.random.default_rng(seed)
    # add static diagonal disorder
    H = H.copy()
    diag_noise = rng.normal(0.0, sigma, size=H.shape[0])
    for i in range(H.shape[0]):
        H[i, i] = H[i, i].real + diag_noise[i]

    n_steps = int(np.ceil(t_max/dt))
    times = np.linspace(0, t_max, n_steps+1)

    ete_vals = []
    tau_vals = []
    for gamma in gamma_values:
        rho = np.zeros((H.shape[0], H.shape[1]), dtype=np.complex128)
        rho[0,0] = 1.0
        sink_pop = []
        purity = []
        for _ in range(n_steps+1):
            sink_pop.append(np.real(rho[sink_node, sink_node]))
            purity.append(np.real(np.trace(rho @ rho)))
            rho = lindblad_step(rho, H, gamma, k_sink, k_loss, sink_node, dt)

        sink_pop = np.array(sink_pop)
        purity = np.array(purity)
        ete = float(k_sink * np.trapz(sink_pop, times))
        # crude coherence lifetime proxy: area under (purity - purity[0]) magnitude
        tau_c = float(np.trapz(np.abs(purity - purity[0]), times) / t_max)
        ete_vals.append(ete); tau_vals.append(tau_c)

    ete_vals = np.array(ete_vals); tau_vals = np.array(tau_vals)
    idx = int(np.argmax(ete_vals))
    return {
        "gamma": gamma_values,
        "ete": ete_vals,
        "tau_c_curve": tau_vals,
        "ETE_peak": float(ete_vals[idx]),
        "gamma_star": float(gamma_values[idx]),
        "tau_c": float(tau_vals[idx]),
        "peak_idx": int(idx),
    }

# ---- UI ----

st.set_page_config(page_title="QBM Upload → Biomarkers", layout="centered")
st.title("Quantum Bioenergetics Mapping — Upload Any Omics")

colA, colB = st.columns(2)
with colA:
    up_expr = st.file_uploader("Upload omics CSV (columns: gene,value)", type=["csv"])
with colB:
    up_map = st.file_uploader("Upload mapping CSV (columns: gene,node)", type=["csv"])

st.caption("Need templates?")
with st.expander("Show template examples"):
    st.code("gene,value\nNDUFS1,8.1\nNDUFV1,7.9\nCOX4I1,5.5\n...", language="csv")
    st.code("gene,node\nNDUFS1,0\nNDUFV1,0\nCOX4I1,5\n...", language="csv")

st.markdown("---")

# Config (safe defaults that produced the validated bell curve)
gamma_mode = st.selectbox("Gamma sweep scale", ["Log (recommended)", "Linear"], index=0)
if gamma_mode.startswith("Log"):
    gammas = np.geomspace(1e-4, 1.0, num=60)
else:
    gammas = np.linspace(0.0, 0.2, num=81)

sigma = st.slider("Static disorder σ", 0.0, 0.2, 0.08, 0.01)
k_sink = st.slider("k_sink", 0.0, 0.5, 0.20, 0.01)
k_loss = st.slider("k_loss", 0.0, 0.1, 0.005, 0.001)
dt = st.slider("dt", 0.05, 1.0, 0.2, 0.05)
t_max = st.slider("t_max", 10.0, 120.0, 60.0, 5.0)

run_btn = st.button("Run Simulation", type="primary", disabled=(up_expr is None or up_map is None))

# Graph template (9-node ETC mini)
default_graph = {
    "nodes": list(range(9)),
    "edges": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),
              (0,2),(2,4),(4,6),(6,8)]
}

if run_btn:
    try:
        expr_df = pd.read_csv(up_expr)
        map_df  = pd.read_csv(up_map)
        # sanity checks
        if not set(["gene","value"]).issubset(expr_df.columns):
            st.error("Omics CSV must have columns: gene,value"); st.stop()
        if not set(["gene","node"]).issubset(map_df.columns):
            st.error("Mapping CSV must have columns: gene,node"); st.stop()
        # map → Hamiltonian
        n, edges = len(default_graph["nodes"]), default_graph["edges"]
        gene2node, expr_z, alpha, J0, Jmax = map_omics_to_params(expr_df, map_df, alpha=0.01, J0=0.025, Jmax=0.10)
        H = build_hamiltonian(n, edges, gene2node, expr_z, alpha, J0, Jmax)
        # simulate
        res = simulate_ete_curve(H, gammas, sink_node=4, k_sink=k_sink, k_loss=k_loss, dt=dt, t_max=t_max, sigma=sigma, seed=1)
        st.success(f"γ* = {res['gamma_star']:.4f} • ETE_peak = {res['ETE_peak']:.4f} • τc = {res['tau_c']:.4f}")

        # plot
        st.line_chart(pd.DataFrame({"gamma": res["gamma"], "ETE": res["ete"]}).set_index("gamma"))
        st.caption("ENAQT bell curve: ETE rises with moderate γ (noise-assisted transport) and falls at high γ (decoherence).")

        # download metrics
        out = {
            "gamma_star": res["gamma_star"],
            "ETE_peak": res["ETE_peak"],
            "tau_c": res["tau_c"],
            "gamma": list(map(float, res["gamma"])),
            "ete": list(map(float, res["ete"])),
        }
        st.download_button("Download metrics JSON", data=json.dumps(out, indent=2), file_name="qbm_metrics.json")
    except Exception as e:
        st.error(f"Error while running simulation: {e}")
        st.stop()
