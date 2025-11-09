
#!/usr/bin/env python3
import argparse, json, yaml, pandas as pd, numpy as np, pathlib, sys
from qmnet import OmicsMapper, BatchSimulator
from QB_1_fixed_enaqt_physics import FixedQuantumIntegrator, FixedMetricsCalculator
import json as _json

def load_graph(path):
    with open(path, 'r') as f:
        gspec = _json.load(f)
    import networkx as nx
    G = nx.Graph()
    for n in gspec["nodes"]:
        G.add_node(int(n))
    for i,j in gspec["edges"]:
        G.add_edge(int(i), int(j))
    return G

def load_mapping(path):
    import pandas as pd
    df = pd.read_csv(path)
    return {str(r["gene"]).strip(): int(r["node"]) for _, r in df.iterrows()}

def zscore_dict(d):
    import numpy as np
    vals = np.array(list(d.values()), dtype=float)
    sd = vals.std() or 1.0
    mu = vals.mean()
    return {k:(v - mu)/sd for k,v in d.items()}

def run_single(expr_path, graph_json, mapping_csv, config_yaml, out_json):
    cfg = yaml.safe_load(open(config_yaml, "r"))
    G = load_graph(graph_json)
    gene2node = load_mapping(mapping_csv)

    df = pd.read_csv(expr_path)
    expr = zscore_dict({str(r["gene"]).strip(): float(r["value"]) for _, r in df.iterrows()})

    mapper = OmicsMapper(alpha=cfg['hamiltonian']['alpha'], J_0=cfg['hamiltonian']['J_0'], J_max=cfg['hamiltonian']['J_max'])
    mapper.set_gene_mapping(gene2node)
    sim = BatchSimulator(cfg, mapper)
    res = sim.simulate_single_sample(expr, G, gamma_sweep=True)
    pathlib.Path(out_json).write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))

def run_cohort(cohort_csv, graph_json, mapping_csv, config_yaml, out_parquet):
    import pyarrow as pa, pyarrow.parquet as pq
    cfg = yaml.safe_load(open(config_yaml, "r"))
    G = load_graph(graph_json)
    gene2node = load_mapping(mapping_csv)
    mapper = OmicsMapper(alpha=cfg['hamiltonian']['alpha'], J_0=cfg['hamiltonian']['J_0'], J_max=cfg['hamiltonian']['J_max'])
    mapper.set_gene_mapping(gene2node)
    sim = BatchSimulator(cfg, mapper)

    dfc = pd.read_csv(cohort_csv)
    rows = []
    for _, r in dfc.iterrows():
        expr_df = pd.read_csv(r["filepath"])
        expr = zscore_dict({str(rr["gene"]).strip(): float(rr["value"]) for _, rr in expr_df.iterrows()})
        res = sim.simulate_single_sample(expr, G, gamma_sweep=True)
        rows.append({"sample_id": r["sample_id"], "group": r["group"], **res})
    out_df = pd.DataFrame(rows)
    import pyarrow as pa, pyarrow.parquet as pq
    table = pa.Table.from_pandas(out_df)
    pq.write_table(table, out_parquet)
    print(out_df.to_string(index=False))

def main():
    p = argparse.ArgumentParser(description="Quantum Energy Mapping Diagnostics (QEMD) CLI")
    sub = p.add_subparsers(dest="cmd", required=True)

    s1 = sub.add_parser("simulate", help="Run a single sample and output JSON metrics")
    s1.add_argument("--expr", required=True, help="CSV with columns: gene,value")
    s1.add_argument("--graph", default="./data/templates/graph.json")
    s1.add_argument("--mapping", default="./data/templates/mapping.csv")
    s1.add_argument("--config", default="./config/default_params.yaml")
    s1.add_argument("--out", default="./outputs/metrics.json")

    s2 = sub.add_parser("cohort", help="Run a CSV cohort list and emit a metrics parquet")
    s2.add_argument("--cohort", default="./data/templates/cohort.csv")
    s2.add_argument("--graph", default="./data/templates/graph.json")
    s2.add_argument("--mapping", default="./data/templates/mapping.csv")
    s2.add_argument("--config", default="./config/default_params.yaml")
    s2.add_argument("--out", default="./outputs/cohort_metrics.parquet")

    args = p.parse_args()
    if args.cmd == "simulate":
        run_single(args.expr, args.graph, args.mapping, args.config, args.out)
    else:
        run_cohort(args.cohort, args.graph, args.mapping, args.config, args.out)

if __name__ == "__main__":
    main()
