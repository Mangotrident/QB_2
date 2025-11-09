
import pandas as pd, numpy as np, pyarrow.parquet as pq, matplotlib.pyplot as plt, seaborn as sns
from scipy import stats
from pathlib import Path

def load_parquet(p):
    return pq.read_table(p).to_pandas()

def violin_box(df, feature, out_png):
    plt.figure(figsize=(6,4))
    sns.violinplot(data=df, x="group", y=feature, inner=None)
    sns.boxplot(data=df, x="group", y=feature, width=0.18, showcaps=True, boxprops={'zorder':2})
    plt.title(f"{feature} by Group")
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()

def effect_and_p(df, feature):
    g = df.groupby("group")[feature]
    groups = list(g.groups.keys())
    if len(groups) != 2:
        return None
    a, b = g.get_group(groups[0]).values, g.get_group(groups[1]).values
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = (np.mean(a)-np.mean(b))/np.sqrt(0.5*(np.var(a, ddof=1)+np.var(b, ddof=1)))
    return {"groups": groups, "p": float(u.pvalue), "cohens_d": float(d)}

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", default="./outputs/cohort_metrics.parquet")
    ap.add_argument("--outdir", default="./outputs")
    args = ap.parse_args()

    df = load_parquet(args.parquet)
    Path(args.outdir).mkdir(exist_ok=True, parents=True)

    for feat in ["ETE_peak","gamma_star","tau_c"]:
        violin_box(df, feat, Path(args.outdir)/f"violin_{feat}.png")
        stats_row = effect_and_p(df, feat)
        if stats_row:
            print(feat, stats_row)
