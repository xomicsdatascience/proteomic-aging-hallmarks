import json
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Re-run preprocessing & ComBat (fast, ~2s)
adata = ad.read_h5ad("data/raw (1).h5ad")
adata = adata[adata.obs["drug"] != "Drug"].copy()
X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.copy()
keep = (X == 0).mean(axis=0) <= 0.50
adata = adata[:, keep].copy()
X = X[:, keep]

with np.errstate(divide="ignore"):
    X_log = np.where(X > 0, np.log2(X), np.nan)

med = np.nanmedian(X_log, axis=1, keepdims=True)
X_norm = X_log - med + np.nanmedian(med)

rng = np.random.default_rng(42)
X_imp = X_norm.copy()
for i in range(X_imp.shape[0]):
    m = np.isnan(X_imp[i])
    if m.any():
        mu, sd = np.nanmean(X_imp[i]), np.nanstd(X_imp[i])
        X_imp[i, m] = rng.normal(mu - 1.8 * sd, 0.3 * sd, size=m.sum())

adata_pre = adata.copy()
adata_pre.X = X_imp.copy()

adata_post = adata.copy()
adata_post.X = X_imp.copy()
sc.pp.combat(adata_post, key="batch")

# --- FIGURE 1: Comprehensive ComBat Normalization Overview ---
fig, axes = plt.subplots(2, 2, figsize=(15, 11))

for ax, ad_obj, title in [
    (axes[0, 0], adata_pre, "A. All Samples PCA (n=1,231) — Before ComBat"),
    (axes[0, 1], adata_post, "B. All Samples PCA (n=1,231) — After ComBat"),
]:
    tmp = ad_obj.copy()
    sc.pp.pca(tmp, n_comps=2)
    coords = tmp.obsm["X_pca"]
    ve = tmp.uns["pca"]["variance_ratio"] * 100
    df_p = pd.DataFrame({"PC1": coords[:, 0], "PC2": coords[:, 1], "Batch": tmp.obs["batch"].values})
    sns.scatterplot(data=df_p, x="PC1", y="PC2", hue="Batch", s=18, alpha=0.65, ax=ax, palette="tab20", linewidth=0)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel(f"PC1 ({ve[0]:.1f}% var)")
    ax.set_ylabel(f"PC2 ({ve[1]:.1f}% var)")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7.5, frameon=False, markerscale=1.3)

# Panel C: Per-protein IQR across QC samples before vs after ComBat (Violin + Boxplot + Stripplot of 500 sampled proteins or all)
qc_pre = adata_pre[adata_pre.obs["sample_type"] == "QC"].X
qc_post = adata_post[adata_post.obs["sample_type"] == "QC"].X
sd_pre = np.std(qc_pre, axis=0)
sd_post = np.std(qc_post, axis=0)
df_sd = pd.DataFrame({
    "Stage": ["Before ComBat"] * len(sd_pre) + ["After ComBat"] * len(sd_post),
    "QC Protein SD (log2)": np.concatenate([sd_pre, sd_post]),
})
sns.violinplot(data=df_sd, x="Stage", y="QC Protein SD (log2)", inner=None, color="#e6eef8", ax=axes[1, 0])
sns.boxplot(data=df_sd, x="Stage", y="QC Protein SD (log2)", width=0.22, color="white", fliersize=0, ax=axes[1, 0])
# Subsample 300 proteins for stripplot readability
sub_idx = rng.choice(len(sd_pre), size=300, replace=False)
df_sd_sub = pd.DataFrame({
    "Stage": ["Before ComBat"] * 300 + ["After ComBat"] * 300,
    "QC Protein SD (log2)": np.concatenate([sd_pre[sub_idx], sd_post[sub_idx]]),
})
sns.stripplot(data=df_sd_sub, x="Stage", y="QC Protein SD (log2)", size=2.2, alpha=0.45, color="#1f77b4", ax=axes[1, 0])
axes[1, 0].set_title("C. Technical Variability Across 155 QC Samples", fontsize=12, fontweight="bold")
axes[1, 0].set_ylabel("Per-Protein SD across QC Samples (log2)")
axes[1, 0].set_xlabel("")
axes[1, 0].text(0.5, 0.92, f"Median SD: {np.median(sd_pre):.3f} -> {np.median(sd_post):.3f}",
                transform=axes[1, 0].transAxes, ha="center", fontsize=10, bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.8))

# Panel D: Sample median log2 intensity per batch after ComBat
post_meds = np.median(adata_post.X, axis=1)
df_med = pd.DataFrame({"Batch": adata_post.obs["batch"].values, "Median Log2": post_meds})
batch_order = sorted(df_med["Batch"].unique())
sns.boxplot(data=df_med, x="Batch", y="Median Log2", order=batch_order, color="lightgray", fliersize=0, ax=axes[1, 1])
sns.stripplot(data=df_med, x="Batch", y="Median Log2", order=batch_order, size=2.2, alpha=0.5, color="#2ca02c", ax=axes[1, 1])
axes[1, 1].set_title("D. Sample Median Log2 Intensity by Batch (After ComBat)", fontsize=12, fontweight="bold")
axes[1, 1].set_ylabel("Median Log2 Intensity")
axes[1, 1].set_xlabel("Batch")
axes[1, 1].tick_params(axis="x", rotation=45)

plt.tight_layout()
plt.savefig("results/fig1_combat_overview.png", dpi=250, bbox_inches="tight")
plt.close()

# --- FIGURE 4: Pathway Proteomic Coverage (Lollipop / Cleveland Dot Plot, No Barplots!) ---
df_gs = pd.read_csv("results/geneset_summary.csv")
df_gs["Pct"] = 100.0 * df_gs["Proteomics_Overlap"] / df_gs["Total_Genes"]
df_gs_plot = df_gs.iloc[::-1].reset_index(drop=True)  # reverse so 1 is at top

fig, ax = plt.subplots(figsize=(10.5, 6.5))
y_pos = np.arange(len(df_gs_plot))
ax.hlines(y=y_pos, xmin=df_gs_plot["Proteomics_Overlap"], xmax=df_gs_plot["Total_Genes"], color="#bbbbbb", linewidth=2, zorder=1)
ax.scatter(df_gs_plot["Total_Genes"], y_pos, color="#7f7f7f", s=70, label="Total Genes in Gene Set", zorder=2)
ax.scatter(df_gs_plot["Proteomics_Overlap"], y_pos, color="#1f77b4", s=95, label="Proteins Quantified in Proteomics", zorder=3)

for i, row in df_gs_plot.iterrows():
    ax.text(row["Total_Genes"] + 4, i, f"{int(row['Proteomics_Overlap'])}/{int(row['Total_Genes'])} ({row['Pct']:.0f}%)",
            va="center", fontsize=8.5, color="#333333")

ax.set_yticks(y_pos)
ax.set_yticklabels(df_gs_plot["GeneSet"], fontsize=9.5)
ax.set_xlabel("Number of Proteins / Genes", fontsize=11, fontweight="bold")
ax.set_title("Proteomic Coverage of the 12 Hallmarks of Aging & 2 Aggregate Aging Gene Sets", fontsize=12, fontweight="bold")
ax.set_xlim(0, df_gs_plot["Total_Genes"].max() + 55)
ax.legend(loc="lower right", frameon=True)
ax.grid(axis="x", linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("results/fig4_pathway_coverage.png", dpi=250, bbox_inches="tight")
plt.close()

# --- FIGURE 5: Deep-Dive into Kinase Inhibitors + DMSO (Violin + Stripplot) ---
with open("results/aging_hallmarks_genesets.json") as f:
    gene_sets = json.load(f)

ad_trt = adata_post[adata_post.obs["sample_type"] != "QC"].copy()
ctrl_X = ad_trt[ad_trt.obs["drug"] == "DMSO"].X
prots = ad_trt.var_names

# Compare Binimetinib, Mobocertinib (succinate), Erlotinib (Hydrochloride) vs DMSO
focus_drugs = ["DMSO", "Binimetinib", "Mobocertinib (succinate)", "Erlotinib (Hydrochloride)"]
oxphos_genes = [g for g in gene_sets["7. Mitochondrial Dysfunction (OXPHOS)"] if g in prots]
inflam_genes = [g for g in gene_sets["11. Chronic Inflammation"] if g in prots]

def get_top_sig_proteins(drug_name, gene_list, top_n=6, direction="up"):
    d_X = ad_trt[ad_trt.obs["drug"] == drug_name].X
    fc = d_X.mean(axis=0) - ctrl_X.mean(axis=0)
    _, pval = stats.ttest_ind(d_X, ctrl_X, axis=0, equal_var=False)
    _, qval, _, _ = multipletests(np.nan_to_num(pval, nan=1.0), method="fdr_bh")
    df_p = pd.DataFrame({"protein": prots, "log2fc": fc, "qval": qval})
    df_p = df_p[df_p["protein"].isin(gene_list)]
    if direction == "up":
        df_p = df_p[df_p["log2fc"] > 0].sort_values(["qval", "log2fc"], ascending=[True, False])
    else:
        df_p = df_p[df_p["log2fc"] < 0].sort_values(["qval", "log2fc"], ascending=[True, True])
    return df_p.head(top_n)

top_oxphos = get_top_sig_proteins("Binimetinib", oxphos_genes, top_n=6, direction="up")
top_inflam = get_top_sig_proteins("Binimetinib", inflam_genes, top_n=6, direction="down")
print("Top Binimetinib OXPHOS proteins (BH FDR):")
print(top_oxphos.to_string(index=False))
print("\nTop Binimetinib Inflammation proteins suppressed (BH FDR):")
print(top_inflam.to_string(index=False))

def build_long_df(prot_list):
    sub = ad_trt[ad_trt.obs["drug"].isin(focus_drugs), prot_list].copy()
    df_w = pd.DataFrame(sub.X, columns=prot_list, index=sub.obs_names)
    df_w["Drug"] = sub.obs["drug"].values
    df_l = df_w.melt(id_vars="Drug", var_name="Protein", value_name="Normalized Log2 Intensity")
    df_l["Drug"] = pd.Categorical(df_l["Drug"], categories=focus_drugs, ordered=True)
    return df_l

fig, axes = plt.subplots(2, 1, figsize=(13, 10))
palette = {"DMSO": "#7f7f7f", "Binimetinib": "#2ca02c", "Mobocertinib (succinate)": "#1f77b4", "Erlotinib (Hydrochloride)": "#d62728"}

for ax, df_top, title in [
    (axes[0], top_oxphos, "A. Mitochondrial Dysfunction (OXPHOS) — Top Upregulated Proteins by Binimetinib (BH FDR q < 0.05)"),
    (axes[1], top_inflam, "B. Chronic Inflammation — Top Suppressed Proteins by Binimetinib vs. Erlotinib (BH FDR q < 0.05)"),
]:
    plist = df_top["protein"].tolist()
    df_l = build_long_df(plist)
    sns.violinplot(data=df_l, x="Protein", y="Normalized Log2 Intensity", hue="Drug", palette=palette, inner=None, cut=0, alpha=0.35, ax=ax)
    sns.stripplot(data=df_l, x="Protein", y="Normalized Log2 Intensity", hue="Drug", palette=palette, dodge=True, size=4.2, alpha=0.9, jitter=0.12, ax=ax)
    # Clean duplicate legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:len(focus_drugs)], labels[:len(focus_drugs)], title="Treatment", loc="upper right", frameon=True, fontsize=8.5)
    ax.set_title(title, fontsize=11.5, fontweight="bold")
    ax.set_xlabel("Protein", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("ComBat-Normalized Log2 Intensity", fontsize=10.5, fontweight="bold")

plt.tight_layout()
plt.savefig("results/fig5_kinase_inhibitor_deepdive.png", dpi=250, bbox_inches="tight")
plt.close()
print("Saved fig1_combat_overview.png, fig4_pathway_coverage.png, and fig5_kinase_inhibitor_deepdive.png!")
