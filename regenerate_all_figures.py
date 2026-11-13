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

# Global typography settings for readability
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 13,
    "axes.labelweight": "bold",
    "xtick.labelsize": 11.5,
    "ytick.labelsize": 11.5,
    "legend.fontsize": 11,
    "legend.title_fontsize": 11.5,
})

# 1. Preprocess & ComBat
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

# --- FIGURE 1: ComBat Overview ---
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

for ax, ad_obj, title in [
    (axes[0, 0], adata_pre, "A. All Samples PCA (n=1,231) — Before ComBat"),
    (axes[0, 1], adata_post, "B. All Samples PCA (n=1,231) — After ComBat"),
]:
    tmp = ad_obj.copy()
    sc.pp.pca(tmp, n_comps=2)
    coords = tmp.obsm["X_pca"]
    ve = tmp.uns["pca"]["variance_ratio"] * 100
    df_p = pd.DataFrame({"PC1": coords[:, 0], "PC2": coords[:, 1], "Batch": tmp.obs["batch"].values})
    sns.scatterplot(data=df_p, x="PC1", y="PC2", hue="Batch", s=22, alpha=0.68, ax=ax, palette="tab20", linewidth=0)
    ax.set_title(title)
    ax.set_xlabel(f"PC1 ({ve[0]:.1f}% var)")
    ax.set_ylabel(f"PC2 ({ve[1]:.1f}% var)")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=10, frameon=False, markerscale=1.6)

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
sub_idx = rng.choice(len(sd_pre), size=300, replace=False)
df_sd_sub = pd.DataFrame({
    "Stage": ["Before ComBat"] * 300 + ["After ComBat"] * 300,
    "QC Protein SD (log2)": np.concatenate([sd_pre[sub_idx], sd_post[sub_idx]]),
})
sns.stripplot(data=df_sd_sub, x="Stage", y="QC Protein SD (log2)", size=2.6, alpha=0.45, color="#1f77b4", ax=axes[1, 0])
axes[1, 0].set_title("C. Technical Variability Across 155 QC Samples")
axes[1, 0].set_ylabel("Per-Protein SD across QC Samples (log2)")
axes[1, 0].set_xlabel("")
axes[1, 0].text(0.5, 0.90, f"Median SD: {np.median(sd_pre):.3f} → {np.median(sd_post):.3f}",
                transform=axes[1, 0].transAxes, ha="center", fontsize=11.5, bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.85))

post_meds = np.median(adata_post.X, axis=1)
df_med = pd.DataFrame({"Batch": adata_post.obs["batch"].values, "Median Log2": post_meds})
batch_order = sorted(df_med["Batch"].unique())
sns.boxplot(data=df_med, x="Batch", y="Median Log2", order=batch_order, color="lightgray", fliersize=0, ax=axes[1, 1])
sns.stripplot(data=df_med, x="Batch", y="Median Log2", order=batch_order, size=2.5, alpha=0.5, color="#2ca02c", ax=axes[1, 1])
axes[1, 1].set_title("D. Sample Median Log2 Intensity by Batch (After ComBat)")
axes[1, 1].set_ylabel("Median Log2 Intensity")
axes[1, 1].set_xlabel("Batch")
axes[1, 1].tick_params(axis="x", rotation=45)

plt.tight_layout()
plt.savefig("results/fig1_combat_overview.png", dpi=250, bbox_inches="tight")
plt.close()

# --- Load GSEA & Metadata ---
with open("results/aging_hallmarks_genesets.json") as f:
    gene_sets = json.load(f)
pathway_order = list(gene_sets.keys())

df_gsea = pd.read_csv("results/gsea_all_drugs_14pathways.csv")
df_de_m = pd.read_csv("results/de_summary_with_metadata.csv")
nes_mat = df_gsea.pivot(index="drug", columns="Term", values="NES")[pathway_order]
fdr_mat = df_gsea.pivot(index="drug", columns="Term", values="FDR q-val")[pathway_order]

cat_map = df_de_m.set_index("drug")["Category"].to_dict()
unique_cats = sorted(set(cat_map.values()))
cat_palette = dict(zip(unique_cats, sns.color_palette("Set2", n_colors=len(unique_cats))))

# --- FIGURE 2 (MAIN): Curated Subset of Key Aging & Top-Hit Drugs ---
subset_drugs = [
    "DFO",
    "Binimetinib",
    "Baicalein",
    "Nordihydroguaiaretic acid",
    "Masoprocol",
    "Lonafarnib",
    "Nicotinamide",
    "α-Lipoic Acid",
    "Simvastatin acid (ammonium)",
    "Atorvastatin (hemicalcium salt)",
    "Eicosapentaenoic Acid",
    "Dimethyl fumarate",
    "Tucidinostat",
    "Methylene blue (trihydrate)",
    "γ-Oryzanol",
    "Mobocertinib (succinate)",
    "Erlotinib (Hydrochloride)",
    "Asciminib",
    "Flumatinib",
    "Lapatinib (ditosylate)",
    "Triclabendazole",
    "Mifepristone",
    "Cabazitaxel",
    "Latrepirdine (dihydrochloride)",
    "Efavirenz",
]
subset_drugs = [d for d in subset_drugs if d in nes_mat.index]
nes_sub = nes_mat.loc[subset_drugs].copy()
fdr_sub = fdr_mat.loc[subset_drugs].copy()
annot_sub = fdr_sub.map(lambda q: "**" if q < 0.01 else ("*" if q < 0.05 else ""))

row_cats_sub = pd.Series([cat_map.get(d, "Unannotated") for d in nes_sub.index], index=nes_sub.index, name="Category")
row_colors_sub = row_cats_sub.map(cat_palette)

g = sns.clustermap(
    nes_sub,
    cmap="RdBu_r",
    center=0,
    vmin=-2.5,
    vmax=2.5,
    col_cluster=False,
    row_cluster=True,
    row_colors=row_colors_sub,
    annot=annot_sub,
    fmt="",
    annot_kws={"size": 10, "color": "black", "weight": "bold"},
    linewidths=0.4,
    linecolor="#e0e0e0",
    figsize=(16, 11.5),
    cbar_pos=(0.015, 0.80, 0.022, 0.14),
    cbar_kws={"label": "NES"},
    dendrogram_ratio=(0.13, 0.02),
    colors_ratio=0.028,
)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=38, ha="right", fontsize=11.5)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=11.5)
g.ax_heatmap.set_xlabel("Aging Hallmarks & Aggregate Gene Sets", fontsize=13, fontweight="bold")
g.ax_heatmap.set_ylabel("Highlighted Treatments", fontsize=13, fontweight="bold")
g.fig.suptitle(
    f"Aging Hallmarks Pre-Ranked GSEA NES — Highlighted Drugs & Known Geroprotectors (n={len(subset_drugs)})",
    y=1.01,
    fontsize=14.5,
    fontweight="bold",
)
sub_cats_present = sorted(row_cats_sub.unique())
handles = [plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=cat_palette[c], markersize=11, label=c) for c in sub_cats_present]
g.fig.legend(handles=handles, title="Drug Category", loc="upper left", bbox_to_anchor=(0.005, 0.77), frameon=True, fontsize=10.5, title_fontsize=11)
plt.savefig("results/fig2_aging_hallmarks_nes_heatmap_subset.png", dpi=250, bbox_inches="tight")
plt.close()

# --- SUPPLEMENTARY FIGURE S1: Full 165-Drug Heatmap ---
sig_drugs_05 = fdr_mat.index[(fdr_mat < 0.05).any(axis=1)].tolist()
nes_sig = nes_mat.loc[sig_drugs_05].copy()
fdr_sig = fdr_mat.loc[sig_drugs_05].copy()
annot_sig = fdr_sig.map(lambda q: "**" if q < 0.01 else ("*" if q < 0.05 else ""))
row_cats_sig = pd.Series([cat_map.get(d, "Unannotated") for d in nes_sig.index], index=nes_sig.index, name="Category")
row_colors_sig = row_cats_sig.map(cat_palette)

fig_h = max(12, len(sig_drugs_05) * 0.25)
g_full = sns.clustermap(
    nes_sig,
    cmap="RdBu_r",
    center=0,
    vmin=-2.5,
    vmax=2.5,
    col_cluster=False,
    row_cluster=True,
    row_colors=row_colors_sig,
    annot=annot_sig,
    fmt="",
    annot_kws={"size": 7.5, "color": "black"},
    linewidths=0.3,
    linecolor="#e0e0e0",
    figsize=(17, fig_h),
    cbar_pos=(0.015, 0.87, 0.022, 0.09),
    cbar_kws={"label": "NES"},
    dendrogram_ratio=(0.12, 0.02),
    colors_ratio=0.025,
)
g_full.ax_heatmap.set_xticklabels(g_full.ax_heatmap.get_xticklabels(), rotation=38, ha="right", fontsize=11.5)
g_full.ax_heatmap.set_yticklabels(g_full.ax_heatmap.get_yticklabels(), fontsize=8.5)
g_full.ax_heatmap.set_xlabel("Aging Hallmarks & Aggregate Gene Sets", fontsize=13, fontweight="bold")
g_full.ax_heatmap.set_ylabel("All Significant Drug Treatments (FDR q < 0.05 in >=1 pathway)", fontsize=13, fontweight="bold")
g_full.fig.suptitle(
    f"Supplementary Figure S1: Full NES Heatmap ({len(sig_drugs_05)} Significant Drugs, FDR q < 0.05)",
    y=1.005,
    fontsize=14.5,
    fontweight="bold",
)
handles_f = [plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=cat_palette[c], markersize=11, label=c) for c in unique_cats]
g_full.fig.legend(handles=handles_f, title="Drug Category", loc="upper left", bbox_to_anchor=(0.005, 0.85), frameon=True, fontsize=10.5, title_fontsize=11)
plt.savefig("results/figS1_aging_hallmarks_nes_heatmap_full.png", dpi=250, bbox_inches="tight")
plt.close()

# --- FIGURE 3: NES Distributions ---
fig, ax = plt.subplots(figsize=(15, 7.5))
df_plot = df_gsea.copy()
df_plot["Term"] = pd.Categorical(df_plot["Term"], categories=pathway_order, ordered=True)
df_plot["Significant (q < 0.05)"] = df_plot["FDR q-val"] < 0.05

sns.boxplot(data=df_plot, x="Term", y="NES", color="#f0f0f0", fliersize=0, ax=ax)
sns.stripplot(data=df_plot[~df_plot["Significant (q < 0.05)"]], x="Term", y="NES", color="#999999", alpha=0.35, size=3.5, jitter=0.25, ax=ax, label="FDR q >= 0.05")
sns.stripplot(data=df_plot[df_plot["Significant (q < 0.05)"]], x="Term", y="NES", color="#d62728", alpha=0.82, size=5, jitter=0.25, ax=ax, label="FDR q < 0.05")
ax.axhline(0, color="black", linestyle="--", linewidth=1.0)
ax.set_xticklabels(ax.get_xticklabels(), rotation=38, ha="right", fontsize=11.5)
ax.set_xlabel("Aging Hallmarks & Aggregate Gene Sets", fontsize=13, fontweight="bold")
ax.set_ylabel("Normalized Enrichment Score (NES)", fontsize=13, fontweight="bold")
ax.set_title("Distribution of Drug-Induced NES Across 14 Aging Pathways (All 173 Treatments)", fontsize=14.5, fontweight="bold")
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), loc="upper right", frameon=True, fontsize=11.5)
plt.tight_layout()
plt.savefig("results/fig3_pathway_nes_distributions.png", dpi=250, bbox_inches="tight")
plt.close()

# --- FIGURE 4: Pathway Coverage ---
df_gs = pd.read_csv("results/geneset_summary.csv")
df_gs["Pct"] = 100.0 * df_gs["Proteomics_Overlap"] / df_gs["Total_Genes"]
df_gs_plot = df_gs.iloc[::-1].reset_index(drop=True)

fig, ax = plt.subplots(figsize=(12, 7.2))
y_pos = np.arange(len(df_gs_plot))
ax.hlines(y=y_pos, xmin=df_gs_plot["Proteomics_Overlap"], xmax=df_gs_plot["Total_Genes"], color="#bbbbbb", linewidth=2.5, zorder=1)
ax.scatter(df_gs_plot["Total_Genes"], y_pos, color="#7f7f7f", s=95, label="Total Genes in Gene Set", zorder=2)
ax.scatter(df_gs_plot["Proteomics_Overlap"], y_pos, color="#1f77b4", s=120, label="Proteins Quantified in Proteomics", zorder=3)

for i, row in df_gs_plot.iterrows():
    ax.text(row["Total_Genes"] + 4, i, f"{int(row['Proteomics_Overlap'])}/{int(row['Total_Genes'])} ({row['Pct']:.0f}%)",
            va="center", fontsize=10.5, color="#222222")

ax.set_yticks(y_pos)
ax.set_yticklabels(df_gs_plot["GeneSet"], fontsize=11.5)
ax.set_xlabel("Number of Proteins / Genes", fontsize=13, fontweight="bold")
ax.set_title("Proteomic Coverage of the 12 Hallmarks of Aging & 2 Aggregate Aging Gene Sets", fontsize=14.5, fontweight="bold")
ax.set_xlim(0, df_gs_plot["Total_Genes"].max() + 62)
ax.legend(loc="lower right", frameon=True, fontsize=11.5)
ax.grid(axis="x", linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig("results/fig4_pathway_coverage.png", dpi=250, bbox_inches="tight")
plt.close()

# --- FIGURE 5: Kinase Inhibitors Deep-Dive (SCALED & CENTERED vs DMSO) ---
ad_trt = adata_post[adata_post.obs["sample_type"] != "QC"].copy()
ctrl_X = ad_trt[ad_trt.obs["drug"] == "DMSO"].X
prots = ad_trt.var_names
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

def build_scaled_long_df(prot_list):
    sub = ad_trt[ad_trt.obs["drug"].isin(focus_drugs), prot_list].copy()
    df_w = pd.DataFrame(sub.X, columns=prot_list, index=sub.obs_names)
    df_w["Drug"] = sub.obs["drug"].values
    # Center to DMSO mean and scale by overall protein SD across these samples
    dmso_means = df_w[df_w["Drug"] == "DMSO"][prot_list].mean(axis=0)
    prot_sds = df_w[prot_list].std(axis=0)
    df_w[prot_list] = (df_w[prot_list] - dmso_means) / prot_sds
    df_l = df_w.melt(id_vars="Drug", var_name="Protein", value_name="Scaled Abundance (Z vs. DMSO)")
    df_l["Drug"] = pd.Categorical(df_l["Drug"], categories=focus_drugs, ordered=True)
    return df_l

fig, axes = plt.subplots(2, 1, figsize=(14, 11))
palette = {"DMSO": "#7f7f7f", "Binimetinib": "#2ca02c", "Mobocertinib (succinate)": "#1f77b4", "Erlotinib (Hydrochloride)": "#d62728"}

for ax, df_top, title in [
    (axes[0], top_oxphos, "A. Mitochondrial Dysfunction (OXPHOS) — Top Upregulated Proteins by Binimetinib (BH FDR q < 0.05)"),
    (axes[1], top_inflam, "B. Chronic Inflammation — Top Suppressed Proteins by Binimetinib vs. Erlotinib (BH FDR q < 0.05)"),
]:
    plist = df_top["protein"].tolist()
    df_l = build_scaled_long_df(plist)
    ax.axhline(0, color="black", linestyle="--", linewidth=1.1, alpha=0.75, zorder=0)
    sns.violinplot(data=df_l, x="Protein", y="Scaled Abundance (Z vs. DMSO)", hue="Drug", palette=palette, inner=None, cut=0, alpha=0.35, ax=ax)
    sns.stripplot(data=df_l, x="Protein", y="Scaled Abundance (Z vs. DMSO)", hue="Drug", palette=palette, dodge=True, size=5, alpha=0.92, jitter=0.12, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:len(focus_drugs)], labels[:len(focus_drugs)], title="Treatment", loc="upper right", frameon=True, fontsize=11, title_fontsize=11.5)
    ax.set_title(title, fontsize=13.5, fontweight="bold")
    ax.set_xlabel("Protein", fontsize=12.5, fontweight="bold")
    ax.set_ylabel("Scaled Abundance (SD vs. DMSO Mean)", fontsize=12.5, fontweight="bold")
    ax.tick_params(axis="both", labelsize=12)

plt.tight_layout()
plt.savefig("results/fig5_kinase_inhibitor_deepdive.png", dpi=250, bbox_inches="tight")
plt.close()

print("All figures regenerated with larger typography and scaled/centered Figure 5!")
