import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs("reanalysis/A08/figures", exist_ok=True)
os.makedirs("reanalysis/A08/source_tables", exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelsize": 11,
    "axes.labelweight": "bold",
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 9.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
})

def save_fig(fig, name):
    fig.savefig(f"reanalysis/A08/figures/{name}.png", bbox_inches="tight")
    fig.savefig(f"reanalysis/A08/figures/{name}.svg", bbox_inches="tight")
    plt.close(fig)

# Load core A01-A07 tables (160 library compounds + DFO across 12 plates)
sm = pd.read_csv("reanalysis/A01/sample_manifest.tsv", sep="\t")
sm_12 = sm[sm["plate"] != "RedMix"].copy()
sa = pd.read_csv("reanalysis/A01/sample_accounting.tsv", sep="\t")
sa_12 = sa[sa["plate"] != "RedMix"].copy()
qc = pd.read_csv("reanalysis/A04/qc_metrics.tsv", sep="\t")
qc_12 = qc[qc["plate"] != "RedMix"].copy()
rep = pd.read_csv("reanalysis/A04/replication_metrics.tsv", sep="\t")
null_cal = pd.read_csv("reanalysis/A04/null_calibration.tsv", sep="\t")
pman = pd.read_csv("reanalysis/A03/program_manifest.tsv", sep="\t")
pov = pd.read_csv("reanalysis/A03/program_overlap.tsv", sep="\t")
prog_df = pd.read_csv("reanalysis/A02/program_effects.tsv", sep="\t")
spec_df = pd.read_csv("reanalysis/A06/candidate_specificity.tsv", sep="\t")
bench_df = pd.read_csv("reanalysis/A05/drug_reference_comparisons.tsv", sep="\t")
sens_df = pd.read_csv("reanalysis/A04/sensitivity_effects.tsv", sep="\t")
comp_sens = pd.read_csv("reanalysis/A06/composition_sensitivity.tsv", sep="\t")

aging_ids = [
    "PROG_DNA_REPAIR", "PROG_TELOMERE_MAINT", "PROG_CHROMATIN_ORG",
    "PROG_PROTEOSTASIS_UPR", "PROG_PROTEASOME", "PROG_MACROAUTOPHAGY",
    "PROG_MTORC1_SIGNALING", "PROG_OXPHOS", "PROG_SENESCENCE_CORE",
    "PROG_CELL_CELL_COMM", "PROG_INFLAMMATORY", "PROG_INNATE_TLR",
    "PRIOR_DEMAGALHAES_UP", "PRIOR_DEMAGALHAES_DN",
]

# ==============================================================================
# MAIN FIGURE 1: Number of proteins per hallmark/program vs all genes in set
# ==============================================================================
pman_14 = pman[pman["program_id"].isin(aging_ids)].set_index("program_id").loc[aging_ids].reset_index()
fig, axes = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw={"width_ratios": [1.15, 1.0]})

y_pos = np.arange(len(pman_14))
axes[0].hlines(y=y_pos, xmin=pman_14["measured_retained_7871"], xmax=pman_14["total_genes"], color="#bdbdbd", linewidth=2.5)
axes[0].scatter(pman_14["total_genes"], y_pos, color="#7f7f7f", s=90, label="Total Genes in Set", zorder=3)
axes[0].scatter(pman_14["measured_retained_7871"], y_pos, color="#1f77b4", s=90, label="Measured Proteins (7,871 retained)", zorder=4)
for i, row in pman_14.iterrows():
    axes[0].text(row["total_genes"] + 5, i, f"{int(row['measured_retained_7871'])}/{int(row['total_genes'])} ({row['coverage_pct_retained']:.0f}%)", va="center", fontsize=8.5)
axes[0].set_yticks(y_pos)
axes[0].set_yticklabels(pman_14["program_id"])
axes[0].invert_yaxis()
axes[0].set_xlim(0, 390)
axes[0].set_xlabel("Gene / Protein Count")
axes[0].set_title("A. Measured Proteins vs. Total Annotated Genes per Program")
axes[0].legend(loc="lower right", frameon=True)

pov_14 = pd.DataFrame(0.0, index=aging_ids, columns=aging_ids)
np.fill_diagonal(pov_14.values, 1.0)
for _, r in pov.iterrows():
    if r["program_1"] in aging_ids and r["program_2"] in aging_ids:
        pov_14.loc[r["program_1"], r["program_2"]] = r["jaccard"]
        pov_14.loc[r["program_2"], r["program_1"]] = r["jaccard"]
sns.heatmap(pov_14, cmap="Blues", vmin=0, vmax=0.35, annot=True, fmt=".2f", annot_kws={"size": 7.5},
            ax=axes[1], cbar_kws={"label": "Jaccard Overlap"})
axes[1].set_title("B. Pairwise Protein Membership Overlap (Jaccard)")
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha="right", fontsize=8.5)
axes[1].set_yticklabels(axes[1].get_yticklabels(), fontsize=8.5)

plt.tight_layout()
save_fig(fig, "fig1_program_coverage")
pman_14.to_csv("reanalysis/A08/source_tables/fig1_program_coverage.tsv", sep="\t", index=False)

# ==============================================================================
# MAIN FIGURE 2: Subset heatmap of highlighted compounds across 14 programs
# ==============================================================================
prim = prog_df[prog_df["program_family"] == "primary_aging"].copy()
piv_t = prim.pivot(index="compound_id", columns="program_id", values="camera_t")[aging_ids]
piv_nes = prim.pivot(index="compound_id", columns="program_id", values="nes_gsea")[aging_ids]
piv_q = prim.pivot(index="compound_id", columns="program_id", values="q_BH_global_camera")[aging_ids]
piv_q_gsea = prim.pivot(index="compound_id", columns="program_id", values="q_BH_global_gsea")[aging_ids]

highlight_drugs = [
    "Binimetinib", "Baicalein", "Nordihydroguaiaretic acid", "Lonafarnib",
    "Nicotinamide", "Simvastatin acid (ammonium)", "Atorvastatin (hemicalcium salt)",
    "Eicosapentaenoic Acid", "Dimethyl fumarate", "Tucidinostat",
    "Methylene blue (trihydrate)", "Mobocertinib (succinate)", "Erlotinib (Hydrochloride)",
    "Asciminib", "Flumatinib", "Lapatinib (ditosylate)", "Triclabendazole",
    "Mifepristone", "6-Mercaptopurine", "DFO",
]
highlight_drugs = [d for d in highlight_drugs if d in piv_t.index]
sub_t = piv_t.loc[highlight_drugs]
sub_q = piv_q.loc[highlight_drugs]
annot_sub = sub_q.map(lambda q: "**" if q < 0.01 else ("*" if q < 0.05 else ""))

fig, ax = plt.subplots(figsize=(13.5, 8.2))
sns.heatmap(sub_t, cmap="RdBu_r", center=0, vmin=-6, vmax=6, annot=annot_sub, fmt="", annot_kws={"size": 10},
            linewidths=0.4, linecolor="white", ax=ax, cbar_kws={"label": "Correlation-Aware CAMERA t-statistic"})
ax.set_title("Highlighted Compounds Across 14 Aging-Associated Programs (* Global BH q < 0.05, ** q < 0.01)")
ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right", fontsize=9.5)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=9.5)
ax.set_xlabel("")
ax.set_ylabel("Compound")
plt.tight_layout()
save_fig(fig, "fig2_highlighted_heatmap")

# Main Supplemental Table S1 (all 160 library compounds + DFO)
supp_t1 = prim.sort_values(["compound_id", "program_id"]).copy()
supp_t1.to_csv("reanalysis/A08/source_tables/Supplementary_Table_S1_all_compounds.tsv", sep="\t", index=False)

# ==============================================================================
# MAIN FIGURE 3: Boxplot + stripplot of drug response distributions in NES
# ==============================================================================
prim["sig_gsea_q05"] = prim["q_BH_global_gsea"] < 0.05
fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2), gridspec_kw={"width_ratios": [1.3, 1.0]})

sns.boxplot(data=prim, x="program_id", y="nes_gsea", order=aging_ids, color="#f0f0f0", fliersize=0, ax=axes[0])
sns.stripplot(data=prim[~prim["sig_gsea_q05"]], x="program_id", y="nes_gsea", order=aging_ids, color="#999999", alpha=0.35, size=2.5, ax=axes[0])
sns.stripplot(data=prim[prim["sig_gsea_q05"]], x="program_id", y="nes_gsea", order=aging_ids, color="#d62728", alpha=0.75, size=3.5, ax=axes[0])
axes[0].axhline(0, color="black", linestyle="--", linewidth=0.9)
axes[0].set_xticks(range(len(aging_ids)))
axes[0].set_xticklabels(aging_ids, rotation=45, ha="right", fontsize=8.5)
axes[0].set_title("A. Distribution of GSEA NES Across 160 Compounds (+ DFO)")
axes[0].set_ylabel("Normalized Enrichment Score (NES)")
axes[0].set_xlabel("")

sc_p = axes[1].scatter(spec_df["NUIS_E2F_PROLIF"], spec_df["PROG_INFLAMMATORY"], c=spec_df["PROG_OXPHOS"],
                       cmap="RdBu_r", vmin=-6, vmax=6, s=48, alpha=0.85, edgecolors="k", linewidths=0.4)
plt.colorbar(sc_p, ax=axes[1], label="PROG_OXPHOS (CAMERA t)")
offsets_3b = {"Binimetinib": (-12, 8), "Baicalein": (8, -10), "DFO": (6, 4), "Erlotinib (Hydrochloride)": (6, 4)}
for d_lbl in ["Binimetinib", "Baicalein", "DFO", "Erlotinib (Hydrochloride)"]:
    r_d = spec_df[spec_df["compound_id"] == d_lbl].iloc[0]
    axes[1].annotate(d_lbl, (r_d["NUIS_E2F_PROLIF"], r_d["PROG_INFLAMMATORY"]), xytext=offsets_3b[d_lbl], textcoords="offset points", fontsize=8.5, fontweight="bold")
axes[1].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes[1].axvline(0, color="gray", linestyle="--", linewidth=0.8)
axes[1].set_xlabel("Proliferation Proxy (NUIS_E2F_PROLIF CAMERA t)")
axes[1].set_ylabel("Intracellular Inflammation (PROG_INFLAMMATORY CAMERA t)")
axes[1].set_title("B. Specificity vs. Proliferation & OXPHOS Axes")

plt.tight_layout()
save_fig(fig, "fig3_nes_distributions")
piv_nes.to_csv("reanalysis/A08/source_tables/fig3_nes_all_compounds.tsv", sep="\t")

# ==============================================================================
# MAIN FIGURE 4: Kinase inhibitor panel (Zero-centered & unit-variance scaled)
# ==============================================================================
prot_map = pd.read_csv("reanalysis/A01/protein_mapping.tsv", sep="\t")
ret_prots = prot_map.loc[prot_map["retained_50pct_filter"], "feature_id"].tolist()
adata = ad.read_h5ad("data/raw (1).h5ad")[:, ret_prots].copy()
X_raw = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.copy()
with np.errstate(divide="ignore"):
    X_log = np.where(X_raw > 0, np.log2(X_raw), np.nan)
meds = np.nanmedian(X_log, axis=1, keepdims=True)
X_norm = X_log - meds + np.nanmedian(meds)

kinase_drugs = [
    "DMSO",
    "Binimetinib",
    "Erlotinib (Hydrochloride)",
    "Mobocertinib (succinate)",
    "Asciminib",
    "Lapatinib (ditosylate)",
]
ox_sentinels = [g for g in ["NDUFS1", "IDH2", "ACO2", "SDHA", "COX4I1", "ATP5F1A"] if g in ret_prots]
inf_sentinels = [g for g in ["ICAM1", "LDLR", "NFKB1", "RELA", "STAT3", "ATP2B1"] if g in ret_prots]

pal_k = {
    "DMSO": "#7f7f7f",
    "Binimetinib": "#2ca02c",
    "Erlotinib (Hydrochloride)": "#d62728",
    "Mobocertinib (succinate)": "#1f77b4",
    "Asciminib": "#9467bd",
    "Lapatinib (ditosylate)": "#ff7f0e",
}

m_k = sm["compound_id"].isin(kinase_drugs).values & (~sm["excluded"].values) & (sm["plate"].values != "RedMix")
fig, axes = plt.subplots(2, 1, figsize=(14.5, 10.5))

for ax, plist, title in [
    (axes[0], ox_sentinels, "A. Kinase Inhibitor Panel: Respiratory / OXPHOS Sentinels (Z-Scored Abundance)"),
    (axes[1], inf_sentinels, "B. Kinase Inhibitor Panel: Inflammatory Sentinels (Z-Scored Abundance)"),
]:
    p_idx = [ret_prots.index(g) for g in plist]
    sub_mat = X_norm[np.ix_(m_k, p_idx)]
    # Zero-center and unit-variance scale per protein across the plotted samples
    mu = np.nanmean(sub_mat, axis=0, keepdims=True)
    sd = np.nanstd(sub_mat, axis=0, ddof=1, keepdims=True)
    sub_z = (sub_mat - mu) / sd
    df_w = pd.DataFrame(sub_z, columns=plist)
    df_w["Drug"] = sm.loc[m_k, "compound_id"].values
    df_l = df_w.melt(id_vars="Drug", var_name="Protein", value_name="Scaled Abundance (Z-score)")
    df_l["Drug"] = pd.Categorical(df_l["Drug"], categories=kinase_drugs, ordered=True)
    sns.violinplot(data=df_l, x="Protein", y="Scaled Abundance (Z-score)", hue="Drug", palette=pal_k, inner=None, cut=0, alpha=0.32, ax=ax)
    sns.stripplot(data=df_l, x="Protein", y="Scaled Abundance (Z-score)", hue="Drug", palette=pal_k, dodge=True, size=4.0, alpha=0.9, jitter=0.1, ax=ax)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[:len(kinase_drugs)], l[:len(kinase_drugs)], title="Treatment (n=6 plates; DMSO n=23)", loc="lower right", frameon=True, ncol=3, fontsize=8.5)
    ax.set_title(title)
    ax.set_ylabel("Scaled Abundance (Z-score)")

plt.tight_layout()
save_fig(fig, "fig4_kinase_inhibitor_panel")

# ==============================================================================
# SUPPLEMENTARY FIGURE S1: QC, Design, Cross-Plate Replication & Null Calibration
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))

sa_melt = sa_12.melt(id_vars="plate", value_vars=["analytical_qc", "library_compound", "positive_assay_control", "unlabeled_excluded", "vehicle_control"],
                     var_name="Role", value_name="Count")
sa_melt = sa_melt[sa_melt["Count"] > 0]
sns.scatterplot(data=sa_melt, x="plate", y="Count", hue="Role", style="Role", s=110, ax=axes[0, 0], palette="Set2")
axes[0, 0].set_title("A. Sample Composition Across 12 Assay Plates")
axes[0, 0].set_xlabel("Plate")
axes[0, 0].set_ylabel("Sample Count")
axes[0, 0].tick_params(axis="x", rotation=45)

sns.scatterplot(data=qc_12, x="plate", y="qc_median_protein_sd", s=110, color="#1f77b4", ax=axes[0, 1])
axes[0, 1].axhline(qc_12["qc_median_protein_sd"].median(), color="gray", linestyle="--", linewidth=1)
axes[0, 1].set_title("B. Median Per-Protein SD Across QC Samples by Plate")
axes[0, 1].set_xlabel("Plate")
axes[0, 1].set_ylabel("Median QC Protein SD (log2)")
axes[0, 1].tick_params(axis="x", rotation=45)

rep_s = rep.sort_values("mean_lopo_pearson_r", ascending=True)
y_r = np.arange(len(rep_s))
axes[1, 0].hlines(y=y_r, xmin=rep_s["min_lopo_pearson_r"], xmax=rep_s["max_lopo_pearson_r"], color="#bbbbbb", linewidth=3)
axes[1, 0].scatter(rep_s["mean_lopo_pearson_r"], y_r, color="#2ca02c", s=100, zorder=3)
axes[1, 0].set_yticks(y_r)
axes[1, 0].set_yticklabels(rep_s["compound_id"])
axes[1, 0].set_xlim(0, 1.0)
axes[1, 0].set_xlabel("Leave-One-Plate-Out Pearson r (log2FC)")
axes[1, 0].set_title("C. Cross-Plate Replication of Proteome Perturbations")
axes[1, 0].grid(axis="x", linestyle="--", alpha=0.4)

null_s = null_cal.sort_values("null_camera_t", key=np.abs, ascending=True)
y_n = np.arange(len(null_s))
axes[1, 1].axvline(0, color="black", linestyle="--", linewidth=1)
axes[1, 1].scatter(null_s["null_camera_t"], y_n, color="#7f7f7f", s=85)
axes[1, 1].set_yticks(y_n)
axes[1, 1].set_yticklabels(null_s["program_id"], fontsize=8.5)
axes[1, 1].set_xlabel("Null CAMERA t-statistic (DMSO Split A vs. Split B)")
axes[1, 1].set_title("D. Null Calibration Across 14 Aging Programs (0 BH Hits)")

plt.tight_layout()
save_fig(fig, "figS1_qc_and_null_calibration")

# ==============================================================================
# SUPPLEMENTARY FIGURE S2: Full 160-Compound (+ DFO) Atlas Heatmap
# ==============================================================================
fig_s2 = plt.figure(figsize=(15, 28))
ax_s2 = fig_s2.add_subplot(111)
annot_s2 = piv_q.map(lambda q: "*" if q < 0.05 else "")
sns.heatmap(piv_t, cmap="RdBu_r", center=0, vmin=-6, vmax=6, annot=annot_s2, fmt="", annot_kws={"size": 6.5},
            linewidths=0.2, linecolor="#eeeeee", ax=ax_s2, cbar_kws={"label": "CAMERA t-statistic"})
ax_s2.set_title("Supplementary Figure S2: Complete 160 Compound (+ DFO) x 14 Program Atlas (* Global BH q < 0.05)")
ax_s2.set_xticklabels(ax_s2.get_xticklabels(), rotation=40, ha="right", fontsize=9.5)
ax_s2.set_yticklabels(ax_s2.get_yticklabels(), fontsize=7)
save_fig(fig_s2, "figS2_full_160_compound_atlas")

# ==============================================================================
# SUPPLEMENTARY FIGURE S3: Independent Benchmarks, Sensitivity & Composition
# ==============================================================================
fig_s3, axes_s3 = plt.subplots(1, 3, figsize=(18, 5.6))

bench_df["any_sig"] = (bench_df["demagalhaes_up_q_BH"] < 0.05) | (bench_df["demagalhaes_dn_q_BH"] < 0.05)
sns.scatterplot(data=bench_df, x="demagalhaes_up_camera_t", y="demagalhaes_dn_camera_t", hue="any_sig",
                palette={False: "#999999", True: "#d62728"}, s=45, alpha=0.8, ax=axes_s3[0])
axes_s3[0].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes_s3[0].axvline(0, color="gray", linestyle="--", linewidth=0.8)
offsets_s3a = {"Binimetinib": (8, 8), "Baicalein": (6, 4), "DFO": (6, -8), "Simvastatin acid (ammonium)": (-110, 6)}
for d_lbl in ["Binimetinib", "Baicalein", "DFO", "Simvastatin acid (ammonium)"]:
    r_d = bench_df[bench_df["compound_id"] == d_lbl].iloc[0]
    axes_s3[0].annotate(d_lbl, (r_d["demagalhaes_up_camera_t"], r_d["demagalhaes_dn_camera_t"]), xytext=offsets_s3a[d_lbl], textcoords="offset points", fontsize=8)
axes_s3[0].set_title("A. RNA-Derived Aging Priors (de Magalhaes)")
axes_s3[0].set_xlabel("PRIOR_DEMAGALHAES_UP (CAMERA t)")
axes_s3[0].set_ylabel("PRIOR_DEMAGALHAES_DN (CAMERA t)")
axes_s3[0].legend(title="Global BH q < 0.05", frameon=True)

sns.scatterplot(data=sens_df, x="nes_combat_imp", y="camera_t", alpha=0.45, s=20, color="#1f77b4", ax=axes_s3[1])
axes_s3[1].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes_s3[1].axvline(0, color="gray", linestyle="--", linewidth=0.8)
axes_s3[1].set_xlabel("ComBat + Imputation GSEA NES")
axes_s3[1].set_ylabel("Plate-Blocked Observed CAMERA t")
axes_s3[1].set_title("B. Sensitivity: Plate-Blocked vs. ComBat + Imputation")

sns.scatterplot(data=comp_sens, x="cytosolic_ribo_mean_log2fc", y="oxphos_mean_log2fc", alpha=0.65, s=35, color="#2ca02c", ax=axes_s3[2])
axes_s3[2].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes_s3[2].axvline(0, color="gray", linestyle="--", linewidth=0.8)
for d_lbl in ["Binimetinib", "Baicalein", "DFO"]:
    r_d = comp_sens[comp_sens["compound_id"] == d_lbl].iloc[0]
    axes_s3[2].annotate(d_lbl, (r_d["cytosolic_ribo_mean_log2fc"], r_d["oxphos_mean_log2fc"]), xytext=(5, 4), textcoords="offset points", fontsize=8.5, fontweight="bold")
axes_s3[2].set_xlabel("Cytosolic Ribosome Mean log2FC (RPL/RPS)")
axes_s3[2].set_ylabel("OXPHOS Mean log2FC")
axes_s3[2].set_title("C. Composition Check: OXPHOS vs. Cytosolic Ribosomes")

plt.tight_layout()
save_fig(fig_s3, "figS3_benchmarks_and_sensitivity")

print("A08 updated: 4 main figures + 3 supplementary figures + Supplementary Table S1 generated.")
