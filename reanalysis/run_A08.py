import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

os.makedirs("reanalysis/A08/figures", exist_ok=True)
os.makedirs("reanalysis/A08/source_tables", exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 9.0,
    "ytick.labelsize": 9.0,
    "legend.fontsize": 8.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
})

def save_fig(fig, name):
    fig.savefig(f"reanalysis/A08/figures/{name}.png", bbox_inches="tight")
    fig.savefig(f"reanalysis/A08/figures/{name}.svg", bbox_inches="tight")
    plt.close(fig)

# Also update A04 null_calibration and candidate_stability with vif_camera_std
prog_df = pd.read_csv("reanalysis/A02/program_effects.tsv", sep="\t")
null_cal = pd.read_csv("reanalysis/A04/null_calibration.tsv", sep="\t")
prim = prog_df[prog_df["program_family"] == "primary_aging"].copy()
m_map = prim.groupby("program_id")["n_genes_measured"].first().to_dict()
vif_emp_map = prim.groupby("program_id")["vif_emp"].first().to_dict()
vif_std_map = prim.groupby("program_id")["vif_camera_std"].first().to_dict()

for idx, r in null_cal.iterrows():
    pid = r["program_id"]
    m_in = m_map[pid]
    m_out = 7871 - m_in
    ve = vif_emp_map[pid]
    vs = vif_std_map[pid]
    t_un = r["null_camera_t"] * np.sqrt((ve / m_in + 1.0 / m_out) / (1.0 / m_in + 1.0 / m_out))
    t_std = t_un * np.sqrt((1.0 / m_in + 1.0 / m_out) / (vs / m_in + 1.0 / m_out))
    p_std = 2.0 * stats.t.sf(abs(t_std), df=7869)
    null_cal.loc[idx, "null_camera_t_emp"] = r["null_camera_t"]
    null_cal.loc[idx, "null_camera_t"] = round(float(t_std), 4)
    null_cal.loc[idx, "null_p_camera"] = float(p_std)
_, null_cal["null_q_BH"], _, _ = multipletests(null_cal["null_p_camera"], method="fdr_bh")
null_cal.to_csv("reanalysis/A04/null_calibration.tsv", sep="\t", index=False)

# Load core A01-A07 tables
sm = pd.read_csv("reanalysis/A01/sample_manifest.tsv", sep="\t")
sm_12 = sm[sm["plate"] != "RedMix"].copy()
sa = pd.read_csv("reanalysis/A01/sample_accounting.tsv", sep="\t")
sa_12 = sa[sa["plate"] != "RedMix"].copy()
qc = pd.read_csv("reanalysis/A04/qc_metrics.tsv", sep="\t")
qc_12 = qc[qc["plate"] != "RedMix"].copy()
rep = pd.read_csv("reanalysis/A04/replication_metrics.tsv", sep="\t")
pman = pd.read_csv("reanalysis/A03/program_manifest.tsv", sep="\t")
pov = pd.read_csv("reanalysis/A03/program_overlap.tsv", sep="\t")
spec_df = pd.read_csv("reanalysis/A06/candidate_specificity.tsv", sep="\t")
bench_df = pd.read_csv("reanalysis/A05/drug_reference_comparisons.tsv", sep="\t")
sens_df = pd.read_csv("reanalysis/A04/sensitivity_effects.tsv", sep="\t")
comp_sens = pd.read_csv("reanalysis/A06/composition_sensitivity.tsv", sep="\t")

aging_ids = [
    "PROG_DNA_REPAIR", "PROG_TELOMERE_MAINT", "PROG_CHROMATIN_ORG",
    "PROG_PROTEOSTASIS_UPR", "PROG_PROTEASOME", "PROG_MACROAUTOPHAGY",
    "PROG_OXPHOS", "PRIOR_DEMAGALHAES_DN",
    "PROG_MTORC1_SIGNALING", "PROG_SENESCENCE_CORE",
    "PROG_CELL_CELL_COMM", "PROG_INFLAMMATORY", "PROG_INNATE_TLR",
    "PRIOR_DEMAGALHAES_UP",
]

signs_14 = {
    "PROG_DNA_REPAIR": +1, "PROG_TELOMERE_MAINT": +1, "PROG_CHROMATIN_ORG": +1,
    "PROG_PROTEOSTASIS_UPR": +1, "PROG_PROTEASOME": +1, "PROG_MACROAUTOPHAGY": +1,
    "PROG_OXPHOS": +1, "PRIOR_DEMAGALHAES_DN": +1,
    "PROG_MTORC1_SIGNALING": -1, "PROG_SENESCENCE_CORE": -1,
    "PROG_CELL_CELL_COMM": -1, "PROG_INFLAMMATORY": -1, "PROG_INNATE_TLR": -1,
    "PRIOR_DEMAGALHAES_UP": -1,
}
dir_labels = {k: f"{k} ({'Up [+1]' if v > 0 else 'Down [-1]'})" for k, v in signs_14.items()}

# ==============================================================================
# MAIN FIGURE 1: Number of proteins per hallmark/program vs all genes in set
# ==============================================================================
pman_14 = pman[pman["program_id"].isin(aging_ids)].set_index("program_id").loc[aging_ids].reset_index()
fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.4), gridspec_kw={"width_ratios": [1.18, 1.0]})

y_pos = np.arange(len(pman_14))
axes[0].hlines(y=y_pos, xmin=pman_14["measured_retained_7871"], xmax=pman_14["total_genes"], color="#bdbdbd", linewidth=2.5)
axes[0].scatter(pman_14["total_genes"], y_pos, color="#7f7f7f", s=90, label="Total Annotated Genes in Set", zorder=3)
colors_dir = ["#1f77b4" if signs_14[pid] > 0 else "#d62728" for pid in pman_14["program_id"]]
axes[0].scatter(pman_14["measured_retained_7871"], y_pos, color=colors_dir, s=95, label="Measured Proteins (Blue: Up [+1], Red: Down [-1])", zorder=4)
for i, row in pman_14.iterrows():
    axes[0].text(row["total_genes"] + 5, i, f"{int(row['measured_retained_7871'])}/{int(row['total_genes'])} ({row['coverage_pct_retained']:.0f}%)", va="center", fontsize=8.5)
axes[0].set_yticks(y_pos)
axes[0].set_yticklabels([dir_labels[pid] for pid in pman_14["program_id"]], fontsize=8.8)
axes[0].invert_yaxis()
axes[0].set_xlim(0, 400)
axes[0].set_xlabel("Gene / Protein Count")
axes[0].set_title("A. Measured Proteins vs. Annotated Genes & Target Direction")
axes[0].legend(loc="lower right", frameon=True, fontsize=8.2)

pov_14 = pd.DataFrame(0.0, index=aging_ids, columns=aging_ids)
np.fill_diagonal(pov_14.values, 1.0)
for _, r in pov.iterrows():
    if r["program_1"] in aging_ids and r["program_2"] in aging_ids:
        pov_14.loc[r["program_1"], r["program_2"]] = r["jaccard"]
        pov_14.loc[r["program_2"], r["program_1"]] = r["jaccard"]
sns.heatmap(pov_14, cmap="Blues", vmin=0, vmax=0.35, annot=True, fmt=".2f", annot_kws={"size": 7.2},
            ax=axes[1], cbar_kws={"label": "Jaccard Overlap"})
axes[1].set_title("B. Pairwise Protein Membership Overlap (Jaccard)")
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha="right", fontsize=8.2)
axes[1].set_yticklabels(axes[1].get_yticklabels(), fontsize=8.2)

plt.tight_layout()
save_fig(fig, "fig1_program_coverage")
pman_14.to_csv("reanalysis/A08/source_tables/fig1_program_coverage.tsv", sep="\t", index=False)

# ==============================================================================
# MAIN FIGURE 2: Highlighted compounds across 14 programs + Composite Score
# ==============================================================================
piv_t = prim.pivot(index="compound_id", columns="program_id", values="camera_t")[aging_ids]
piv_nes = prim.pivot(index="compound_id", columns="program_id", values="nes_gsea")[aging_ids]
piv_q = prim.pivot(index="compound_id", columns="program_id", values="q_BH_global_camera")[aging_ids]

highlight_drugs = [
    # Top multi-pathway geroprotective candidates (by CGS_combined)
    "Binimetinib",
    "Cyclothiazide",
    "Zanubrutinib",
    "Baicalein",
    "Almonertinib (hydrochloride)",
    "Ferulic acid",
    "α-Lipoic Acid",
    "Digitoxin",
    "Osimertinib",
    "Orelabrutinib",
    # Canonical longevity / metabolic / epigenetic comparators
    "Atorvastatin (hemicalcium salt)",
    "Eicosapentaenoic Acid",
    "Simvastatin acid (ammonium)",
    "Methylene blue (trihydrate)",
    "Lonafarnib",
    "Nordihydroguaiaretic acid",
    "Dimethyl fumarate",
    # Negative / pro-aging / stress comparators
    "Erlotinib (Hydrochloride)",
    "Asciminib",
    "Lapatinib (ditosylate)",
    "Icotinib",
    "Desloratadine",
    "Amiodarone (hydrochloride)",
    "DFO",
]
highlight_drugs = [d for d in highlight_drugs if d in piv_t.index]
# Sort by CGS_combined descending
spec_idx = spec_df.set_index("compound_id")
highlight_drugs = sorted(highlight_drugs, key=lambda d: spec_idx.loc[d, "CGS_combined"], reverse=True)

sub_t = piv_t.loc[highlight_drugs]
sub_q = piv_q.loc[highlight_drugs]
annot_sub = sub_q.map(lambda q: "**" if q < 0.01 else ("*" if q < 0.05 else ("+" if q < 0.20 else "")))

from matplotlib.lines import Line2D

fig, axes = plt.subplots(1, 2, figsize=(17.0, 9.6), gridspec_kw={"width_ratios": [2.25, 1.0]}, sharey=False)

sns.heatmap(sub_t, cmap="RdBu_r", center=0, vmin=-6.5, vmax=6.5, annot=annot_sub, fmt="", annot_kws={"size": 9.0, "weight": "bold"},
            linewidths=0.4, linecolor="white", ax=axes[0], cbar_kws={"label": "Correlation-Aware CAMERA t-statistic", "pad": 0.02})
axes[0].axvline(8.0, color="black", linestyle="-", linewidth=1.8)
axes[0].set_title("A. Compound Responses Across 14 Programs (** q<0.01, * q<0.05, + q<0.20)")
axes[0].set_xticklabels([dir_labels[c] for c in aging_ids], rotation=40, ha="right", fontsize=8.5)
axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=9.0, rotation=0)
axes[0].set_xlabel("")
axes[0].set_ylabel("Compound (Ordered by Composite Geroprotective Score)")

# Panel B: Aligned Cleveland dot/lollipop plot of Composite Geroprotective Scores (CGS)
y_cgs = np.arange(len(highlight_drugs)) + 0.5
sub_spec = spec_idx.loc[highlight_drugs]

axes[1].axvline(0, color="black", linestyle="-", linewidth=1.0)
# Symmetric tiered BH FDR threshold lines on both positive and negative sides
axes[1].axvline(2.38, color="#555555", linestyle="--", linewidth=0.95, alpha=0.75)
axes[1].axvline(-2.38, color="#555555", linestyle="--", linewidth=0.95, alpha=0.75)
axes[1].axvline(2.85, color="#222222", linestyle=":", linewidth=1.15, alpha=0.85)
axes[1].axvline(-2.85, color="#222222", linestyle=":", linewidth=1.15, alpha=0.85)

axes[1].hlines(y=y_cgs, xmin=np.minimum(0, np.minimum(sub_spec["CGS_unadj"], sub_spec["CGS_combined"])),
               xmax=np.maximum(0, np.maximum(sub_spec["CGS_adj"], sub_spec["CGS_combined"])),
               color="#cccccc", linewidth=2.0, zorder=2)
axes[1].scatter(sub_spec["CGS_unadj"], y_cgs, color="#1f77b4", marker="o", s=50, alpha=0.88, zorder=3)
axes[1].scatter(sub_spec["CGS_adj"], y_cgs, color="#ff7f0e", marker="^", s=54, alpha=0.88, zorder=4)
c_comb = ["#2ca02c" if v >= 0 else "#d62728" for v in sub_spec["CGS_combined"]]
axes[1].scatter(sub_spec["CGS_combined"], y_cgs, color=c_comb, marker="D", s=66, edgecolors="k", linewidths=0.7, zorder=5)

for idx_y, (d_name, r_s) in enumerate(sub_spec.iterrows()):
    val = r_s["CGS_combined"]
    q_c = r_s["q_BH_CGS_combined"]
    sig_str = "**" if q_c < 0.01 else ("*" if q_c < 0.05 else ("+" if q_c < 0.10 else ("~" if q_c < 0.20 else "")))
    n_b = int(max(r_s["n_beneficial_unadj"], r_s["n_beneficial_adj"]))
    if val >= 0:
        rightmost = max(0.0, float(r_s["CGS_unadj"]), float(r_s["CGS_adj"]), float(val))
        axes[1].text(rightmost + 0.22, idx_y + 0.5, f"{val:+.2f}{sig_str} ({n_b}/14)", va="center", ha="left", fontsize=7.8, fontweight="bold")
    else:
        leftmost = min(0.0, float(r_s["CGS_unadj"]), float(r_s["CGS_adj"]), float(val))
        axes[1].text(leftmost - 0.22, idx_y + 0.5, f"{val:+.2f}{sig_str} ({n_b}/14)", va="center", ha="right", fontsize=7.8, fontweight="bold")

axes[1].set_ylim(len(highlight_drugs), 0)
axes[1].set_yticks(y_cgs)
axes[1].set_yticklabels([])
axes[1].set_xlim(-9.3, 6.6)
axes[1].set_xlabel("14-Pathway Composite Geroprotective Score (Stouffer Z)")
axes[1].set_title("B. Composite Geroprotective Score (CGS)")
axes[1].grid(axis="x", linestyle="--", alpha=0.35)

legend_handles_2b = [
    Line2D([0], [0], color="#555555", linestyle="--", linewidth=0.95, label="BH q < 0.10 (+/-2.38)"),
    Line2D([0], [0], color="#222222", linestyle=":", linewidth=1.15, label="BH q < 0.05 (+/-2.85)"),
    Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="#1f77b4", markeredgecolor="#1f77b4", markersize=7.0, label="CGS Unadjusted"),
    Line2D([0], [0], marker="^", linestyle="None", markerfacecolor="#ff7f0e", markeredgecolor="#ff7f0e", markersize=7.5, label="CGS Prolif-Adjusted"),
    Line2D([0], [0], marker="D", linestyle="None", markerfacecolor="none", markeredgecolor="black", markeredgewidth=1.1, markersize=7.5, label="CGS Combined (Consensus)"),
]
axes[1].legend(handles=legend_handles_2b, loc="lower right", frameon=True, fontsize=7.8)

plt.tight_layout()
save_fig(fig, "fig2_highlighted_heatmap")

# Main Supplemental Table S1 (all 160 library compounds + DFO)
supp_t1 = prim.merge(
    spec_df[["compound_id", "CGS_unadj", "q_BH_CGS_unadj", "CGS_adj", "q_BH_CGS_adj", "CGS_combined", "q_BH_CGS_combined", "n_beneficial_unadj", "n_beneficial_adj"]],
    on="compound_id", how="left"
).sort_values(["CGS_combined", "compound_id", "program_id"], ascending=[False, True, True]).copy()
supp_t1.to_csv("reanalysis/A08/source_tables/Supplementary_Table_S1_all_compounds.tsv", sep="\t", index=False)

# ==============================================================================
# MAIN FIGURE 3: NES Distributions, Composite Score Landscape & Specificity
# ==============================================================================
prim["sig_gsea_q05"] = prim["q_BH_global_gsea"] < 0.05
fig, axes = plt.subplots(1, 3, figsize=(19.0, 6.0), gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]})

sns.boxplot(data=prim, x="program_id", y="nes_gsea", order=aging_ids, color="#f0f0f0", fliersize=0, ax=axes[0])
sns.stripplot(data=prim[~prim["sig_gsea_q05"]], x="program_id", y="nes_gsea", order=aging_ids, color="#999999", alpha=0.35, size=2.5, ax=axes[0])
sns.stripplot(data=prim[prim["sig_gsea_q05"]], x="program_id", y="nes_gsea", order=aging_ids, color="#d62728", alpha=0.75, size=3.5, ax=axes[0])
axes[0].axhline(0, color="black", linestyle="--", linewidth=0.9)
axes[0].axvline(7.5, color="black", linestyle=":", linewidth=1.2)
axes[0].set_xticks(range(len(aging_ids)))
axes[0].set_xticklabels([dir_labels[c] for c in aging_ids], rotation=45, ha="right", fontsize=7.8)
axes[0].set_title("A. GSEA NES Distributions Across 161 Contrasts")
axes[0].set_ylabel("Normalized Enrichment Score (NES)")
axes[0].set_xlabel("")

# Panel B: CGS_unadj vs CGS_adj across all 161 compounds
sc_b = axes[1].scatter(spec_df["CGS_unadj"], spec_df["CGS_adj"], c=spec_df["CGS_combined"],
                       cmap="RdYlGn", vmin=-5.0, vmax=4.0, s=46, alpha=0.88, edgecolors="k", linewidths=0.4)
plt.colorbar(sc_b, ax=axes[1], label="Combined CGS (Stouffer Z)")
axes[1].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes[1].axvline(0, color="gray", linestyle="--", linewidth=0.8)
x_grid = np.linspace(-6.5, 4.0, 100)
axes[1].plot(x_grid, 2.38 * np.sqrt(2) - x_grid, color="#006d2c", linestyle="--", linewidth=0.9, alpha=0.7, label="Combined BH q < 0.10")
axes[1].plot(x_grid, -2.85 * np.sqrt(2) - x_grid, color="#a50f15", linestyle="--", linewidth=0.9, alpha=0.7, label="Pro-Aging BH q < 0.05")
axes[1].set_xlim(-6.5, 4.5)
axes[1].set_ylim(-3.9, 4.3)

offsets_3b = {
    "Binimetinib": (-60, 12),
    "Cyclothiazide": (-72, 15),
    "Zanubrutinib": (-75, -18),
    "Baicalein": (-62, 10),
    "Almonertinib (hydrochloride)": (-82, 10),
    "Ferulic acid": (-68, -28),
    "DFO": (10, -4),
    "Amiodarone (hydrochloride)": (12, -14),
    "Erlotinib (Hydrochloride)": (-58, 14),
}
for d_lbl, xy_off in offsets_3b.items():
    if d_lbl in spec_idx.index:
        r_d = spec_idx.loc[d_lbl]
        short_lbl = d_lbl.split(" (")[0]
        axes[1].annotate(short_lbl, (r_d["CGS_unadj"], r_d["CGS_adj"]), xytext=xy_off, textcoords="offset points",
                         fontsize=7.6, fontweight="bold", arrowprops=dict(arrowstyle="-", color="#555555", lw=0.6))
axes[1].set_xlabel("Unadjusted 14-Pathway CGS (vs. DMSO)")
axes[1].set_ylabel("Proliferation-Adjusted 14-Pathway CGS")
axes[1].set_title("B. Multi-Pathway Composite Geroprotective Score")
axes[1].legend(loc="upper left", frameon=True, fontsize=7.6)

# Panel C: Specificity vs Proliferation & OXPHOS Axes
sc_c = axes[2].scatter(spec_df["NUIS_E2F_PROLIF"], spec_df["PROG_INFLAMMATORY"], c=spec_df["PROG_OXPHOS"],
                       cmap="RdBu_r", vmin=-6.5, vmax=6.5, s=46, alpha=0.85, edgecolors="k", linewidths=0.4)
plt.colorbar(sc_c, ax=axes[2], label="PROG_OXPHOS (CAMERA t)")
offsets_3c = {
    "Binimetinib": (-18, 10),
    "Baicalein": (8, -6),
    "Zanubrutinib": (-68, -14),
    "DFO": (8, 5),
    "Erlotinib (Hydrochloride)": (8, 5),
}
for d_lbl, xy_off in offsets_3c.items():
    if d_lbl in spec_idx.index:
        r_d = spec_idx.loc[d_lbl]
        short_lbl = d_lbl.split(" (")[0]
        axes[2].annotate(short_lbl, (r_d["NUIS_E2F_PROLIF"], r_d["PROG_INFLAMMATORY"]), xytext=xy_off, textcoords="offset points", fontsize=8.0, fontweight="bold")
axes[2].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes[2].axvline(0, color="gray", linestyle="--", linewidth=0.8)
axes[2].set_xlim(-13.2, 3.2)
axes[2].set_xlabel("Proliferation Proxy (NUIS_E2F_PROLIF CAMERA t)")
axes[2].set_ylabel("Intracellular Inflammation (PROG_INFLAMMATORY CAMERA t)")
axes[2].set_title("C. Specificity vs. Proliferation & OXPHOS Axes")

plt.tight_layout()
save_fig(fig, "fig3_nes_distributions")
piv_nes.to_csv("reanalysis/A08/source_tables/fig3_nes_all_compounds.tsv", sep="\t")
piv_t.to_csv("reanalysis/A08/source_tables/fig3_camera_t_all_161.tsv", sep="\t")

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
    "Zanubrutinib",
    "Almonertinib (hydrochloride)",
    "Erlotinib (Hydrochloride)",
    "Mobocertinib (succinate)",
    "Lapatinib (ditosylate)",
]
ox_sentinels = [g for g in ["NDUFS1", "IDH2", "ACO2", "SDHA", "COX4I1", "ATP5F1A"] if g in ret_prots]
inf_sentinels = [g for g in ["ICAM1", "LDLR", "NFKB1", "RELA", "STAT3", "ATP2B1"] if g in ret_prots]

pal_k = {
    "DMSO": "#7f7f7f",
    "Binimetinib": "#2ca02c",
    "Zanubrutinib": "#17becf",
    "Almonertinib (hydrochloride)": "#9467bd",
    "Erlotinib (Hydrochloride)": "#d62728",
    "Mobocertinib (succinate)": "#1f77b4",
    "Lapatinib (ditosylate)": "#ff7f0e",
}

m_k = sm["compound_id"].isin(kinase_drugs).values & (~sm["excluded"].values) & (sm["plate"].values != "RedMix")
fig, axes = plt.subplots(2, 1, figsize=(15.0, 10.6))

for ax, plist, title in [
    (axes[0], ox_sentinels, "A. Kinase Inhibitor Panel: Respiratory / OXPHOS Sentinels (Zero-Centered, Unit-Variance Z-Score)"),
    (axes[1], inf_sentinels, "B. Kinase Inhibitor Panel: Inflammatory Sentinels (Zero-Centered, Unit-Variance Z-Score)"),
]:
    p_idx = [ret_prots.index(g) for g in plist]
    sub_mat = X_norm[np.ix_(m_k, p_idx)]
    mu = np.nanmean(sub_mat, axis=0, keepdims=True)
    sd = np.nanstd(sub_mat, axis=0, ddof=1, keepdims=True)
    sub_z = (sub_mat - mu) / sd
    df_w = pd.DataFrame(sub_z, columns=plist)
    df_w["Drug"] = sm.loc[m_k, "compound_id"].values
    df_l = df_w.melt(id_vars="Drug", var_name="Protein", value_name="Scaled Abundance (Z-score)")
    df_l["Drug"] = pd.Categorical(df_l["Drug"], categories=kinase_drugs, ordered=True)
    sns.violinplot(data=df_l, x="Protein", y="Scaled Abundance (Z-score)", hue="Drug", palette=pal_k, inner="box", cut=0, alpha=0.30, ax=ax)
    sns.stripplot(data=df_l, x="Protein", y="Scaled Abundance (Z-score)", hue="Drug", palette=pal_k, dodge=True, size=4.2, alpha=0.9, jitter=0.1, ax=ax)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[:len(kinase_drugs)], l[:len(kinase_drugs)], title="Treatment (n=5-6 plates; DMSO n=23)", loc="lower right", frameon=True, ncol=4, fontsize=8.2)
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
axes[1, 1].set_xlim(-1.0, 1.0)
axes[1, 1].set_xlabel("Null CAMERA t-statistic (DMSO Split A vs. Split B)")
axes[1, 1].set_title("D. Null Calibration Across 14 Aging Programs (0 BH Hits)")

plt.tight_layout()
save_fig(fig, "figS1_qc_and_null_calibration")

# ==============================================================================
# SUPPLEMENTARY FIGURE S2: Full 160-Compound (+ DFO) Atlas Heatmap (Ordered by CGS)
# ==============================================================================
order_all = spec_df.sort_values("CGS_combined", ascending=False)["compound_id"].tolist()
piv_t_ord = piv_t.loc[order_all]
piv_q_ord = piv_q.loc[order_all]
fig_s2 = plt.figure(figsize=(15, 28))
ax_s2 = fig_s2.add_subplot(111)
annot_s2 = piv_q_ord.map(lambda q: "*" if q < 0.05 else "")
sns.heatmap(piv_t_ord, cmap="RdBu_r", center=0, vmin=-6.5, vmax=6.5, annot=annot_s2, fmt="", annot_kws={"size": 6.5},
            linewidths=0.2, linecolor="#eeeeee", ax=ax_s2, cbar_kws={"label": "CAMERA t-statistic"})
ax_s2.set_title("Supplementary Figure S2: Complete 160 Compound (+ DFO) x 14 Program Atlas Ordered by CGS (* Global BH q < 0.05)")
ax_s2.set_xticklabels([dir_labels[c] for c in aging_ids], rotation=40, ha="right", fontsize=9.0)
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
