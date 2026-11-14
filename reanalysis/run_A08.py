import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs("reanalysis/A08/figures", exist_ok=True)
os.makedirs("reanalysis/A08/source_tables", exist_ok=True)

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.labelsize": 11.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9.5,
})

# Load frozen tables
acct = pd.read_csv("reanalysis/A01/sample_accounting.tsv", sep="\t")
elig = pd.read_csv("reanalysis/A01/contrast_eligibility.tsv", sep="\t")
prog_man = pd.read_csv("reanalysis/A03/program_manifest.tsv", sep="\t")
prog_ov = pd.read_csv("reanalysis/A03/program_overlap.tsv", sep="\t")
prog_df = pd.read_csv("reanalysis/A02/program_effects.tsv", sep="\t")
prot_df = pd.read_csv("reanalysis/A02/protein_effects.tsv", sep="\t")
bench_df = pd.read_csv("reanalysis/A05/drug_reference_comparisons.tsv", sep="\t")
spec_df = pd.read_csv("reanalysis/A06/candidate_specificity.tsv", sep="\t")
comp_sens = pd.read_csv("reanalysis/A06/composition_sensitivity.tsv", sep="\t")
ev_df = pd.read_csv("reanalysis/A07/candidate_evidence.tsv", sep="\t")
qc_df = pd.read_csv("reanalysis/A04/qc_metrics.tsv", sep="\t")
rep_df = pd.read_csv("reanalysis/A04/replication_metrics.tsv", sep="\t")
sens_df = pd.read_csv("reanalysis/A04/sensitivity_effects.tsv", sep="\t")
null_df = pd.read_csv("reanalysis/A04/null_calibration.tsv", sep="\t")

def save_fig(fig, name):
    fig.savefig(f"reanalysis/A08/figures/{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"reanalysis/A08/figures/{name}.svg", bbox_inches="tight")
    plt.close(fig)

# ==============================================================================
# FIGURE 1: Dataset provenance and performance
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(15, 11))

# Panel A: Sample accounting per plate (dot/lollipop per role, no barplot)
acct_m = acct.melt(id_vars=["plate", "total"], var_name="role", value_name="count")
acct_m = acct_m[acct_m["count"] > 0]
sns.scatterplot(data=acct_m, x="plate", y="count", hue="role", style="role", s=90, ax=axes[0, 0], palette="Set2")
axes[0, 0].set_title("A. Sample Composition Across 13 Plates (Highlighting RedMix)")
axes[0, 0].set_ylabel("Sample Count")
axes[0, 0].set_xlabel("Plate")
axes[0, 0].tick_params(axis="x", rotation=45)
axes[0, 0].legend(title="Role", loc="center right", frameon=True)

# Panel B: QC per-protein SD by plate (dot/lollipop)
sns.scatterplot(data=qc_df, x="plate", y="qc_median_protein_sd", s=100, color="#1f77b4", ax=axes[0, 1])
axes[0, 1].axhline(qc_df["qc_median_protein_sd"].median(), color="gray", linestyle="--", linewidth=1)
axes[0, 1].set_title("B. Median Per-Protein SD Across QC Samples by Plate")
axes[0, 1].set_ylabel("Median QC Protein SD (log2)")
axes[0, 1].set_xlabel("Plate")
axes[0, 1].tick_params(axis="x", rotation=45)

# Panel C: Leave-one-plate-out replication Pearson r
rep_plot = rep_df.sort_values("mean_lopo_pearson_r", ascending=True)
y_p = np.arange(len(rep_plot))
axes[1, 0].hlines(y=y_p, xmin=rep_plot["min_lopo_pearson_r"], xmax=rep_plot["max_lopo_pearson_r"], color="#bbbbbb", linewidth=2.5)
axes[1, 0].scatter(rep_plot["mean_lopo_pearson_r"], y_p, color="#2ca02c", s=90, zorder=3, label="Mean LOPO r")
axes[1, 0].set_yticks(y_p)
axes[1, 0].set_yticklabels(rep_plot["compound_id"])
axes[1, 0].set_xlim(0, 1.0)
axes[1, 0].set_xlabel("Leave-One-Plate-Out Pearson r (log2FC)")
axes[1, 0].set_title("C. Cross-Plate Replication of Proteome Perturbations")
axes[1, 0].grid(axis="x", linestyle="--", alpha=0.4)

# Panel D: Null calibration (DMSO_A vs DMSO_B across 11 plates)
null_s = null_df.sort_values("null_camera_t")
y_n = np.arange(len(null_s))
axes[1, 1].axvline(0, color="black", linestyle="--", linewidth=1)
axes[1, 1].scatter(null_s["null_camera_t"], y_n, color="#7f7f7f", s=75)
axes[1, 1].set_yticks(y_n)
axes[1, 1].set_yticklabels(null_s["program_id"], fontsize=8.5)
axes[1, 1].set_xlabel("Null CAMERA t-statistic (DMSO Split A vs. Split B)")
axes[1, 1].set_title("D. Null Calibration Across 14 Aging Programs (0 BH Hits)")

plt.tight_layout()
save_fig(fig, "fig1_provenance_and_qc")
acct.to_csv("reanalysis/A08/source_tables/fig1_sample_accounting.tsv", sep="\t", index=False)
null_df.to_csv("reanalysis/A08/source_tables/fig1_null_calibration.tsv", sep="\t", index=False)

# ==============================================================================
# FIGURE 2: Evidence and coverage
# ==============================================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 6.8), gridspec_kw={"width_ratios": [1.15, 1.0]})

pm_aging = prog_man[~prog_man["program_id"].str.startswith("NUIS_")].iloc[::-1].reset_index(drop=True)
y_g = np.arange(len(pm_aging))
axes[0].hlines(y=y_g, xmin=pm_aging["measured_retained_7871"], xmax=pm_aging["total_genes"], color="#bbbbbb", linewidth=2)
axes[0].scatter(pm_aging["total_genes"], y_g, color="#7f7f7f", s=75, label="Total Genes")
axes[0].scatter(pm_aging["measured_retained_7871"], y_g, color="#1f77b4", s=95, label="Measured (7,871 retained)")
for idx, r in pm_aging.iterrows():
    axes[0].text(r["total_genes"] + 4, idx, f"{r['measured_retained_7871']}/{r['total_genes']} ({r['coverage_pct_retained']:.0f}%) [{r['aging_evidence_grade']}]", va="center", fontsize=8)
axes[0].set_yticks(y_g)
axes[0].set_yticklabels(pm_aging["program_id"], fontsize=9)
axes[0].set_xlim(0, pm_aging["total_genes"].max() + 115)
axes[0].set_xlabel("Gene / Protein Count")
axes[0].set_title("A. Program Coverage and Evidence Grade")
axes[0].legend(loc="lower right", frameon=True)

# Panel B: Jaccard overlap matrix
aging_ids = prog_man.loc[~prog_man["program_id"].str.startswith("NUIS_"), "program_id"].tolist()
j_mat = pd.DataFrame(0.0, index=aging_ids, columns=aging_ids)
np.fill_diagonal(j_mat.values, 1.0)
for _, r in prog_ov.iterrows():
    if r["program_1"] in aging_ids and r["program_2"] in aging_ids:
        j_mat.loc[r["program_1"], r["program_2"]] = r["jaccard"]
        j_mat.loc[r["program_2"], r["program_1"]] = r["jaccard"]
sns.heatmap(j_mat, cmap="Blues", vmin=0, vmax=0.35, annot=True, fmt=".2f", annot_kws={"size": 7.5}, ax=axes[1], cbar_kws={"label": "Jaccard Overlap"})
axes[1].set_title("B. Pairwise Protein Membership Overlap (Jaccard)")
axes[1].tick_params(axis="x", rotation=45)

plt.tight_layout()
save_fig(fig, "fig2_evidence_and_coverage")
pm_aging.to_csv("reanalysis/A08/source_tables/fig2_program_coverage.tsv", sep="\t", index=False)

# ==============================================================================
# FIGURE 3: Compound response atlas (Main subset + distribution + specificity)
# ==============================================================================
prim = prog_df[prog_df["program_family"] == "primary_aging"].copy()
piv_t = prim.pivot(index="compound_id", columns="program_id", values="camera_t")[aging_ids]
piv_q = prim.pivot(index="compound_id", columns="program_id", values="q_BH_global_camera")[aging_ids]

# Highlighted panel of 20 directly estimable drugs (literature + top hits)
hi_drugs = [d for d in ["DFO", "Binimetinib", "Baicalein", "Nordihydroguaiaretic acid", "Lonafarnib", "Nicotinamide",
                        "Simvastatin acid (ammonium)", "Atorvastatin (hemicalcium salt)", "Eicosapentaenoic Acid",
                        "Dimethyl fumarate", "Tucidinostat", "Methylene blue (trihydrate)", "Mobocertinib (succinate)",
                        "Erlotinib (Hydrochloride)", "Asciminib", "Flumatinib", "Lapatinib (ditosylate)",
                        "Triclabendazole", "Mifepristone", "6-Mercaptopurine"] if d in piv_t.index]

fig = plt.figure(figsize=(17, 11))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], width_ratios=[1.1, 1.0])

ax_hm = fig.add_subplot(gs[0, :])
annot_hm = piv_q.loc[hi_drugs].map(lambda q: "**" if q < 0.01 else ("*" if q < 0.05 else ""))
sns.heatmap(piv_t.loc[hi_drugs], cmap="RdBu_r", center=0, vmin=-6, vmax=6, annot=annot_hm, fmt="",
            annot_kws={"size": 9, "weight": "bold"}, linewidths=0.3, linecolor="#e0e0e0", ax=ax_hm,
            cbar_kws={"label": "Correlation-Aware CAMERA t-statistic"})
ax_hm.set_title("A. Highlighted Directly Estimable Compounds Across 14 Aging Programs (* Global BH q < 0.05, ** q < 0.01)")
ax_hm.set_xticklabels(ax_hm.get_xticklabels(), rotation=35, ha="right", fontsize=9.5)
ax_hm.set_xlabel("")
ax_hm.set_ylabel("Compound")

# Panel B: Distributions across all 161 compounds
ax_dist = fig.add_subplot(gs[1, 0])
prim["sig_q05"] = prim["q_BH_global_camera"] < 0.05
sns.boxplot(data=prim, x="program_id", y="camera_t", order=aging_ids, color="#f0f0f0", fliersize=0, ax=ax_dist)
sns.stripplot(data=prim[~prim["sig_q05"]], x="program_id", y="camera_t", order=aging_ids, color="#999999", alpha=0.35, size=2.5, ax=ax_dist)
sns.stripplot(data=prim[prim["sig_q05"]], x="program_id", y="camera_t", order=aging_ids, color="#d62728", alpha=0.8, size=4, ax=ax_dist)
ax_dist.axhline(0, color="black", linestyle="--", linewidth=0.9)
ax_dist.set_xticklabels(ax_dist.get_xticklabels(), rotation=45, ha="right", fontsize=8.5)
ax_dist.set_title("B. CAMERA t-statistic Distribution Across All 161 Compounds")
ax_dist.set_ylabel("CAMERA t-statistic")
ax_dist.set_xlabel("")

# Panel C: Specificity scatter (NUIS_E2F_PROLIF vs PROG_INFLAMMATORY, colored by PROG_OXPHOS)
ax_sp = fig.add_subplot(gs[1, 1])
sc_p = ax_sp.scatter(spec_df["NUIS_E2F_PROLIF"], spec_df["PROG_INFLAMMATORY"], c=spec_df["PROG_OXPHOS"],
                     cmap="RdBu_r", vmin=-6, vmax=6, s=45, alpha=0.85, edgecolors="k", linewidths=0.4)
plt.colorbar(sc_p, ax=ax_sp, label="PROG_OXPHOS (CAMERA t)")
offsets_3c = {"Binimetinib": (-12, 8), "Baicalein": (8, -10), "DFO": (6, 4), "Erlotinib (Hydrochloride)": (6, 4)}
for d_lbl in ["Binimetinib", "Baicalein", "DFO", "Erlotinib (Hydrochloride)"]:
    r_d = spec_df[spec_df["compound_id"] == d_lbl].iloc[0]
    ax_sp.annotate(d_lbl, (r_d["NUIS_E2F_PROLIF"], r_d["PROG_INFLAMMATORY"]), xytext=offsets_3c[d_lbl], textcoords="offset points", fontsize=8.5, fontweight="bold")
ax_sp.axhline(0, color="gray", linestyle="--", linewidth=0.8)
ax_sp.axvline(0, color="gray", linestyle="--", linewidth=0.8)
ax_sp.set_xlabel("Proliferation Proxy (NUIS_E2F_PROLIF CAMERA t)")
ax_sp.set_ylabel("Intracellular Inflammation (PROG_INFLAMMATORY CAMERA t)")
ax_sp.set_title("C. Specificity vs. Proliferation & Shared Stress Axes")

plt.tight_layout()
save_fig(fig, "fig3_compound_response_atlas")
piv_t.to_csv("reanalysis/A08/source_tables/fig3_camera_t_all_161.tsv", sep="\t")

# ==============================================================================
# FIGURE 4: Independent aging benchmarks
# ==============================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

bench_df["any_sig"] = (bench_df["demagalhaes_up_q_BH"] < 0.05) | (bench_df["demagalhaes_dn_q_BH"] < 0.05)
sns.scatterplot(data=bench_df, x="demagalhaes_up_camera_t", y="demagalhaes_dn_camera_t", hue="any_sig",
                palette={False: "#999999", True: "#d62728"}, s=50, alpha=0.8, ax=axes[0])
axes[0].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes[0].axvline(0, color="gray", linestyle="--", linewidth=0.8)
offsets_4a = {"Binimetinib": (8, 8), "Baicalein": (6, 4), "Nordihydroguaiaretic acid": (6, -8), "DFO": (6, -8), "Simvastatin acid (ammonium)": (-115, 6)}
for d_lbl in ["Binimetinib", "Baicalein", "Nordihydroguaiaretic acid", "DFO", "Simvastatin acid (ammonium)"]:
    r_d = bench_df[bench_df["compound_id"] == d_lbl].iloc[0]
    axes[0].annotate(d_lbl, (r_d["demagalhaes_up_camera_t"], r_d["demagalhaes_dn_camera_t"]), xytext=offsets_4a[d_lbl], textcoords="offset points", fontsize=8)
axes[0].set_title("A. RNA-Derived Aging Priors (de Magalhaes UP vs. DN)")
axes[0].set_xlabel("PRIOR_DEMAGALHAES_UP (CAMERA t)")
axes[0].set_ylabel("PRIOR_DEMAGALHAES_DN (CAMERA t)")
axes[0].legend(title="Global BH q < 0.05", frameon=True)

sns.scatterplot(data=bench_df, x="demagalhaes_opposition_score", y="ox_sen_camera_t", color="#1f77b4", s=50, alpha=0.75, ax=axes[1])
axes[1].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes[1].axvline(0, color="gray", linestyle="--", linewidth=0.8)
for d_lbl in ["Binimetinib", "Baicalein", "DFO", "Erlotinib (Hydrochloride)"]:
    r_d = bench_df[bench_df["compound_id"] == d_lbl].iloc[0]
    axes[1].annotate(d_lbl, (r_d["demagalhaes_opposition_score"], r_d["ox_sen_camera_t"]), xytext=(5, 4), textcoords="offset points", fontsize=8)
axes[1].set_title("B. Aging Opposition Score vs. Oxidative Senescence")
axes[1].set_xlabel("de Magalhaes Opposition Score (DN - UP CAMERA t)")
axes[1].set_ylabel("Reactome Oxidative Senescence (CAMERA t)")

plt.tight_layout()
save_fig(fig, "fig4_independent_aging_benchmarks")
bench_df.to_csv("reanalysis/A08/source_tables/fig4_benchmarks.tsv", sep="\t", index=False)

# ==============================================================================
# FIGURE 5: Candidate mechanisms (OBSERVED log2 abundance, NOT DMSO-standardized)
# ==============================================================================
sm = pd.read_csv("reanalysis/A01/sample_manifest.tsv", sep="\t")
prot_map = pd.read_csv("reanalysis/A01/protein_mapping.tsv", sep="\t")
ret_prots = prot_map.loc[prot_map["retained_50pct_filter"], "feature_id"].tolist()
adata = ad.read_h5ad("data/raw (1).h5ad")[:, ret_prots].copy()
X_raw = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.copy()
with np.errstate(divide="ignore"):
    X_log = np.where(X_raw > 0, np.log2(X_raw), np.nan)
meds = np.nanmedian(X_log, axis=1, keepdims=True)
X_norm = X_log - meds + np.nanmedian(meds)

focus_drugs = ["DMSO", "Binimetinib", "Baicalein", "Mobocertinib (succinate)", "Erlotinib (Hydrochloride)"]
ox_sentinels = [g for g in ["NDUFS1", "IDH2", "ACO2", "SDHA", "COX4I1", "ATP5F1A"] if g in ret_prots]
inf_sentinels = [g for g in ["ICAM1", "LDLR", "NFKB1", "RELA", "STAT3", "ATP2B1"] if g in ret_prots]

fig, axes = plt.subplots(2, 1, figsize=(14, 10.5))
pal = {"DMSO": "#7f7f7f", "Binimetinib": "#2ca02c", "Baicalein": "#9467bd", "Mobocertinib (succinate)": "#1f77b4", "Erlotinib (Hydrochloride)": "#d62728"}

for ax, plist, title in [
    (axes[0], ox_sentinels, "A. Prespecified Respiratory/OXPHOS Sentinels (Observed Sample-Median Normalized Log2 Abundance)"),
    (axes[1], inf_sentinels, "B. Prespecified Inflammatory Sentinels (Observed Sample-Median Normalized Log2 Abundance)"),
]:
    m_s = sm["compound_id"].isin(focus_drugs).values & (~sm["excluded"].values) & (sm["plate"].values != "RedMix")
    p_idx = [ret_prots.index(g) for g in plist]
    df_w = pd.DataFrame(X_norm[np.ix_(m_s, p_idx)], columns=plist)
    df_w["Drug"] = sm.loc[m_s, "compound_id"].values
    df_l = df_w.melt(id_vars="Drug", var_name="Protein", value_name="Observed Log2 Abundance")
    df_l["Drug"] = pd.Categorical(df_l["Drug"], categories=focus_drugs, ordered=True)
    sns.violinplot(data=df_l, x="Protein", y="Observed Log2 Abundance", hue="Drug", palette=pal, inner=None, cut=0, alpha=0.35, ax=ax)
    sns.stripplot(data=df_l, x="Protein", y="Observed Log2 Abundance", hue="Drug", palette=pal, dodge=True, size=4.2, alpha=0.9, jitter=0.1, ax=ax)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[:len(focus_drugs)], l[:len(focus_drugs)], title="Treatment (n=6 plates; DMSO n=23)", loc="lower right", frameon=True, ncol=2)
    ax.set_title(title)
    ax.set_ylabel("Observed Log2 Abundance")

plt.tight_layout()
save_fig(fig, "fig5_candidate_mechanisms")

# ==============================================================================
# SUPPLEMENTARY FIGURES S1 & S2
# ==============================================================================
# S1: Full 161-compound heatmap
fig_s1 = plt.figure(figsize=(15, 28))
ax_s1 = fig_s1.add_subplot(111)
annot_s1 = piv_q.map(lambda q: "*" if q < 0.05 else "")
sns.heatmap(piv_t, cmap="RdBu_r", center=0, vmin=-6, vmax=6, annot=annot_s1, fmt="", annot_kws={"size": 6.5},
            linewidths=0.2, linecolor="#eeeeee", ax=ax_s1, cbar_kws={"label": "CAMERA t-statistic"})
ax_s1.set_title("Supplementary Figure S1: Complete 161 Directly Estimable Compound x 14 Program Atlas (* Global BH q < 0.05)")
ax_s1.set_xticklabels(ax_s1.get_xticklabels(), rotation=40, ha="right", fontsize=9.5)
ax_s1.set_yticklabels(ax_s1.get_yticklabels(), fontsize=7)
save_fig(fig_s1, "figS1_full_161_compound_atlas")

# S2: Sensitivity (Primary CAMERA t vs ComBat+Imputation NES & Composition Check)
fig_s2, axes_s2 = plt.subplots(1, 2, figsize=(14, 6))
sns.scatterplot(data=sens_df, x="nes_combat_imp", y="camera_t", alpha=0.45, s=22, color="#1f77b4", ax=axes_s2[0])
axes_s2[0].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes_s2[0].axvline(0, color="gray", linestyle="--", linewidth=0.8)
axes_s2[0].set_xlabel("Preliminary ComBat + Imputation GSEA NES")
axes_s2[0].set_ylabel("Primary Plate-Blocked Observed CAMERA t")
axes_s2[0].set_title("A. Sensitivity: Primary Model vs. ComBat + Imputation")

sns.scatterplot(data=comp_sens, x="cytosolic_ribo_mean_log2fc", y="oxphos_mean_log2fc", alpha=0.65, s=35, color="#2ca02c", ax=axes_s2[1])
axes_s2[1].axhline(0, color="gray", linestyle="--", linewidth=0.8)
axes_s2[1].axvline(0, color="gray", linestyle="--", linewidth=0.8)
for d_lbl in ["Binimetinib", "Baicalein", "DFO"]:
    r_d = comp_sens[comp_sens["compound_id"] == d_lbl].iloc[0]
    axes_s2[1].annotate(d_lbl, (r_d["cytosolic_ribo_mean_log2fc"], r_d["oxphos_mean_log2fc"]), xytext=(5, 4), textcoords="offset points", fontsize=8.5, fontweight="bold")
axes_s2[1].set_xlabel("Cytosolic Ribosome Mean log2FC (RPL/RPS)")
axes_s2[1].set_ylabel("OXPHOS Mean log2FC")
axes_s2[1].set_title("B. Relative-Composition Check: OXPHOS vs. Cytosolic Ribosomes")
plt.tight_layout()
save_fig(fig_s2, "figS2_sensitivity_and_composition")

# ==============================================================================
# CLAIM-TO-RESULT LEDGER & RESULTS SUMMARY (NO EM DASHES)
# ==============================================================================
ledger = [
    {"claim_id": "CLM_01", "claim_text": "161 contrasts (160 library drugs + DFO) are directly estimable across 12 plates; 12 RedMix drugs lack within-plate DMSO.", "source_table": "reanalysis/A01/contrast_eligibility.tsv", "figure_panel": "Fig 1A", "statistical_family": "Design audit", "limitations": "RedMix excluded from primary inference", "manuscript_location": "Results 3.1"},
    {"claim_id": "CLM_02", "claim_text": "Null calibration (DMSO split A vs B across 11 plates) yields 0 BH-significant proteins and 0 BH-significant aging programs.", "source_table": "reanalysis/A04/null_calibration.tsv", "figure_panel": "Fig 1D", "statistical_family": "Paired t + CAMERA BH FDR (n=14)", "limitations": "11 plates with 2 DMSOs", "manuscript_location": "Results 3.1"},
    {"claim_id": "CLM_03", "claim_text": "Under correlation-aware CAMERA and global BH FDR (2,254 tests), 218/2,254 pairs across 142/161 compounds are significant at q < 0.05.", "source_table": "reanalysis/A02/program_effects.tsv", "figure_panel": "Fig 3A-B, Fig S1", "statistical_family": "Global BH FDR (161 x 14 = 2,254)", "limitations": "Relative proteome composition only", "manuscript_location": "Results 3.3"},
    {"claim_id": "CLM_04", "claim_text": "Binimetinib and baicalein remain outliers for OXPHOS elevation and inflammatory suppression after adjusting for proliferation and perturbation size.", "source_table": "reanalysis/A06/candidate_specificity.tsv", "figure_panel": "Fig 3C, Fig 5", "statistical_family": "OLS residual Z-score", "limitations": "No measured cell counts/viability; single MEK inhibitor", "manuscript_location": "Results 3.4-3.5"},
    {"claim_id": "CLM_05", "claim_text": "Under correlation-aware testing, 1/161 (de Magalhaes UP) and 5/161 (de Magalhaes DN) reach global BH q < 0.05; donor-age proteomes await local download.", "source_table": "reanalysis/A05/drug_reference_comparisons.tsv", "figure_panel": "Fig 4", "statistical_family": "Global BH FDR (n=161)", "limitations": "RNA-derived priors; primary donor-age tables not local", "manuscript_location": "Results 3.4"},
]
pd.DataFrame(ledger).to_csv("reanalysis/A08/claim_result_ledger.tsv", sep="\t", index=False)

with open("reanalysis/A08/results_summary.md", "w", encoding="utf-8") as f:
    f.write("""# Reanalysis Results Summary (A01-A08)

- **Title**: Proteomic profiling maps pharmacological responses across aging associated cellular programs
- **Primary Estimable Universe**: 161 contrasts (160 library drugs + DFO) across 12 plates (`Blue1-6`, `Red1-6`); 12 `RedMix` drugs excluded from primary inference due to 0 within-plate `DMSO` controls.
- **Recomputed Significance (Global BH FDR q < 0.05, 2,254 tests)**:
  - Correlation-aware CAMERA: 218 / 2,254 pairs across 142 / 161 compounds.
  - Pre-ranked GSEA (2,000 permutations): 826 / 2,254 pairs across 161 / 161 compounds.
- **Null Calibration**: 0 / 7,872 proteins and 0 / 14 aging programs significant at BH q < 0.05 in within-plate `DMSO` split.
- **Exploratory Candidates**: `Binimetinib` (`PROG_OXPHOS_resid_z = +1.88`, `PROG_INFLAMMATORY_resid_z = -4.15`) and `Baicalein` (`PROG_OXPHOS_resid_z = +1.37`, `PROG_INFLAMMATORY_resid_z = -4.28`) are nominated for orthogonal functional validation.
""")

print("A08 complete: 5 main figures + 2 supplementary figures (PNG+SVG), source tables, ledger, and summary generated.")
