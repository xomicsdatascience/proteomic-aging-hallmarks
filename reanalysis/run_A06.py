import os
import numpy as np
import pandas as pd
from scipy import stats

os.makedirs("reanalysis/A06", exist_ok=True)
print("=== Running A06: Shared Stress, Growth & Relative Composition Specificity ===")

prog_df = pd.read_csv("reanalysis/A02/program_effects.tsv", sep="\t")
prot_df = pd.read_csv("reanalysis/A02/protein_effects.tsv", sep="\t")

# 1. Shared response scores per compound (nuisance proliferation, apoptosis, global perturbation)
piv_cam = prog_df.pivot(index="compound_id", columns="program_id", values="camera_t")

pert_stats = prot_df.groupby("compound_id").agg(
    mean_abs_t=("t_mod", lambda s: np.nanmean(np.abs(s))),
    mean_abs_log2fc=("log2fc", lambda s: np.nanmean(np.abs(s))),
    n_de_q05=("q_BH_contrast", lambda s: int((s < 0.05).sum())),
).reset_index()

shared = pert_stats.merge(
    piv_cam[["NUIS_E2F_PROLIF", "NUIS_G2M_PROLIF", "NUIS_APOPTOSIS",
             "PROG_OXPHOS", "PROG_INFLAMMATORY", "PROG_MTORC1_SIGNALING",
             "PROG_CHROMATIN_ORG", "PROG_MACROAUTOPHAGY"]].reset_index(),
    on="compound_id",
)
shared.to_csv("reanalysis/A06/shared_response_scores.tsv", sep="\t", index=False)

# 2. Composition sensitivity: compare mean log2FC of OXPHOS vs cytosolic ribosomal proteins (RPL/RPS)
ribo_prots = [p for p in prot_df["feature_id"].unique() if p.startswith(("RPL", "RPS"))]
with open("reanalysis/A03/programs.gmt") as f:
    gmt = {l.split("\t")[0]: set(l.strip().split("\t")[2:]) for l in f}
ox_set = gmt["PROG_OXPHOS"]

comp_rows = []
for comp, grp in prot_df.groupby("compound_id"):
    ox_fc = grp.loc[grp["feature_id"].isin(ox_set), "log2fc"].mean()
    ribo_fc = grp.loc[grp["feature_id"].isin(ribo_prots), "log2fc"].mean()
    non_ox_fc = grp.loc[~grp["feature_id"].isin(ox_set), "log2fc"].mean()
    comp_rows.append({
        "compound_id": comp,
        "oxphos_mean_log2fc": round(float(ox_fc), 5),
        "cytosolic_ribo_mean_log2fc": round(float(ribo_fc), 5),
        "non_oxphos_mean_log2fc": round(float(non_ox_fc), 5),
        "oxphos_minus_ribo_log2fc": round(float(ox_fc - ribo_fc), 5),
    })
comp_df = pd.DataFrame(comp_rows)
comp_df.to_csv("reanalysis/A06/composition_sensitivity.tsv", sep="\t", index=False)

# 3. Specificity analysis: regress PROG_OXPHOS and PROG_INFLAMMATORY on NUIS_E2F_PROLIF + mean_abs_t
spec = shared.copy()
X_cov = np.column_stack([np.ones(len(spec)), spec["NUIS_E2F_PROLIF"].values, spec["mean_abs_t"].values])

for target in ["PROG_OXPHOS", "PROG_INFLAMMATORY", "PROG_MTORC1_SIGNALING", "PROG_MACROAUTOPHAGY"]:
    y = spec[target].values
    beta, *_ = np.linalg.lstsq(X_cov, y, rcond=None)
    resid = y - X_cov @ beta
    spec[f"{target}_resid_z"] = (resid - resid.mean()) / resid.std(ddof=1)

# 4. 14-Pathway Composite Geroprotective Score (CGS) taking into account desired direction for reduced aging
from statsmodels.stats.multitest import multipletests

signs_14 = {
    "PROG_OXPHOS": +1,
    "PRIOR_DEMAGALHAES_DN": +1,
    "PROG_MACROAUTOPHAGY": +1,
    "PROG_PROTEOSTASIS_UPR": +1,
    "PROG_PROTEASOME": +1,
    "PROG_DNA_REPAIR": +1,
    "PROG_TELOMERE_MAINT": +1,
    "PROG_CHROMATIN_ORG": +1,
    "PROG_INFLAMMATORY": -1,
    "PROG_INNATE_TLR": -1,
    "PROG_SENESCENCE_CORE": -1,
    "PROG_MTORC1_SIGNALING": -1,
    "PRIOR_DEMAGALHAES_UP": -1,
    "PROG_CELL_CELL_COMM": -1,
}
aging_cols = list(signs_14.keys())
s_vec = np.array([signs_14[c] for c in aging_cols], dtype=float)

ret_prots_set = set(prot_df["feature_id"].unique())
Sigma_0 = np.eye(len(aging_cols))
for i_k, k1 in enumerate(aging_cols):
    g1 = gmt[k1] & ret_prots_set
    for j_k, k2 in enumerate(aging_cols):
        if i_k != j_k:
            g2 = gmt[k2] & ret_prots_set
            Sigma_0[i_k, j_k] = len(g1 & g2) / np.sqrt(len(g1) * len(g2))
se_null = np.sqrt(s_vec @ Sigma_0 @ s_vec)

piv_cam_ordered = piv_cam.loc[spec["compound_id"]]
z_unadj = (piv_cam_ordered[aging_cols].values @ s_vec) / se_null
p_unadj = 2.0 * stats.norm.sf(np.abs(z_unadj))
_, q_unadj, _, _ = multipletests(p_unadj, method="fdr_bh")

X_prolif = np.column_stack([
    np.ones(len(piv_cam_ordered)),
    piv_cam_ordered["NUIS_E2F_PROLIF"].values,
    piv_cam_ordered["NUIS_G2M_PROLIF"].values,
])
T_adj = pd.DataFrame(index=piv_cam_ordered.index, columns=aging_cols, dtype=float)
for c in aging_cols:
    y = piv_cam_ordered[c].values
    b, *_ = np.linalg.lstsq(X_prolif, y, rcond=None)
    T_adj[c] = y - X_prolif @ b

z_adj = (T_adj.values @ s_vec) / se_null
p_adj = 2.0 * stats.norm.sf(np.abs(z_adj))
_, q_adj, _, _ = multipletests(p_adj, method="fdr_bh")

z_comb = (z_unadj + z_adj) / np.sqrt(2.0)
p_comb = 2.0 * stats.norm.sf(np.abs(z_comb))
_, q_comb, _, _ = multipletests(p_comb, method="fdr_bh")

piv_nes = prog_df.pivot(index="compound_id", columns="program_id", values="nes_gsea").loc[spec["compound_id"]]
Z_nes = (piv_nes[aging_cols] - piv_nes[aging_cols].mean()) / piv_nes[aging_cols].std(ddof=1)
cgs_nes = (Z_nes.values @ s_vec) / np.sqrt(s_vec @ Z_nes.corr().values @ s_vec)

spec["CGS_unadj"] = np.round(z_unadj, 4)
spec["q_BH_CGS_unadj"] = np.round(q_unadj, 5)
spec["n_beneficial_unadj"] = ((piv_cam_ordered[aging_cols].values * s_vec) > 0).sum(axis=1)
spec["CGS_adj"] = np.round(z_adj, 4)
spec["q_BH_CGS_adj"] = np.round(q_adj, 5)
spec["n_beneficial_adj"] = ((T_adj.values * s_vec) > 0).sum(axis=1)
spec["CGS_combined"] = np.round(z_comb, 4)
spec["q_BH_CGS_combined"] = np.round(q_comb, 5)
spec["CGS_NES_screen"] = np.round(cgs_nes, 4)

spec.to_csv("reanalysis/A06/candidate_specificity.tsv", sep="\t", index=False)

r_ox_e2f, p_ox_e2f = stats.pearsonr(spec["PROG_OXPHOS"], spec["NUIS_E2F_PROLIF"])
r_inf_e2f, p_inf_e2f = stats.pearsonr(spec["PROG_INFLAMMATORY"], spec["NUIS_E2F_PROLIF"])

rep_a06 = f"""# A06 Shared Stress, Growth & Relative Composition Specificity Report

- **Missing Measured Growth Data**: Cell counts, viability, and total protein per cell are not present in `raw (1).h5ad` or `drug_details (1).xlsx`. Cytostasis or toxicity cannot be ruled out; molecular proxies (`NUIS_E2F_PROLIF`, `NUIS_G2M_PROLIF`, `NUIS_APOPTOSIS`) are used for sensitivity only.
- **Correlation with Proliferation Proxy (`NUIS_E2F_PROLIF`)**:
  - `PROG_OXPHOS` vs `NUIS_E2F_PROLIF`: $r = {r_ox_e2f:.3f}$ ($p = {p_ox_e2f:.2e}$)
  - `PROG_INFLAMMATORY` vs `NUIS_E2F_PROLIF`: $r = {r_inf_e2f:.3f}$ ($p = {p_inf_e2f:.2e}$)
- **Candidate Specificity Residuals** (`candidate_specificity.tsv`, Z-scored after adjusting for `NUIS_E2F_PROLIF` and `mean_abs_t`):
  - `Binimetinib`: `PROG_OXPHOS_resid_z = {spec.loc[spec['compound_id']=='Binimetinib', 'PROG_OXPHOS_resid_z'].values[0]:+.2f}`, `PROG_INFLAMMATORY_resid_z = {spec.loc[spec['compound_id']=='Binimetinib', 'PROG_INFLAMMATORY_resid_z'].values[0]:+.2f}`
  - `Baicalein`: `PROG_OXPHOS_resid_z = {spec.loc[spec['compound_id']=='Baicalein', 'PROG_OXPHOS_resid_z'].values[0]:+.2f}`, `PROG_INFLAMMATORY_resid_z = {spec.loc[spec['compound_id']=='Baicalein', 'PROG_INFLAMMATORY_resid_z'].values[0]:+.2f}`
- **Relative-Composition Check**: Sample-median normalization measures relative proteome composition. `composition_sensitivity.tsv` reports `oxphos_mean_log2fc` alongside `cytosolic_ribo_mean_log2fc`.
"""
with open("reanalysis/A06/specificity_report.md", "w", encoding="utf-8") as f:
    f.write(rep_a06)
print(rep_a06)
