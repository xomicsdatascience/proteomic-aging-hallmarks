import os
import json
import yaml
import time
import platform
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests
import gseapy as gp

os.makedirs("reanalysis/A02", exist_ok=True)
print("=== Running A02: Batch-Aware Moderated Linear Model & Pathway Inference ===")

# Load frozen A01 & A03 manifests
sample_manifest = pd.read_csv("reanalysis/A01/sample_manifest.tsv", sep="\t")
elig_df = pd.read_csv("reanalysis/A01/contrast_eligibility.tsv", sep="\t")
prot_map = pd.read_csv("reanalysis/A01/protein_mapping.tsv", sep="\t")
prog_manifest = pd.read_csv("reanalysis/A03/program_manifest.tsv", sep="\t")

# Load GMT
gene_sets = {}
with open("reanalysis/A03/programs.gmt", "r", encoding="utf-8") as f:
    for line in f:
        parts = line.strip().split("\t")
        gene_sets[parts[0]] = parts[2:]

aging_progs = [p for p in prog_manifest["program_id"] if not p.startswith("NUIS_")]
nuis_progs = [p for p in prog_manifest["program_id"] if p.startswith("NUIS_")]
all_progs = aging_progs + nuis_progs

# Load raw data & filter to retained 7,871 proteins
adata = ad.read_h5ad("data/raw (1).h5ad")
retained_prots = prot_map.loc[prot_map["retained_50pct_filter"], "feature_id"].tolist()
adata = adata[:, retained_prots].copy()

X_raw = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.copy()
with np.errstate(divide="ignore"):
    X_log = np.where(X_raw > 0, np.log2(X_raw), np.nan)

# Sample-wise median normalization (relative proteome composition)
meds = np.nanmedian(X_log, axis=1, keepdims=True)
X_norm = X_log - meds + np.nanmedian(meds)

# ==============================================================================
# 1. PRIMARY MODEL: Plate-blocked observed-data OLS + Empirical Bayes (161 contrasts)
# ==============================================================================
direct_compounds = elig_df.loc[elig_df["eligibility_status"] == "directly_estimable", "compound_id"].tolist()
prim_mask = (
    (~sample_manifest["excluded"].values)
    & (sample_manifest["plate"].values != "RedMix")
    & (sample_manifest["condition_role"].values != "analytical_qc")
)
sm_prim = sample_manifest.loc[prim_mask].reset_index(drop=True)
Y_prim = X_norm[prim_mask, :]  # (n_samples, n_prots)

# Design matrix: DMSO is reference level; sum-to-zero or Blue1 reference for plate
# Use pandas get_dummies with drop_first=False for compound, and plate effects (Blue2..6, Red1..6)
# Note: Blue1..6 and Red1..6 are disjoint sets of compounds!
# Within Blue plates: Blue1..6; within Red plates: Red1..6.
# Fitting ~ 1 + C(compound, Treatment('DMSO')) + C(plate)
# Wait: since DMSO is on ALL 12 plates, DMSO connects Blue and Red plates!
# So the full design matrix [Intercept, compound_dummies (161), plate_dummies (11)] is FULL RANK!
plates_sorted = sorted(sm_prim["plate"].unique())
comp_cols = direct_compounds  # 161 contrasts vs DMSO
plate_cols = plates_sorted[1:]  # 11 plate dummies (Blue1 is ref)

n_samp = len(sm_prim)
p_design = 1 + len(comp_cols) + len(plate_cols)
D = np.zeros((n_samp, p_design), dtype=float)
D[:, 0] = 1.0  # Intercept = DMSO on Blue1
comp_idx = {c: i + 1 for i, c in enumerate(comp_cols)}
plate_idx = {pl: 1 + len(comp_cols) + i for i, pl in enumerate(plate_cols)}

for row_i, r in sm_prim.iterrows():
    if r["compound_id"] in comp_idx:
        D[row_i, comp_idx[r["compound_id"]]] = 1.0
    if r["plate"] in plate_idx:
        D[row_i, plate_idx[r["plate"]]] = 1.0

rank_D = np.linalg.matrix_rank(D)
print(f"Primary design matrix D: shape={D.shape}, rank={rank_D} (full rank={rank_D == p_design})")

# Fit observed-data OLS per protein
n_prots = Y_prim.shape[1]
betas = np.full((n_prots, len(comp_cols)), np.nan)
se_unmod = np.full((n_prots, len(comp_cols)), np.nan)
df_resid = np.zeros(n_prots)
s2_resid = np.full(n_prots, np.nan)
resid_mat = np.full_like(Y_prim, np.nan)

# Fast vectorized path for proteins with 0 missingness in primary samples, fallback for proteins with NaNs
no_nan_mask = ~np.isnan(Y_prim).any(axis=0)
print(f"Proteins with complete observations in primary set: {no_nan_mask.sum()} / {n_prots}")

# Complete-case vectorized OLS
DtD_inv = np.linalg.inv(D.T @ D)
P_proj = DtD_inv @ D.T
Y_cc = Y_prim[:, no_nan_mask]
coef_cc = P_proj @ Y_cc  # (p_design, n_cc)
res_cc = Y_cc - D @ coef_cc
df_cc = n_samp - p_design
s2_cc = (res_cc ** 2).sum(axis=0) / df_cc
diag_inv_cc = np.diag(DtD_inv)[1:1 + len(comp_cols)]

betas[no_nan_mask, :] = coef_cc[1:1 + len(comp_cols), :].T
se_unmod[no_nan_mask, :] = np.sqrt(np.outer(s2_cc, diag_inv_cc))
df_resid[no_nan_mask] = df_cc
s2_resid[no_nan_mask] = s2_cc
resid_mat[:, no_nan_mask] = res_cc

# Per-protein OLS for proteins with some missing values
nan_indices = np.where(~no_nan_mask)[0]
for j in nan_indices:
    y = Y_prim[:, j]
    obs_m = ~np.isnan(y)
    if obs_m.sum() <= p_design + 5:
        continue
    Dj = D[obs_m, :]
    yj = y[obs_m]
    # Drop any columns with all zeros (if a compound had 0 obs for this protein)
    col_ok = (Dj != 0).any(axis=0)
    Dj_sub = Dj[:, col_ok]
    if np.linalg.matrix_rank(Dj_sub) < Dj_sub.shape[1]:
        continue
    DtD_j_inv = np.linalg.inv(Dj_sub.T @ Dj_sub)
    coef_j = DtD_j_inv @ (Dj_sub.T @ yj)
    res_j = yj - Dj_sub @ coef_j
    df_j = len(yj) - Dj_sub.shape[1]
    s2_j = (res_j ** 2).sum() / df_j
    df_resid[j] = df_j
    s2_resid[j] = s2_j
    resid_mat[obs_m, j] = res_j
    # Map back to comp_cols
    sub_col_map = np.where(col_ok)[0]
    diag_j = np.diag(DtD_j_inv)
    for k_sub, k_full in enumerate(sub_col_map):
        if 1 <= k_full <= len(comp_cols):
            c_i = k_full - 1
            betas[j, c_i] = coef_j[k_sub]
            se_unmod[j, c_i] = np.sqrt(s2_j * diag_j[k_sub])

# Empirical Bayes variance moderation (Smyth 2004 method of moments on log(s2))
valid_s2 = (~np.isnan(s2_resid)) & (s2_resid > 0) & (df_resid > 5)
z = np.log(s2_resid[valid_s2])
df_v = df_resid[valid_s2]
from scipy.special import digamma, polygamma
e_z = z - digamma(df_v / 2.0) + np.log(df_v / 2.0)
var_e = np.var(e_z, ddof=1)
mean_trigamma = np.mean(polygamma(1, df_v / 2.0))
target_trigamma = max(var_e - mean_trigamma, 1e-4)

# Solve polygamma(1, d0/2) = target_trigamma via grid + binary search
lo, hi = 0.1, 200.0
for _ in range(50):
    mid = 0.5 * (lo + hi)
    if polygamma(1, mid / 2.0) > target_trigamma:
        lo = mid
    else:
        hi = mid
d0 = 0.5 * (lo + hi)
s0_2 = np.exp(np.mean(e_z) + digamma(d0 / 2.0) - np.log(d0 / 2.0))
print(f"Empirical Bayes prior: d0 = {d0:.2f}, s0^2 = {s0_2:.4f}")

# Moderated variance & t-statistics
s2_post = (df_resid[:, None] * s2_resid[:, None] + d0 * s0_2) / (df_resid[:, None] + d0)
# Scale unmoderated SE by sqrt(s2_post / s2_resid)
se_mod = se_unmod * np.sqrt(s2_post / s2_resid[:, None])
t_mod = betas / se_mod
df_tot = df_resid[:, None] + d0
p_mod = 2.0 * stats.t.sf(np.abs(t_mod), df=df_tot)

# Estimate inter-protein correlation per gene set for CAMERA-style competitive test
# Use complete-case residual matrix or pairwise correlation of residuals
res_df = pd.DataFrame(resid_mat, columns=retained_prots)
vif_map = {}
for pid in all_progs:
    gs = [g for g in gene_sets[pid] if g in retained_prots]
    if len(gs) >= 5:
        sub_r = res_df[gs].dropna(axis=0)
        if len(sub_r) > 20:
            corr_m = np.corrcoef(sub_r.values, rowvar=False)
            m_g = corr_m.shape[0]
            rho_bar = (corr_m.sum() - m_g) / (m_g * (m_g - 1))
            vif = max(1.0 + (m_g - 1) * rho_bar, 1.0)
        else:
            rho_bar, vif = 0.0, 1.0
    else:
        rho_bar, vif = 0.0, 1.0
    vif_map[pid] = (len(gs), rho_bar, vif)
print("CAMERA VIF estimates computed for all 17 programs.")

# Save protein-level effects (compact summary + full table for primary contrasts)
# Also compute global protein-level BH FDR and per-contrast BH FDR
dmso_mask = (sm_prim["compound_id"] == "DMSO").values
prot_eff_rows = []
prog_eff_rows = []

t0 = time.time()
for c_i, comp in enumerate(comp_cols):
    c_mask = (sm_prim["compound_id"] == comp).values
    n_c = int(c_mask.sum())
    t_vec = t_mod[:, c_i]
    b_vec = betas[:, c_i]
    se_vec = se_mod[:, c_i]
    p_vec = p_mod[:, c_i]
    valid_p = ~np.isnan(p_vec)
    q_per_contrast = np.full(n_prots, np.nan)
    if valid_p.any():
        _, q_per_contrast[valid_p], _, _ = multipletests(p_vec[valid_p], method="fdr_bh")

    # Save protein effects for significant or all? 161 * 7871 = 1.26M rows (~100 MB TSV) — totally fine to save!
    df_c = pd.DataFrame({
        "pipeline": "primary_plate_blocked_eb",
        "compound_id": comp,
        "contrast_id": f"{comp}_vs_DMSO",
        "feature_id": retained_prots,
        "n_treated_obs": (~np.isnan(Y_prim[c_mask, :])).sum(axis=0),
        "n_dmso_obs": (~np.isnan(Y_prim[dmso_mask, :])).sum(axis=0),
        "log2fc": np.round(b_vec, 5),
        "se_mod": np.round(se_vec, 5),
        "t_mod": np.round(t_vec, 4),
        "p_nominal": p_vec,
        "q_BH_contrast": q_per_contrast,
    })
    prot_eff_rows.append(df_c)

    # CAMERA-style competitive test on t_vec (using non-NaN proteins)
    t_valid = t_vec[valid_p]
    prots_valid = np.array(retained_prots)[valid_p]
    t_ser = pd.Series(t_valid, index=prots_valid)

    for pid in all_progs:
        gs_in = [g for g in gene_sets[pid] if g in t_ser.index]
        m_in = len(gs_in)
        if m_in < 5:
            continue
        t_in = t_ser.loc[gs_in].values
        t_out = t_ser.drop(gs_in).values
        m_out = len(t_out)
        n_tot = m_in + m_out
        delta = t_in.mean() - t_out.mean()
        # Pooled variance of t-stats
        s_pool2 = ((m_in - 1) * t_in.var(ddof=1) + (m_out - 1) * t_out.var(ddof=1)) / (n_tot - 2)
        _, rho_bar, vif = vif_map[pid]
        se_cam_emp = np.sqrt(s_pool2 * (vif / m_in + 1.0 / m_out))
        t_cam_emp = delta / se_cam_emp if se_cam_emp > 0 else 0.0
        p_cam_emp = 2.0 * stats.t.sf(np.abs(t_cam_emp), df=n_tot - 2)

        # Standard regularized limma::camera (Wu & Smyth 2012, inter.gene.cor = 0.01)
        rho_std = 0.01
        vif_std = 1.0 + rho_std * (m_in - 1.0)
        se_cam = np.sqrt(s_pool2 * (vif_std / m_in + 1.0 / m_out))
        t_cam = delta / se_cam if se_cam > 0 else 0.0
        p_cam = 2.0 * stats.t.sf(np.abs(t_cam), df=n_tot - 2)

        prog_eff_rows.append({
            "pipeline": "primary_plate_blocked_eb",
            "compound_id": comp,
            "contrast_id": f"{comp}_vs_DMSO",
            "program_id": pid,
            "program_family": "nuisance_control" if pid.startswith("NUIS_") else "primary_aging",
            "n_genes_measured": m_in,
            "rho_inter_protein": round(rho_bar, 4),
            "vif": round(vif, 3),
            "vif_emp": round(vif, 3),
            "camera_t_emp": round(t_cam_emp, 4),
            "p_camera_emp": p_cam_emp,
            "rho_camera_std": rho_std,
            "vif_camera_std": round(vif_std, 3),
            "camera_t": round(t_cam, 4),
            "p_camera": p_cam,
            "direction": "up" if t_cam > 0 else "down",
        })

    # Pre-ranked GSEA (2000 permutations for better resolution)
    rnk = t_ser.sort_values(ascending=False)
    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets={k: gene_sets[k] for k in all_progs},
        min_size=5,
        max_size=500,
        permutation_num=2000,
        seed=42,
        threads=4,
        verbose=False,
    )
    res_g = pre_res.res2d.set_index("Term")
    for row in prog_eff_rows[-len(all_progs):]:
        pid = row["program_id"]
        if pid in res_g.index:
            row["nes_gsea"] = float(res_g.loc[pid, "NES"])
            # Finite-sampling nonzero convention: max(p, 1/(2000+1))
            p_nom_g = float(res_g.loc[pid, "NOM p-val"])
            row["p_gsea_nominal"] = max(p_nom_g, 1.0 / 2001.0)
            row["q_gsea_native"] = float(res_g.loc[pid, "FDR q-val"])
            row["lead_genes"] = str(res_g.loc[pid, "Lead_genes"])
        else:
            row["nes_gsea"] = np.nan
            row["p_gsea_nominal"] = np.nan
            row["q_gsea_native"] = np.nan
            row["lead_genes"] = ""

    if (c_i + 1) % 25 == 0 or (c_i + 1) == len(comp_cols):
        print(f"   Primary A02: {c_i + 1}/{len(comp_cols)} contrasts ({time.time() - t0:.1f}s)")

df_prot_all = pd.concat(prot_eff_rows, ignore_index=True)
df_prot_all.to_csv("reanalysis/A02/protein_effects.tsv", sep="\t", index=False)

df_prog = pd.DataFrame(prog_eff_rows)

# Compute global BH FDR across the PRIMARY family (161 contrasts x 14 aging programs)
prim_fam_mask = df_prog["program_family"] == "primary_aging"
df_prog["q_BH_global_camera"] = np.nan
df_prog["q_BH_global_gsea"] = np.nan
_, df_prog.loc[prim_fam_mask, "q_BH_global_camera"], _, _ = multipletests(
    df_prog.loc[prim_fam_mask, "p_camera"], method="fdr_bh"
)
_, df_prog.loc[prim_fam_mask, "q_BH_global_gsea"], _, _ = multipletests(
    df_prog.loc[prim_fam_mask, "p_gsea_nominal"], method="fdr_bh"
)

# Also compute BH FDR for the secondary nuisance family
nuis_mask = df_prog["program_family"] == "nuisance_control"
_, df_prog.loc[nuis_mask, "q_BH_global_camera"], _, _ = multipletests(
    df_prog.loc[nuis_mask, "p_camera"], method="fdr_bh"
)
_, df_prog.loc[nuis_mask, "q_BH_global_gsea"], _, _ = multipletests(
    df_prog.loc[nuis_mask, "p_gsea_nominal"], method="fdr_bh"
)

df_prog.to_csv("reanalysis/A02/program_effects.tsv", sep="\t", index=False)

# Save config & session info
with open("reanalysis/A02/contrast_definitions.json", "w") as f:
    json.dump({"primary_contrasts": [f"{c}_vs_DMSO" for c in comp_cols], "reference": "DMSO (n=23 across 12 plates)"}, f, indent=2)

config = {
    "primary_model": "plate_blocked_observed_OLS_empirical_Bayes",
    "normalization": "sample_median_log2_relative_composition",
    "missingness_handling": "observed_data_only_no_imputation",
    "d0_prior": float(d0),
    "s0_2_prior": float(s0_2),
    "gsea_permutations": 2000,
    "gsea_nonzero_floor": 1.0 / 2001.0,
    "primary_bh_family_size": int(prim_fam_mask.sum()),
}
with open("reanalysis/A02/analysis_config.yaml", "w") as f:
    yaml.dump(config, f, sort_keys=False)

with open("reanalysis/A02/session_info.txt", "w") as f:
    f.write(f"Python: {platform.python_version()}\nOS: {platform.platform()}\n")
    for mod in [np, pd, ad, sc, gp]:
        f.write(f"{mod.__name__}: {getattr(mod, '__version__', 'unknown')}\n")

# Summary stats
n_sig_cam = (df_prog.loc[prim_fam_mask, "q_BH_global_camera"] < 0.05).sum()
n_sig_gsea = (df_prog.loc[prim_fam_mask, "q_BH_global_gsea"] < 0.05).sum()
n_drugs_cam = df_prog.loc[prim_fam_mask & (df_prog["q_BH_global_camera"] < 0.05), "compound_id"].nunique()
n_drugs_gsea = df_prog.loc[prim_fam_mask & (df_prog["q_BH_global_gsea"] < 0.05), "compound_id"].nunique()

rep_md = f"""# A02 Inference Report

- **Primary Family**: `{prim_fam_mask.sum()}` hypotheses (`{len(comp_cols)}` directly estimable contrasts x `14` primary aging programs). `12` `RedMix` contrasts were excluded from primary inference due to `0` within-plate `DMSO` controls.
- **Model**: Plate-blocked observed-data OLS (`~ 1 + compound + plate`) with empirical Bayes variance moderation (`d0 = {d0:.2f}`, `s0^2 = {s0_2:.4f}`). Normalization is sample-median $\\log_2$, measuring **relative proteome composition** rather than absolute per-cell protein mass.
- **Recomputed Significance Counts (Global BH FDR $q < 0.05$)**:
  - **Correlation-aware competitive test (CAMERA)**: `{n_sig_cam}` / `{prim_fam_mask.sum()}` pairs across `{n_drugs_cam}` / `{len(comp_cols)}` compounds.
  - **Pre-Ranked GSEA (`q_BH_global_gsea < 0.05`)**: `{n_sig_gsea}` / `{prim_fam_mask.sum()}` pairs across `{n_drugs_gsea}` / `{len(comp_cols)}` compounds.
  - *(Compare to preliminary unblocked ComBat + native-q report of 165/173 drugs and 999 pairs.)*
"""
with open("reanalysis/A02/inference_report.md", "w", encoding="utf-8") as f:
    f.write(rep_md)

print(rep_md)
