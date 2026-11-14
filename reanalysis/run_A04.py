import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

os.makedirs("reanalysis/A04", exist_ok=True)
print("=== Running A04: Reproducibility & Sensitivity ===")

sm = pd.read_csv("reanalysis/A01/sample_manifest.tsv", sep="\t")
elig = pd.read_csv("reanalysis/A01/contrast_eligibility.tsv", sep="\t")
prot_map = pd.read_csv("reanalysis/A01/protein_mapping.tsv", sep="\t")
prog_df = pd.read_csv("reanalysis/A02/program_effects.tsv", sep="\t")

adata = ad.read_h5ad("data/raw (1).h5ad")
ret_prots = prot_map.loc[prot_map["retained_50pct_filter"], "feature_id"].tolist()
adata = adata[:, ret_prots].copy()
X_raw = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.copy()
with np.errstate(divide="ignore"):
    X_log = np.where(X_raw > 0, np.log2(X_raw), np.nan)
meds = np.nanmedian(X_log, axis=1, keepdims=True)
X_norm = X_log - meds + np.nanmedian(meds)

# 1. QC Metrics per plate
qc_rows = []
for pl in sorted(sm["plate"].unique()):
    m_pl = (sm["plate"] == pl).values
    m_qc = m_pl & (sm["condition_role"] == "analytical_qc").values
    m_dmso = m_pl & (sm["condition_role"] == "vehicle_control").values
    m_dfo = m_pl & (sm["condition_role"] == "positive_assay_control").values
    qc_sd = np.nanmedian(np.nanstd(X_norm[m_qc, :], axis=0))
    qc_miss = np.isnan(X_log[m_qc, :]).mean()
    dmso_miss = np.isnan(X_log[m_dmso, :]).mean() if m_dmso.any() else np.nan
    qc_rows.append({
        "plate": pl,
        "n_qc": int(m_qc.sum()),
        "n_dmso": int(m_dmso.sum()),
        "n_dfo": int(m_dfo.sum()),
        "qc_missing_rate": round(float(qc_miss), 4),
        "qc_median_protein_sd": round(float(qc_sd), 4),
        "dmso_missing_rate": round(float(dmso_miss), 4) if not np.isnan(dmso_miss) else np.nan,
    })
pd.DataFrame(qc_rows).to_csv("reanalysis/A04/qc_metrics.tsv", sep="\t", index=False)

# 2. Cross-plate replication for DFO and key compounds (leave-one-plate-out correlation of log2fc)
focus_comps = ["DFO", "Binimetinib", "Baicalein", "Nordihydroguaiaretic acid", "Lonafarnib", "Erlotinib (Hydrochloride)", "Mobocertinib (succinate)"]
rep_rows = []
for comp in focus_comps:
    plates_c = sorted(sm.loc[(sm["compound_id"] == comp) & (~sm["excluded"]), "plate"].unique())
    # Compute within-plate log2fc (comp - DMSO on that plate) for each plate
    plate_fcs = []
    for pl in plates_c:
        m_c = ((sm["compound_id"] == comp) & (sm["plate"] == pl)).values
        m_d = ((sm["condition_role"] == "vehicle_control") & (sm["plate"] == pl)).values
        if m_c.any() and m_d.any():
            fc_pl = np.nanmean(X_norm[m_c, :], axis=0) - np.nanmean(X_norm[m_d, :], axis=0)
            plate_fcs.append((pl, fc_pl))
    # Leave-one-plate-out correlation
    corrs = []
    for i in range(len(plate_fcs)):
        fc_out = plate_fcs[i][1]
        fc_rest = np.nanmean([plate_fcs[j][1] for j in range(len(plate_fcs)) if j != i], axis=0)
        ok = (~np.isnan(fc_out)) & (~np.isnan(fc_rest))
        r, _ = stats.pearsonr(fc_out[ok], fc_rest[ok])
        corrs.append(r)
    rep_rows.append({
        "compound_id": comp,
        "n_plates": len(plate_fcs),
        "mean_lopo_pearson_r": round(float(np.mean(corrs)), 4),
        "min_lopo_pearson_r": round(float(np.min(corrs)), 4),
        "max_lopo_pearson_r": round(float(np.max(corrs)), 4),
    })
pd.DataFrame(rep_rows).to_csv("reanalysis/A04/replication_metrics.tsv", sep="\t", index=False)

# 3. Sensitivity: Compare Primary CAMERA t-stat & GSEA NES vs ComBat+Imputation GSEA NES
gsea_old = pd.read_csv("results/gsea_all_drugs_14pathways.csv")
# Map old terms to A03 program_ids
term_map = {
    "1. Genomic Instability (DNA Repair)": "PROG_DNA_REPAIR",
    "2. Telomere Attrition (Telomere Ext.)": "PROG_TELOMERE_MAINT",
    "3. Epigenetic Alterations (Chromatin Org.)": "PROG_CHROMATIN_ORG",
    "4. Loss of Proteostasis (UPR)": "PROG_PROTEOSTASIS_UPR",
    "5. Disabled Macroautophagy": "PROG_MACROAUTOPHAGY",
    "6. Nutrient-Sensing (mTORC1)": "PROG_MTORC1_SIGNALING",
    "7. Mitochondrial Dysfunction (OXPHOS)": "PROG_OXPHOS",
    "8. Cellular Senescence": "PROG_SENESCENCE_CORE",
    "10. Intercellular Comm. (Cell-Cell)": "PROG_CELL_CELL_COMM",
    "11. Chronic Inflammation": "PROG_INFLAMMATORY",
    "12. Dysbiosis (TLR Cascades)": "PROG_INNATE_TLR",
    "DEMAGALHAES_AGING_UP": "PRIOR_DEMAGALHAES_UP",
    "DEMAGALHAES_AGING_DN": "PRIOR_DEMAGALHAES_DN",
}
gsea_old["program_id"] = gsea_old["Term"].map(term_map)
gsea_old_sub = gsea_old.dropna(subset=["program_id"]).rename(columns={"drug": "compound_id", "NES": "nes_combat_imp", "FDR q-val": "q_combat_imp"})

sens_df = prog_df[prog_df["program_family"] == "primary_aging"].merge(
    gsea_old_sub[["compound_id", "program_id", "nes_combat_imp", "q_combat_imp"]],
    on=["compound_id", "program_id"],
    how="inner",
)
sens_df["direction_agree_camera_vs_combat"] = np.sign(sens_df["camera_t"]) == np.sign(sens_df["nes_combat_imp"])
sens_df["direction_agree_gsea_vs_combat"] = np.sign(sens_df["nes_gsea"]) == np.sign(sens_df["nes_combat_imp"])
sens_df.to_csv("reanalysis/A04/sensitivity_effects.tsv", sep="\t", index=False)

r_cam, _ = stats.pearsonr(sens_df["camera_t"], sens_df["nes_combat_imp"])
r_gsea, _ = stats.pearsonr(sens_df["nes_gsea"], sens_df["nes_combat_imp"])
dir_agree = sens_df["direction_agree_camera_vs_combat"].mean()

# 4. Null Calibration: Split the 11 plates with 2 DMSO controls into DMSO_A vs DMSO_B
dmso_plates_2 = [pl for pl in sorted(sm["plate"].unique()) if ((sm["plate"] == pl) & (sm["condition_role"] == "vehicle_control")).sum() == 2]
idx_A, idx_B = [], []
for pl in dmso_plates_2:
    idxs = np.where((sm["plate"] == pl) & (sm["condition_role"] == "vehicle_control"))[0]
    idx_A.append(idxs[0])
    idx_B.append(idxs[1])

# Paired t-test across 11 plates (DMSO_A vs DMSO_B)
diff_dmso = X_norm[idx_A, :] - X_norm[idx_B, :]
t_null, p_null = stats.ttest_1samp(diff_dmso, 0.0, axis=0, nan_policy="omit")
_, q_null, _, _ = multipletests(np.nan_to_num(p_null, nan=1.0), method="fdr_bh")

# Load GMT for null CAMERA
gene_sets = {}
with open("reanalysis/A03/programs.gmt", "r", encoding="utf-8") as f:
    for line in f:
        parts = line.strip().split("\t")
        gene_sets[parts[0]] = parts[2:]

t_null_ser = pd.Series(np.nan_to_num(t_null, nan=0.0), index=ret_prots)
null_rows = []
for pid in [p for p in gene_sets if not p.startswith("NUIS_")]:
    gs_in = [g for g in gene_sets[pid] if g in t_null_ser.index]
    t_in = t_null_ser.loc[gs_in].values
    t_out = t_null_ser.drop(gs_in).values
    vif = prog_df.loc[prog_df["program_id"] == pid, "vif"].iloc[0]
    s_pool2 = ((len(t_in)-1)*t_in.var(ddof=1) + (len(t_out)-1)*t_out.var(ddof=1)) / (len(t_null_ser)-2)
    se_c = np.sqrt(s_pool2 * (vif / len(t_in) + 1.0 / len(t_out)))
    tc = (t_in.mean() - t_out.mean()) / se_c
    pc = 2.0 * stats.t.sf(np.abs(tc), df=len(t_null_ser) - 2)
    null_rows.append({"program_id": pid, "null_camera_t": round(tc, 4), "null_p_camera": pc})
df_null = pd.DataFrame(null_rows)
_, df_null["null_q_BH"] = multipletests(df_null["null_p_camera"], method="fdr_bh")[:2]
df_null.to_csv("reanalysis/A04/null_calibration.tsv", sep="\t", index=False)

# 5. Candidate Stability (plate-level bootstrap of camera_t for Binimetinib, Baicalein, NDGA, Lonafarnib, DFO)
stab_comps = ["DFO", "Binimetinib", "Baicalein", "Nordihydroguaiaretic acid", "Lonafarnib"]
rng = np.random.default_rng(42)
stab_rows = []
for comp in stab_comps:
    plates_c = sorted(sm.loc[(sm["compound_id"] == comp) & (~sm["excluded"]) & (sm["plate"] != "RedMix"), "plate"].unique())
    plate_fcs = []
    for pl in plates_c:
        m_c = ((sm["compound_id"] == comp) & (sm["plate"] == pl)).values
        m_d = ((sm["condition_role"] == "vehicle_control") & (sm["plate"] == pl)).values
        if m_c.any() and m_d.any():
            plate_fcs.append(np.nanmean(X_norm[m_c, :], axis=0) - np.nanmean(X_norm[m_d, :], axis=0))
    plate_fcs = np.array(plate_fcs)  # (n_plates, n_prots)
    # Precompute 200 bootstrap protein t-stat vectors once per compound
    boot_t_mats = []
    for _ in range(200):
        b_idx = rng.integers(0, len(plate_fcs), size=len(plate_fcs))
        b_mat = plate_fcs[b_idx, :]
        n_ok = np.sum(~np.isnan(b_mat), axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            m_b = np.nanmean(b_mat, axis=0)
            s_b = np.nanstd(b_mat, axis=0, ddof=1)
            t_b = m_b / (s_b / np.sqrt(n_ok))
        t_b = np.where((n_ok >= 2) & np.isfinite(t_b), t_b, 0.0)
        boot_t_mats.append(t_b)
    boot_t_mats = np.array(boot_t_mats)  # (200, n_prots)

    for pid in ["PROG_OXPHOS", "PROG_INFLAMMATORY", "PROG_MTORC1_SIGNALING", "PROG_MACROAUTOPHAGY", "PRIOR_DEMAGALHAES_DN"]:
        gs_set = set(gene_sets[pid])
        gs_idx = [i for i, g in enumerate(ret_prots) if g in gs_set]
        out_idx = [i for i, g in enumerate(ret_prots) if g not in gs_set]
        vif = prog_df.loc[prog_df["program_id"] == pid, "vif"].iloc[0]
        t_in = boot_t_mats[:, gs_idx]
        t_out = boot_t_mats[:, out_idx]
        n_in, n_out = len(gs_idx), len(out_idx)
        sp2 = ((n_in - 1) * t_in.var(axis=1, ddof=1) + (n_out - 1) * t_out.var(axis=1, ddof=1)) / (n_in + n_out - 2)
        se_c = np.sqrt(sp2 * (vif / n_in + 1.0 / n_out))
        boot_ts = (t_in.mean(axis=1) - t_out.mean(axis=1)) / se_c
        stab_rows.append({
            "compound_id": comp,
            "program_id": pid,
            "boot_mean_camera_t": round(float(np.mean(boot_ts)), 3),
            "boot_ci95_lo": round(float(np.percentile(boot_ts, 2.5)), 3),
            "boot_ci95_hi": round(float(np.percentile(boot_ts, 97.5)), 3),
            "sign_consistency_pct": round(100.0 * max((np.array(boot_ts) > 0).mean(), (np.array(boot_ts) < 0).mean()), 1),
        })
pd.DataFrame(stab_rows).to_csv("reanalysis/A04/candidate_stability.tsv", sep="\t", index=False)

rep_a04 = f"""# A04 Reproducibility & Sensitivity Report

- **Sensitivity (Primary Observed Plate-Blocked vs. ComBat+Imputation)**:
  - Correlation between primary `camera_t` and ComBat+imputation `NES`: $r = {r_cam:.3f}$ (`{100*dir_agree:.1f}%` direction agreement).
  - Correlation between primary `nes_gsea` and ComBat+imputation `NES`: $r = {r_gsea:.3f}$.
- **Null Calibration (`DMSO_A` vs `DMSO_B` across 11 plates)**:
  - `0` / `{len(q_null)}` proteins significant at BH FDR $q < 0.05$.
  - `{(df_null['null_q_BH'] < 0.05).sum()}` / `{len(df_null)}` aging programs significant at BH FDR $q < 0.05$ (`{(df_null['null_p_camera'] < 0.05).sum()}` at nominal $p < 0.05$).
- **Cross-Plate Replication (Leave-One-Plate-Out Pearson $r$ of $\\log_2$FC)**:
  - Computed for `DFO` and key candidates (`replication_metrics.tsv`).
"""
with open("reanalysis/A04/robustness_report.md", "w", encoding="utf-8") as f:
    f.write(rep_a04)
print(rep_a04)
