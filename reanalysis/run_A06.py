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
