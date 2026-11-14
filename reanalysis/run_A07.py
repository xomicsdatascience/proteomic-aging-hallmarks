import os
import yaml
import numpy as np
import pandas as pd
from scipy import stats

os.makedirs("reanalysis/A07", exist_ok=True)
print("=== Running A07: Candidate Mechanisms & Pharmacology ===")

prog_df = pd.read_csv("reanalysis/A02/program_effects.tsv", sep="\t")
prot_df = pd.read_csv("reanalysis/A02/protein_effects.tsv", sep="\t")
spec_df = pd.read_csv("reanalysis/A06/candidate_specificity.tsv", sep="\t")
comp_df = pd.read_csv("reanalysis/A01/compound_manifest.tsv", sep="\t")
bench_df = pd.read_csv("reanalysis/A05/drug_reference_comparisons.tsv", sep="\t")

# 1. Freeze candidate selection rules
rules = {
    "eligibility": "directly_estimable (within-plate DMSO controls across 5-6 plates)",
    "primary_pathway_threshold": "q_BH_global_camera < 0.05",
    "specificity_criterion": "abs(resid_z) >= 1.5 after adjusting for NUIS_E2F_PROLIF and mean_abs_t",
    "exploratory_candidates": ["Binimetinib", "Baicalein"],
    "literature_comparators": ["Nordihydroguaiaretic acid", "Lonafarnib", "α-Lipoic Acid", "Simvastatin acid (ammonium)", "Atorvastatin (hemicalcium salt)", "DFO"],
    "kinase_comparators": ["Erlotinib (Hydrochloride)", "Mobocertinib (succinate)"],
}
with open("reanalysis/A07/candidate_selection_rules.yaml", "w") as f:
    yaml.dump(rules, f, sort_keys=False)

# 2. Candidate evidence table across all 161 eligible compounds
prim = prog_df[prog_df["program_family"] == "primary_aging"]
sig_cam_counts = prim.groupby("compound_id")["q_BH_global_camera"].apply(lambda s: int((s < 0.05).sum())).rename("n_sig_camera_q05")
sig_gsea_counts = prim.groupby("compound_id")["q_BH_global_gsea"].apply(lambda s: int((s < 0.05).sum())).rename("n_sig_gsea_q05")

ev = spec_df.merge(sig_cam_counts, on="compound_id").merge(sig_gsea_counts, on="compound_id")
ev = ev.merge(comp_df[["compound_id", "Category", "Target", "PathWay"]], on="compound_id", how="left")
ev = ev.merge(bench_df[["compound_id", "demagalhaes_opposition_score", "ox_sen_camera_t", "ox_sen_q_BH"]], on="compound_id", how="left")
ev["passes_specificity_oxphos_or_inflam"] = (ev["PROG_OXPHOS_resid_z"].abs() >= 1.5) | (ev["PROG_INFLAMMATORY_resid_z"].abs() >= 1.5)
ev = ev.sort_values(["CGS_combined", "CGS_adj", "n_sig_camera_q05"], ascending=[False, False, False])
ev.to_csv("reanalysis/A07/candidate_evidence.tsv", sep="\t", index=False)

# 3. Complete prespecified respiratory-complex and inflammatory module effects for key panel
panel_drugs = ["Binimetinib", "Erlotinib (Hydrochloride)", "Mobocertinib (succinate)", "Baicalein", "DFO", "Nordihydroguaiaretic acid", "Lonafarnib"]
with open("reanalysis/A03/programs.gmt") as f:
    gmt = {l.split("\t")[0]: set(l.strip().split("\t")[2:]) for l in f}

mod_genes = sorted((gmt["PROG_OXPHOS"] | gmt["PROG_INFLAMMATORY"]) & set(prot_df["feature_id"].unique()))
sub_prot = prot_df[prot_df["compound_id"].isin(panel_drugs) & prot_df["feature_id"].isin(mod_genes)].copy()
sub_prot["module"] = np.where(sub_prot["feature_id"].isin(gmt["PROG_OXPHOS"]), "PROG_OXPHOS", "PROG_INFLAMMATORY")

# Annotate GSEA leading-edge membership
le_map = {}
for _, r in prim[prim["compound_id"].isin(panel_drugs) & prim["program_id"].isin(["PROG_OXPHOS", "PROG_INFLAMMATORY"])].iterrows():
    le_map[(r["compound_id"], r["program_id"])] = set(str(r["lead_genes"]).split(";"))

sub_prot["in_gsea_leading_edge"] = [
    r["feature_id"] in le_map.get((r["compound_id"], r["module"]), set())
    for _, r in sub_prot.iterrows()
]
sub_prot.to_csv("reanalysis/A07/candidate_protein_effects.tsv", sep="\t", index=False)

# 4. Drug similarity: full proteome vs 14 aging programs
piv_prot = prot_df.pivot(index="compound_id", columns="feature_id", values="t_mod").fillna(0.0)
piv_prog = prim.pivot(index="compound_id", columns="program_id", values="camera_t").fillna(0.0)
corr_prot = np.corrcoef(piv_prot.values)
corr_prog = np.corrcoef(piv_prog.values)
iu = np.triu_indices(len(piv_prot), k=1)
r_sim, p_sim = stats.pearsonr(corr_prot[iu], corr_prog[iu])

# Within-category vs between-category similarity in program space
cats = comp_df.set_index("compound_id")["Category"].reindex(piv_prog.index).fillna("Unannotated").values
same_cat = (cats[:, None] == cats[None, :])[iu] & (cats[iu[0]] != "Unannotated")
diff_cat = (cats[:, None] != cats[None, :])[iu] & (cats[iu[0]] != "Unannotated") & (cats[iu[1]] != "Unannotated")
sim_summary = pd.DataFrame([
    {"comparison": "Proteome_vs_Program_pairwise_r", "value": round(float(r_sim), 4), "p_val": float(p_sim)},
    {"comparison": "Mean_within_category_program_r", "value": round(float(corr_prog[iu][same_cat].mean()), 4), "p_val": np.nan},
    {"comparison": "Mean_between_category_program_r", "value": round(float(corr_prog[iu][diff_cat].mean()), 4), "p_val": np.nan},
])
sim_summary.to_csv("reanalysis/A07/drug_similarity.tsv", sep="\t", index=False)

# 5. Candidate report
bini = ev[ev["compound_id"] == "Binimetinib"].iloc[0]
baic = ev[ev["compound_id"] == "Baicalein"].iloc[0]
ndga = ev[ev["compound_id"] == "Nordihydroguaiaretic acid"].iloc[0]
lona = ev[ev["compound_id"] == "Lonafarnib"].iloc[0]

rep_a07 = f"""# A07 Candidate Mechanisms & Pharmacology Report

## Exploratory Candidates (Post-Specificity & Correlation-Aware Audit)
1. **Binimetinib** (MEK1/2 inhibitor; `Kinase_Inhibitors`):
   - **Favorable evidence**: Strongest `PROG_OXPHOS` relative-abundance elevation (`camera_t = {bini['PROG_OXPHOS']:+.2f}`, specificity `resid_z = {bini['PROG_OXPHOS_resid_z']:+.2f}`) coupled with selective suppression of `PROG_INFLAMMATORY` (`camera_t = {bini['PROG_INFLAMMATORY']:+.2f}`, `resid_z = {bini['PROG_INFLAMMATORY_resid_z']:+.2f}`). Distinct from EGFR inhibitors `Erlotinib` and `Mobocertinib`.
   - **Contradictory / Unresolved**: Does not significantly oppose `PRIOR_DEMAGALHAES_UP/DN` under global BH FDR; single MEK1/2 inhibitor in screen (cannot generalize to class); HepG2 hepatoma context; no measured cell viability/counts or Seahorse OCR.
   - **Falsifiable next experiment**: Dose-response in primary human fibroblasts/hepatocytes measuring Seahorse OCR (normalized to cell count & mtDNA) and secreted IL-6/ICAM1 +/- MEK1/2 rescue.
2. **Baicalein** (Flavonoid; `Metabolism_Detox`):
   - **Favorable evidence**: Elevates `PROG_OXPHOS` (`camera_t = {baic['PROG_OXPHOS']:+.2f}`, `resid_z = {baic['PROG_OXPHOS_resid_z']:+.2f}`) while suppressing `PROG_INFLAMMATORY` (`camera_t = {baic['PROG_INFLAMMATORY']:+.2f}`, `resid_z = {baic['PROG_INFLAMMATORY_resid_z']:+.2f}`).
   - **Contradictory / Unresolved**: Flavonoids are promiscuous/PAINS-prone at 10 uM; no viability or primary-cell validation.
   - **Falsifiable next experiment**: Primary-cell senescence/SASP suppression assay with cell-count viability control and redox-inactive analog comparator.

## Literature Comparators
- **Nordihydroguaiaretic acid (NDGA)**: `PROG_OXPHOS = {ndga['PROG_OXPHOS']:+.2f}`, `PROG_MTORC1_SIGNALING = {ndga['PROG_MTORC1_SIGNALING']:+.2f}`.
- **Lonafarnib**: `PROG_MACROAUTOPHAGY = {lona['PROG_MACROAUTOPHAGY']:+.2f}`, `PROG_OXPHOS = {lona['PROG_OXPHOS']:+.2f}`.
- **DFO**: Assay response control (iron chelator / hypoxia mimetic), not a geroprotector control.
- **Proteome vs. Program Similarity**: Pairwise compound correlation in 14-program space correlates with full 7,871-protein space at $r = {r_sim:.3f}$.
"""
with open("reanalysis/A07/candidate_report.md", "w", encoding="utf-8") as f:
    f.write(rep_a07)
print("A07 complete.")
