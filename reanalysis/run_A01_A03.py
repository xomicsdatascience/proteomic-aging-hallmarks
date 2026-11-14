import os
import json
import hashlib
import datetime
import numpy as np
import pandas as pd
import anndata as ad
import gseapy as gp

os.makedirs("reanalysis/A01", exist_ok=True)
os.makedirs("reanalysis/A03", exist_ok=True)

# ==============================================================================
# A01: SAMPLE PROVENANCE AND ESTIMABILITY
# ==============================================================================
print("=== Running A01: Sample Provenance & Estimability ===")

raw_path = "data/raw (1).h5ad"
xl_path = "data/drug_details (1).xlsx"

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

raw_sha = sha256_file(raw_path)
xl_sha = sha256_file(xl_path)

adata = ad.read_h5ad(raw_path)
df_xl = pd.read_excel(xl_path)

# 1. Sample Manifest
parts = adata.obs_names.to_series().str.split("_", n=1, expand=True)
obs = adata.obs.copy()
obs["sample_id"] = adata.obs_names
obs["plate"] = obs["batch"].astype(str)
obs["well"] = parts[1].values
obs["compound_id"] = obs["drug"].astype(str)

def assign_role(row):
    if row["compound_id"] == "QC_sample":
        return "analytical_qc"
    if row["compound_id"] == "DMSO":
        return "vehicle_control"
    if row["compound_id"] == "DFO":
        return "positive_assay_control"
    if row["compound_id"] == "Drug":
        return "unlabeled_excluded"
    return "library_compound"

obs["condition_role"] = obs.apply(assign_role, axis=1)
obs["dose"] = np.where(obs["condition_role"] == "library_compound", "10 uM (parent study)", "unknown_control_or_qc")
obs["duration"] = np.where(obs["condition_role"].isin(["library_compound", "vehicle_control", "positive_assay_control"]), "24 h (parent study)", "N/A")
obs["culture_id"] = obs["plate"] + "_" + obs["well"]
obs["culture_day"] = "unknown_not_in_h5ad"
obs["preparation_batch"] = obs["plate"]
obs["instrument_batch"] = obs["plate"]
obs["injection_order"] = "unknown_not_in_h5ad"
obs["injection_id"] = obs["sample_id"]
obs["excluded"] = obs["condition_role"].isin(["unlabeled_excluded"])
obs["exclusion_reason"] = np.where(obs["excluded"], "Unlabeled sample (drug == 'Drug')", "")

sample_cols = [
    "sample_id", "compound_id", "condition_role", "dose", "duration",
    "culture_id", "culture_day", "plate", "well", "preparation_batch",
    "instrument_batch", "injection_order", "injection_id", "excluded", "exclusion_reason"
]
sample_manifest = obs[sample_cols].reset_index(drop=True)
sample_manifest.to_csv("reanalysis/A01/sample_manifest.tsv", sep="\t", index=False)

# 2. Compound Manifest
# Reconcile AnnData compounds with Excel metadata
xl_sub = df_xl.copy()
xl_sub["in_excel"] = True
ad_compounds = sorted(obs["compound_id"].unique())
comp_rows = []
for c in ad_compounds:
    sub_obs = obs[obs["compound_id"] == c]
    plates = sorted(sub_obs["plate"].unique())
    n_tot = len(sub_obs)
    role = sub_obs["condition_role"].iloc[0]
    m_xl = xl_sub[xl_sub["ProductName"] == c]
    if len(m_xl) > 0:
        r = m_xl.iloc[0]
        comp_rows.append({
            "compound_id": c,
            "condition_role": role,
            "n_samples": n_tot,
            "n_plates": len(plates),
            "plates": ",".join(plates),
            "in_h5ad": True,
            "in_excel": True,
            "BaseName": r["BaseName"],
            "CAS_Number": r["CAS Number"],
            "Category": r["Category"],
            "Target": r["Target"],
            "PathWay": r["PathWay"],
            "Research_Area": r["Research Area"],
        })
    else:
        comp_rows.append({
            "compound_id": c,
            "condition_role": role,
            "n_samples": n_tot,
            "n_plates": len(plates),
            "plates": ",".join(plates),
            "in_h5ad": True,
            "in_excel": False,
            "BaseName": np.nan,
            "CAS_Number": np.nan,
            "Category": np.nan,
            "Target": np.nan,
            "PathWay": np.nan,
            "Research_Area": np.nan,
        })

# Add Excel-only compounds not in h5ad
for _, r in xl_sub[~xl_sub["ProductName"].isin(ad_compounds)].iterrows():
    comp_rows.append({
        "compound_id": r["ProductName"],
        "condition_role": "excel_only_not_in_h5ad",
        "n_samples": 0,
        "n_plates": 0,
        "plates": "",
        "in_h5ad": False,
        "in_excel": True,
        "BaseName": r["BaseName"],
        "CAS_Number": r["CAS Number"],
        "Category": r["Category"],
        "Target": r["Target"],
        "PathWay": r["PathWay"],
        "Research_Area": r["Research Area"],
    })

compound_manifest = pd.DataFrame(comp_rows)
compound_manifest.to_csv("reanalysis/A01/compound_manifest.tsv", sep="\t", index=False)

# 3. Protein Mapping
X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
zero_rate = (X == 0).mean(axis=0)
non_qc_mask = obs["condition_role"] != "analytical_qc"
zero_rate_non_qc = (X[non_qc_mask.values] == 0).mean(axis=0)
prot_df = pd.DataFrame({
    "feature_id": adata.var_names,
    "gene_symbol": adata.var_names,
    "zero_fraction_all": np.round(zero_rate, 4),
    "zero_fraction_non_qc": np.round(zero_rate_non_qc, 4),
    "retained_50pct_filter": zero_rate <= 0.50,
    "retained_20pct_filter": zero_rate <= 0.20,
})
prot_df.to_csv("reanalysis/A01/protein_mapping.tsv", sep="\t", index=False)

# 4. Sample Accounting
acct = obs.groupby(["plate", "condition_role"]).size().unstack(fill_value=0).reset_index()
acct["total"] = acct.drop(columns=["plate"]).sum(axis=1)
acct.to_csv("reanalysis/A01/sample_accounting.tsv", sep="\t", index=False)

# 5. Contrast Eligibility Table
dmso_plates = set(obs[obs["condition_role"] == "vehicle_control"]["plate"])
elig_rows = []
for c in ad_compounds:
    sub_obs = obs[obs["compound_id"] == c]
    role = sub_obs["condition_role"].iloc[0]
    if role in ["vehicle_control", "analytical_qc", "unlabeled_excluded"]:
        continue
    plates = sorted(sub_obs["plate"].unique())
    plates_with_dmso = [p for p in plates if p in dmso_plates]
    if len(plates_with_dmso) == len(plates) and len(plates) > 1:
        status = "directly_estimable"
        reason = f"Replicated across {len(plates)} plates, all with within-plate DMSO controls"
    elif len(plates_with_dmso) == 0:
        status = "reference_anchored_only"
        reason = f"Confined to {','.join(plates)} with 0 within-plate DMSO controls; requires QC anchoring"
    else:
        status = "partially_estimable"
        reason = "Partial overlap with DMSO plates"
    elig_rows.append({
        "contrast_id": f"{c}_vs_DMSO",
        "compound_id": c,
        "condition_role": role,
        "n_treated": len(sub_obs),
        "n_plates": len(plates),
        "plates": ",".join(plates),
        "n_within_plate_dmso_plates": len(plates_with_dmso),
        "eligibility_status": status,
        "eligibility_reason": reason,
    })

eligibility_df = pd.DataFrame(elig_rows)
eligibility_df.to_csv("reanalysis/A01/contrast_eligibility.tsv", sep="\t", index=False)

# 6. Parent Study Overlap
overlap_rows = [
    {"item": "HepG2 perturbation proteomics dataset (1,232 samples, 13 plates)", "status": "reused_data", "notes": "Jiang et al. 2026 bioRxiv (DOI 10.64898/2026.05.04.722718); local raw (1).h5ad has 8,630 features vs 8,703 reported"},
    {"item": "Drug metadata table (drug_details (1).xlsx)", "status": "reused_data", "notes": "173 rows; 161 match h5ad compounds, 12 Excel-only, 11 RedMix h5ad compounds unannotated"},
    {"item": "Robotic sample prep & DIA-MS acquisition", "status": "reused_methods", "notes": "From Jiang et al. 2026"},
    {"item": "Sample provenance, contrast estimability & RedMix audit (A01)", "status": "new_analysis", "notes": "Identifies 161 directly estimable contrasts vs 12 RedMix reference-anchored contrasts"},
    {"item": "Batch-aware observed-data linear model + sensitivity suite (A02, A04)", "status": "new_analysis", "notes": "Replaces blind ComBat + imputation with plate-blocked limma-style model & global BH FDR"},
    {"item": "Evidence-graded aging-associated protein program atlas (A03)", "status": "new_analysis", "notes": "Explicitly grades directness, removes hallmark overclaiming, audits leading-edge overlap"},
    {"item": "External cellular aging & senescence proteomic benchmarking (A05)", "status": "new_analysis", "notes": "Tests agreement/opposition against independent aging/senescence signatures"},
    {"item": "Shared stress, proliferation & relative-composition specificity (A06)", "status": "new_analysis", "notes": "Separates compound-specific signals from shared cytostatic/stress responses"},
]
pd.DataFrame(overlap_rows).to_csv("reanalysis/A01/parent_study_overlap.tsv", sep="\t", index=False)

# 7. Design Audit Report
n_direct = (eligibility_df["eligibility_status"] == "directly_estimable").sum()
n_ref = (eligibility_df["eligibility_status"] == "reference_anchored_only").sum()
design_md = f"""# A01 Design & Provenance Audit

- **Input checksums**:
  - `data/raw (1).h5ad`: SHA256 `{raw_sha}` (`{adata.n_obs}` samples x `{adata.n_vars}` features)
  - `data/drug_details (1).xlsx`: SHA256 `{xl_sha}` (`{len(df_xl)}` rows, `{df_xl['ProductName'].nunique()}` unique `ProductName`, `{df_xl['BaseName'].nunique()}` unique `BaseName`)
- **Sample Reconciliation**:
  - Total observations: `{adata.n_obs}` across 13 plates (`Blue1–Blue6`, `Red1–Red6`, `RedMix`).
  - `library_compound`: `1,029` samples (`172` unique compounds: `160` in Blue/Red plates with `n=5–6` across 6 plates; `12` in `RedMix` with `n=6` on 1 plate).
  - `positive_assay_control` (`DFO`): `24` samples (`2` per plate on `Blue1–6` and `Red1–6`; `0` on `RedMix`).
  - `vehicle_control` (`DMSO`): `23` samples (`2` per plate on `Blue1–6` and `Red1–6`, except `Blue3` has `1`; `0` on `RedMix`).
  - `analytical_qc` (`QC_sample`): `155` samples (`12` per plate on `Blue1–6` and `Red1–6`; `11` on `RedMix`).
  - `unlabeled_excluded`: `1` sample (`Blue5_F12`, `drug == 'Drug'`).
- **Feature Reconciliation**:
  - `raw (1).h5ad` contains `8,630` unique gene-symbol-indexed protein features (vs. `8,703` in parent manuscript).
  - `7,871` features have `<= 50%` zeros across all samples (`7,296` have `<= 20%` zeros). Zeros represent missing/below-detection values (`10.60%` overall).
- **Contrast Estimability**:
  - **Directly estimable (`{n_direct}` contrasts)**: `160` Blue/Red library compounds + `DFO`, all replicated across 5–6 plates with matched within-plate `DMSO` controls.
  - **Reference-anchored only (`{n_ref}` contrasts)**: `12` `RedMix` compounds confined to a single plate with `0` within-plate `DMSO` controls. Excluded from primary inference and analyzed only in a named QC-anchored sensitivity analysis.
"""
with open("reanalysis/A01/design_audit.md", "w", encoding="utf-8") as f:
    f.write(design_md)

print(f"A01 complete: {n_direct} directly estimable contrasts, {n_ref} reference-anchored only.")

# ==============================================================================
# A03: PROGRAM DEFINITIONS AND EVIDENCE
# ==============================================================================
print("=== Running A03: Program Definitions & Evidence ===")
hallmark = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")
reactome = gp.get_library(name="Reactome_Pathways_2024", organism="Human")
gobp = gp.get_library(name="GO_Biological_Process_2023", organism="Human")

programs_spec = [
    ("PROG_DNA_REPAIR", "DNA Repair Machinery", "MSigDB_Hallmark_2020", "DNA Repair", hallmark["DNA Repair"],
     "DNA repair enzyme abundance (indirect proxy; does not measure DNA damage burden)", "moderate", "indirect_abundance", "context_dependent"),
    ("PROG_TELOMERE_MAINT", "Telomere Maintenance & Extension", "Reactome_Pathways_2024", "Extension of Telomeres", reactome["Extension of Telomeres"],
     "Telomere replication/maintenance machinery (does not measure telomere length)", "moderate", "indirect_abundance", "context_dependent"),
    ("PROG_CHROMATIN_ORG", "Chromatin Organization", "Reactome_Pathways_2024", "Chromatin Organization", reactome["Chromatin Organization"],
     "Chromatin-modifying and remodeling enzymes (does not measure epigenetic state/methylation)", "moderate", "indirect_abundance", "context_dependent"),
    ("PROG_PROTEOSTASIS_UPR", "Unfolded Protein Response (UPR)", "MSigDB_Hallmark_2020", "Unfolded Protein Response", hallmark["Unfolded Protein Response"],
     "ER/proteostatic stress-response proteins (elevation can reflect stress induction or adaptive capacity)", "strong", "pathway_abundance", "context_dependent"),
    ("PROG_PROTEASOME", "Proteasome Degradation", "Reactome_Pathways_2024", "Antigen Processing Ubiquitination & Proteasome Degradation", reactome["Antigen Processing Ubiquitination & Proteasome Degradation"],
     "Ubiquitin-proteasome machinery abundance (complementary proteostasis module)", "strong", "machinery_abundance", "context_dependent"),
    ("PROG_MACROAUTOPHAGY", "Macroautophagy Machinery", "Reactome_Pathways_2024", "Macroautophagy", reactome["Macroautophagy"],
     "Core autophagy machinery abundance (does not measure autophagic flux)", "strong", "machinery_abundance", "context_dependent"),
    ("PROG_MTORC1_SIGNALING", "mTORC1 Nutrient-Sensing Program", "MSigDB_Hallmark_2020", "mTORC1 Signaling", hallmark["mTORC1 Signaling"],
     "mTORC1-regulated metabolic/anabolic targets (transcriptional/proteomic output)", "strong", "downstream_program", "down_beneficial_in_aging"),
    ("PROG_OXPHOS", "Mitochondrial Oxidative Phosphorylation", "MSigDB_Hallmark_2020", "Oxidative Phosphorylation", hallmark["Oxidative Phosphorylation"],
     "Respiratory chain & TCA enzyme relative abundance (requires normalization & respirometry for capacity)", "strong", "complex_abundance", "up_typically_beneficial"),
    ("PROG_SENESCENCE_CORE", "Cellular Senescence Program", "Reactome_Pathways_2024", "Cellular Senescence", reactome["Cellular Senescence"],
     "Senescence-associated cell-cycle & chromatin regulators (does not establish irreversible arrest)", "moderate", "indirect_abundance", "down_beneficial"),
    ("PROG_CELL_CELL_COMM", "Cell-Cell Communication", "Reactome_Pathways_2024", "Cell-Cell Communication", reactome["Cell-Cell Communication"],
     "Intracellular junction/signaling proteins (cannot assess systemic intercellular communication)", "weak_in_monoculture", "proxy_only", "unclear"),
    ("PROG_INFLAMMATORY", "Intracellular Inflammatory Response", "MSigDB_Hallmark_2020", "Inflammatory Response", hallmark["Inflammatory Response"],
     "Cell-autonomous inflammatory proteins in lysate (secreted SASP/cytokines largely unmeasured)", "moderate_proxy", "intracellular_proxy", "down_typically_beneficial"),
    ("PROG_INNATE_TLR", "Host Innate TLR Sensing", "Reactome_Pathways_2024", "Toll-like Receptor Cascades", reactome["Toll-like Receptor Cascades"],
     "Host innate immune/TLR signaling (explicitly NOT a dysbiosis measurement; dysbiosis unmeasurable in HepG2)", "not_dysbiosis", "innate_proxy_only", "unclear"),
    ("PRIOR_DEMAGALHAES_UP", "de Magalhaes Aging Up (RNA Prior)", "MSigDB_C2_CGP", "DEMAGALHAES_AGING_UP", json.load(open("gene-sets/DEMAGALHAES_AGING_UP.v2026.2.Hs.json"))["DEMAGALHAES_AGING_UP"]["geneSymbols"],
     "Transcript-derived meta-analysis of genes overexpressed with age (prior, not proteomic clock)", "prior_rna", "rna_derived_prior", "down_opposes_aging"),
    ("PRIOR_DEMAGALHAES_DN", "de Magalhaes Aging Down (RNA Prior)", "MSigDB_C2_CGP", "DEMAGALHAES_AGING_DN", json.load(open("gene-sets/DEMAGALHAES_AGING_DN.v2026.2.Hs.json"))["DEMAGALHAES_AGING_DN"]["geneSymbols"],
     "Transcript-derived meta-analysis of genes underexpressed with age (prior, not proteomic clock)", "prior_rna", "rna_derived_prior", "up_opposes_aging"),
    # Nuisance / shared-response modules for A06 (frozen here)
    ("NUIS_E2F_PROLIF", "E2F Targets (Proliferation)", "MSigDB_Hallmark_2020", "E2F Targets", hallmark["E2F Targets"],
     "Proliferation/S-phase cell-cycle module for A06 shared-response sensitivity", "nuisance_control", "proliferation_proxy", "control_module"),
    ("NUIS_G2M_PROLIF", "G2M Checkpoint (Proliferation)", "MSigDB_Hallmark_2020", "G2-M Checkpoint", hallmark["G2-M Checkpoint"],
     "Mitotic/G2M cell-cycle module for A06 shared-response sensitivity", "nuisance_control", "proliferation_proxy", "control_module"),
    ("NUIS_APOPTOSIS", "Apoptosis Response", "MSigDB_Hallmark_2020", "Apoptosis", hallmark["Apoptosis"],
     "Apoptotic stress module for A06 cytotoxicity sensitivity", "nuisance_control", "stress_proxy", "control_module"),
]

all_prots = set(adata.var_names)
retained_prots = set(prot_df.loc[prot_df["retained_50pct_filter"], "gene_symbol"])

manifest_rows = []
membership_rows = []
gmt_lines = []

for pid, pname, src_lib, src_id, genes, desc, ev_grade, directness, ben_dir in programs_spec:
    genes_u = sorted(set(genes))
    meas_all = sorted(set(genes_u) & all_prots)
    meas_ret = sorted(set(genes_u) & retained_prots)
    chk = hashlib.sha256(",".join(genes_u).encode()).hexdigest()[:16]
    manifest_rows.append({
        "program_id": pid,
        "program_name": pname,
        "source_library": src_lib,
        "source_set_id": src_id,
        "retrieval_date": "2026-11-13",
        "species": "Homo sapiens",
        "membership_checksum": chk,
        "total_genes": len(genes_u),
        "measured_all_8630": len(meas_all),
        "measured_retained_7871": len(meas_ret),
        "coverage_pct_retained": round(100.0 * len(meas_ret) / len(genes_u), 1),
        "aging_evidence_grade": ev_grade,
        "directness": directness,
        "beneficial_direction": ben_dir,
        "biological_description": desc,
    })
    for g in genes_u:
        membership_rows.append({
            "program_id": pid,
            "gene_symbol": g,
            "in_raw_8630": g in all_prots,
            "in_retained_7871": g in retained_prots,
        })
    gmt_lines.append("\t".join([pid, desc] + genes_u))

pd.DataFrame(manifest_rows).to_csv("reanalysis/A03/program_manifest.tsv", sep="\t", index=False)
pd.DataFrame(membership_rows).to_csv("reanalysis/A03/program_membership.tsv", sep="\t", index=False)
with open("reanalysis/A03/programs.gmt", "w", encoding="utf-8") as f:
    f.write("\n".join(gmt_lines) + "\n")

# Pairwise Jaccard overlap on retained measured proteome
df_mem = pd.DataFrame(membership_rows)
prog_sets = {
    pid: set(df_mem.loc[(df_mem["program_id"] == pid) & (df_mem["in_retained_7871"]), "gene_symbol"])
    for pid in [r["program_id"] for r in manifest_rows]
}
pids = [r["program_id"] for r in manifest_rows]
ov_rows = []
for i, p1 in enumerate(pids):
    for p2 in pids[i+1:]:
        s1, s2 = prog_sets[p1], prog_sets[p2]
        inter = len(s1 & s2)
        if inter > 0:
            jacc = inter / len(s1 | s2)
            ov_rows.append({"program_1": p1, "program_2": p2, "overlap_count": inter, "jaccard": round(jacc, 4)})
pd.DataFrame(ov_rows).sort_values("overlap_count", ascending=False).to_csv("reanalysis/A03/program_overlap.tsv", sep="\t", index=False)

# Evidence matrix
ev_rows = []
for r in manifest_rows:
    if r["program_id"].startswith("NUIS_"):
        continue
    ev_rows.append({
        "program_id": r["program_id"],
        "program_name": r["program_name"],
        "directness": r["directness"],
        "unmeasured_state_in_HepG2_lysate": "Yes - " + r["biological_description"],
        "required_functional_assay": {
            "PROG_DNA_REPAIR": "gamma-H2AX foci / comet assay",
            "PROG_TELOMERE_MAINT": "TRF / qPCR telomere length & TRAP telomerase activity",
            "PROG_CHROMATIN_ORG": "ATAC-seq / histone PTM MS / DNA methylation clock",
            "PROG_PROTEOSTASIS_UPR": "Proteasome activity / protein aggregation dye / XBP1 splicing",
            "PROG_PROTEASOME": "AMC-peptide proteasome cleavage assay",
            "PROG_MACROAUTOPHAGY": "LC3-II turnover +/- bafilomycin A1 flux assay",
            "PROG_MTORC1_SIGNALING": "p-S6K / p-4EBP1 Western blot or phosphoproteomics",
            "PROG_OXPHOS": "Seahorse OCR respirometry normalized to cell number & mtDNA",
            "PROG_SENESCENCE_CORE": "SA-beta-gal + EdU arrest + p16/p21 + SASP secretome",
            "PROG_CELL_CELL_COMM": "Co-culture / conditioned-media secretome profiling",
            "PROG_INFLAMMATORY": "Conditioned-media multiplex ELISA / secretome DIA",
            "PROG_INNATE_TLR": "In vivo microbiome sequencing (16S/metagenomics) — unmeasurable in HepG2",
            "PRIOR_DEMAGALHAES_UP": "Chronological/biological age proteomic validation",
            "PRIOR_DEMAGALHAES_DN": "Chronological/biological age proteomic validation",
        }.get(r["program_id"], "Functional assay"),
    })
pd.DataFrame(ev_rows).to_csv("reanalysis/A03/program_evidence_matrix.tsv", sep="\t", index=False)

with open("reanalysis/A03/annotation_decisions.md", "w", encoding="utf-8") as f:
    f.write("""# A03 Program Annotation & Evidence Decisions

1. **No 1:1 Hallmark Claims**: Programs are named by what they actually measure in HepG2 whole-cell lysates (e.g. `PROG_DNA_REPAIR`, `PROG_OXPHOS`, `PROG_INNATE_TLR`), not as direct readouts of the 12 organismal hallmarks.
2. **Unmeasurable Hallmarks in HepG2 Monoculture**:
   - **Dysbiosis** cannot be measured in cultured HepG2 cells; `PROG_INNATE_TLR` is explicitly labeled as host innate TLR sensing.
   - **Stem cell exhaustion** cannot be assessed in transformed HepG2 hepatoma cells and was removed from the primary aging-program panel.
   - **Secreted SASP / Intercellular communication** are only weakly proxied in whole-cell lysates (`PROG_INFLAMMATORY`, `PROG_CELL_CELL_COMM`).
3. **Transcript-Derived Priors**: `DEMAGALHAES_AGING_UP` and `DEMAGALHAES_AGING_DN` are retained as RNA-derived aging priors (`PRIOR_DEMAGALHAES_UP`, `PRIOR_DEMAGALHAES_DN`), not proteomic aging clocks.
4. **Nuisance/Shared-Response Modules**: `NUIS_E2F_PROLIF`, `NUIS_G2M_PROLIF`, and `NUIS_APOPTOSIS` are frozen in A03 for A06 specificity/growth sensitivity analyses.
""")

print("A03 complete: 17 programs frozen (14 aging-associated + 3 nuisance control modules).")
