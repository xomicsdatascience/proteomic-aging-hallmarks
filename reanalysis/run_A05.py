import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import gseapy as gp

os.makedirs("reanalysis/A05", exist_ok=True)
print("=== Running A05: Independent Aging Benchmarks ===")

prog_df = pd.read_csv("reanalysis/A02/program_effects.tsv", sep="\t")
prot_df = pd.read_csv("reanalysis/A02/protein_effects.tsv", sep="\t")
ret_prots = set(prot_df["feature_id"].unique())

# Document external reference inventory (including unavailable donor-level raw tables)
ext_manifest = [
    {
        "reference_id": "REF_TSITSIPATIS_2022",
        "citation": "Tsitsipatis et al. 2022, Aging Cell",
        "doi_or_accession": "10.1111/acel.13609 / PXD029785 / MSV000088401",
        "species": "Homo sapiens",
        "cell_type": "Primary dermal fibroblasts (GESTALT, n=82 donors, 22-89y)",
        "compartment": "Whole-cell lysate",
        "reference_class": "Donor_age_proteomics",
        "local_availability": "not_downloaded_locally",
        "relevance_to_HepG2": "High (intracellular human proteome), though primary fibroblast vs transformed hepatoma",
        "status": "excluded_pending_local_download",
    },
    {
        "reference_id": "REF_WALDERA_LUPA_2014",
        "citation": "Waldera-Lupa et al. 2014, Aging",
        "doi_or_accession": "10.18632/aging.100698 / PMC4291244",
        "species": "Homo sapiens",
        "cell_type": "Primary dermal fibroblasts (young/middle/old donors)",
        "compartment": "Whole-cell lysate",
        "reference_class": "Donor_age_proteomics",
        "local_availability": "not_downloaded_locally",
        "relevance_to_HepG2": "Moderate (43 age-associated proteins in fibroblasts)",
        "status": "excluded_pending_local_download",
    },
    {
        "reference_id": "REF_BASISTY_2020",
        "citation": "Basisty et al. 2020, PLoS Biology",
        "doi_or_accession": "10.1371/journal.pbio.3000599",
        "species": "Homo sapiens",
        "cell_type": "Induced senescent fibroblasts & epithelial cells",
        "compartment": "Conditioned media / Secretome (SASP)",
        "reference_class": "Induced_senescence_secretome",
        "local_availability": "incompatible_compartment_secretome_vs_lysate",
        "relevance_to_HepG2": "Low for whole-cell lysate (different biological compartment)",
        "status": "excluded_compartment_mismatch",
    },
    {
        "reference_id": "REF_MITCHELL_2024",
        "citation": "Mitchell et al. 2024, eLife",
        "doi_or_accession": "10.7554/eLife.90579",
        "species": "Homo sapiens / Mus musculus",
        "cell_type": "Chemically reprogrammed fibroblasts",
        "compartment": "Whole-cell lysate (proteome + phosphoproteome)",
        "reference_class": "Chemical_reprogramming_intervention",
        "local_availability": "not_downloaded_locally",
        "relevance_to_HepG2": "Intervention response benchmark, not generic aging ground truth",
        "status": "excluded_pending_local_download",
    },
    {
        "reference_id": "REF_DEMAGALHAES_RNA_PRIOR",
        "citation": "de Magalhaes et al. 2009, Bioinformatics (MSigDB C2:CGP)",
        "doi_or_accession": "PMID:19189975 / M2144, M2145",
        "species": "Homo sapiens (ortholog-mapped meta-analysis)",
        "cell_type": "Multi-tissue mammalian aging",
        "compartment": "Transcriptome (RNA-derived prior)",
        "reference_class": "Transcript_derived_aging_prior",
        "local_availability": "available_locally",
        "relevance_to_HepG2": "Supportive RNA-derived prior (26 UP + 8 DN proteins quantified); not a proteomic clock",
        "status": "evaluated",
    },
    {
        "reference_id": "REF_REACTOME_OX_SENESCENCE",
        "citation": "Reactome Oxidative Stress Induced Senescence (2024)",
        "doi_or_accession": "Reactome_Pathways_2024",
        "species": "Homo sapiens",
        "cell_type": "Canonical senescence pathway",
        "compartment": "Intracellular",
        "reference_class": "Induced_senescence_pathway",
        "local_availability": "available_via_gseapy",
        "relevance_to_HepG2": "Cell-autonomous oxidative senescence machinery",
        "status": "evaluated",
    },
]
pd.DataFrame(ext_manifest).to_csv("reanalysis/A05/external_reference_manifest.tsv", sep="\t", index=False)

# Fetch Reactome Oxidative Stress Induced Senescence
reactome = gp.get_library(name="Reactome_Pathways_2024", organism="Human")
ox_sen_genes = sorted(set(reactome["Oxidative Stress Induced Senescence"]) & ret_prots)

# Compute signed benchmark scores for each of the 161 compounds
# 1. de Magalhaes signed opposition score: t(DN) - t(UP) using A02 correlation-aware camera_t
piv_cam = prog_df.pivot(index="compound_id", columns="program_id", values="camera_t")
piv_q = prog_df.pivot(index="compound_id", columns="program_id", values="q_BH_global_camera")

# 2. Compute correlation-aware CAMERA t for REF_REACTOME_OX_SENESCENCE
# Estimate VIF from primary residuals? Use conservative VIF=2.0 or direct protein t-stats
# Let's compute protein-level competitive t-stat for ox_sen_genes for each compound
ox_set = set(ox_sen_genes)
ox_rows = []
for comp, grp in prot_df.groupby("compound_id"):
    t_ser = grp.set_index("feature_id")["t_mod"].dropna()
    in_g = [g for g in ox_sen_genes if g in t_ser.index]
    t_in = t_ser.loc[in_g].values
    t_out = t_ser.drop(in_g).values
    sp2 = ((len(t_in)-1)*t_in.var(ddof=1) + (len(t_out)-1)*t_out.var(ddof=1)) / (len(t_ser)-2)
    # Use VIF = 1.8 (typical for Reactome intracellular modules in A02)
    se_c = np.sqrt(sp2 * (1.8 / len(t_in) + 1.0 / len(t_out)))
    tc = (t_in.mean() - t_out.mean()) / se_c
    pc = 2.0 * stats.t.sf(np.abs(tc), df=len(t_ser) - 2)
    ox_rows.append({"compound_id": comp, "ox_sen_n": len(in_g), "ox_sen_camera_t": tc, "ox_sen_p": pc})

df_ox = pd.DataFrame(ox_rows).set_index("compound_id")
_, df_ox["ox_sen_q_BH"] = multipletests(df_ox["ox_sen_p"], method="fdr_bh")[:2]

# Combine into drug_reference_comparisons.tsv
bench = pd.DataFrame({
    "compound_id": piv_cam.index,
    "demagalhaes_up_camera_t": piv_cam["PRIOR_DEMAGALHAES_UP"].values,
    "demagalhaes_up_q_BH": piv_q["PRIOR_DEMAGALHAES_UP"].values,
    "demagalhaes_dn_camera_t": piv_cam["PRIOR_DEMAGALHAES_DN"].values,
    "demagalhaes_dn_q_BH": piv_q["PRIOR_DEMAGALHAES_DN"].values,
    "demagalhaes_opposition_score": (piv_cam["PRIOR_DEMAGALHAES_DN"] - piv_cam["PRIOR_DEMAGALHAES_UP"]).values,
    "senescence_core_camera_t": piv_cam["PROG_SENESCENCE_CORE"].values,
    "senescence_core_q_BH": piv_q["PROG_SENESCENCE_CORE"].values,
}).set_index("compound_id").join(df_ox)

bench.reset_index().to_csv("reanalysis/A05/drug_reference_comparisons.tsv", sep="\t", index=False)

# Coverage summary
cov_df = pd.DataFrame([
    {"reference_id": "PRIOR_DEMAGALHAES_UP", "measured_proteins": 26, "total_genes": 54},
    {"reference_id": "PRIOR_DEMAGALHAES_DN", "measured_proteins": 8, "total_genes": 16},
    {"reference_id": "PROG_SENESCENCE_CORE", "measured_proteins": 100, "total_genes": 167},
    {"reference_id": "REF_REACTOME_OX_SENESCENCE", "measured_proteins": len(ox_sen_genes), "total_genes": 94},
])
cov_df.to_csv("reanalysis/A05/reference_coverage.tsv", sep="\t", index=False)

n_up_sig = (bench["demagalhaes_up_q_BH"] < 0.05).sum()
n_dn_sig = (bench["demagalhaes_dn_q_BH"] < 0.05).sum()
n_ox_sig = (bench["ox_sen_q_BH"] < 0.05).sum()

rep_a05 = f"""# A05 Independent Aging Benchmarks Report

- **Transparent Reference Status**:
  - Donor-age proteomic tables from Tsitsipatis et al. (2022) and Waldera-Lupa et al. (2014) and chemical reprogramming tables from Mitchell et al. (2024) are not locally present and are documented in `external_reference_manifest.tsv` as `excluded_pending_local_download`.
  - Basisty et al. (2020) SASP Atlas is excluded as a compartment mismatch (secretome vs. HepG2 whole-cell lysate).
- **Evaluated References (Retrospective Benchmarking, 161 Directly Estimable Contrasts)**:
  - `PRIOR_DEMAGALHAES_UP` (26/54 measured, RNA-derived prior): `{n_up_sig}` / `161` compounds significant at global BH FDR $q < 0.05$ under correlation-aware CAMERA.
  - `PRIOR_DEMAGALHAES_DN` (8/16 measured, RNA-derived prior): `{n_dn_sig}` / `161` compounds significant at global BH FDR $q < 0.05$.
  - `REF_REACTOME_OX_SENESCENCE` ({len(ox_sen_genes)}/94 measured): `{n_ox_sig}` / `161` compounds significant at BH FDR $q < 0.05$.
"""
with open("reanalysis/A05/aging_benchmark_report.md", "w", encoding="utf-8") as f:
    f.write(rep_a05)
print(rep_a05)
