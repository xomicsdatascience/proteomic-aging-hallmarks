# Proteomic Profiling of the 12 Hallmarks of Aging Across 173 Pharmacological Perturbations

This repository contains the data, analysis scripts, curated gene sets, figures, and manuscript report demonstrating that high-throughput quantitative proteomics simultaneously captures all **12 Hallmarks of Aging** (López-Otín et al., *Cell* 2023) plus two aggregate aging signatures (`DEMAGALHAES_AGING_UP` and `DEMAGALHAES_AGING_DN`).

## Repository Structure
- `data/`: Raw proteomic dataset (`raw (1).h5ad`, `1,232 × 8,630`) and drug metadata (`drug_details (1).xlsx`).
- `gene-sets/`: Aggregate aging gene sets (`DEMAGALHAES_AGING_UP`, `DEMAGALHAES_AGING_DN`).
- `1-s2.0-S0092867422013770-main.pdf`: López-Otín et al. (2023) *Cell* review defining the 12 Hallmarks of Aging.
- `run_aging_hallmarks_analysis.py`: Full preprocessing, ComBat batch correction, Welch's *t*-test DE (BH FDR), and Pre-Ranked GSEA pipeline.
- `regenerate_all_figures.py`: Generates all publication-formatted main and supplementary figures.
- `md_to_docx.py`: Compiles the Markdown report and figures into a formatted Word document.
- `results/`: Output tables (`nes_matrix_all_drugs.csv`, `fdr_q_matrix_all_drugs.csv`, `significant_drug_pathway_hits_q05.csv`), figures (`fig1`–`fig5`, `figS1`), and reports (`Aging_Hallmarks_Proteomics_Report.md`, `Aging_Hallmarks_Proteomics_Report.docx`).
