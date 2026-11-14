# A07 Candidate Mechanisms & Pharmacology Report

## Exploratory Candidates (Post-Specificity & Correlation-Aware Audit)
1. **Binimetinib** (MEK1/2 inhibitor; `Kinase_Inhibitors`):
   - **Favorable evidence**: Strongest `PROG_OXPHOS` relative-abundance elevation (`camera_t = +2.50`, specificity `resid_z = +1.88`) coupled with selective suppression of `PROG_INFLAMMATORY` (`camera_t = -2.75`, `resid_z = -4.15`). Distinct from EGFR inhibitors `Erlotinib` and `Mobocertinib`.
   - **Contradictory / Unresolved**: Does not significantly oppose `PRIOR_DEMAGALHAES_UP/DN` under global BH FDR; single MEK1/2 inhibitor in screen (cannot generalize to class); HepG2 hepatoma context; no measured cell viability/counts or Seahorse OCR.
   - **Falsifiable next experiment**: Dose-response in primary human fibroblasts/hepatocytes measuring Seahorse OCR (normalized to cell count & mtDNA) and secreted IL-6/ICAM1 +/- MEK1/2 rescue.
2. **Baicalein** (Flavonoid; `Metabolism_Detox`):
   - **Favorable evidence**: Elevates `PROG_OXPHOS` (`camera_t = +1.93`, `resid_z = +1.37`) while suppressing `PROG_INFLAMMATORY` (`camera_t = -2.96`, `resid_z = -4.28`).
   - **Contradictory / Unresolved**: Flavonoids are promiscuous/PAINS-prone at 10 uM; no viability or primary-cell validation.
   - **Falsifiable next experiment**: Primary-cell senescence/SASP suppression assay with cell-count viability control and redox-inactive analog comparator.

## Literature Comparators
- **Nordihydroguaiaretic acid (NDGA)**: `PROG_OXPHOS = +0.92`, `PROG_MTORC1_SIGNALING = +1.62`.
- **Lonafarnib**: `PROG_MACROAUTOPHAGY = +0.23`, `PROG_OXPHOS = +0.88`.
- **DFO**: Assay response control (iron chelator / hypoxia mimetic), not a geroprotector control.
- **Proteome vs. Program Similarity**: Pairwise compound correlation in 14-program space correlates with full 7,871-protein space at $r = 0.847$.
