# A02 Inference Report

- **Primary Family**: `2254` hypotheses (`161` directly estimable contrasts x `14` primary aging programs). `12` `RedMix` contrasts were excluded from primary inference due to `0` within-plate `DMSO` controls.
- **Model**: Plate-blocked observed-data OLS (`~ 1 + compound + plate`) with empirical Bayes variance moderation (`d0 = 2.89`, `s0^2 = 0.0204`). Normalization is sample-median $\log_2$, measuring **relative proteome composition** rather than absolute per-cell protein mass.
- **Recomputed Significance Counts (Global BH FDR $q < 0.05$)**:
  - **Correlation-aware competitive test (CAMERA)**: `218` / `2254` pairs across `142` / `161` compounds.
  - **Pre-Ranked GSEA (`q_BH_global_gsea < 0.05`)**: `826` / `2254` pairs across `161` / `161` compounds.
  - *(Compare to preliminary unblocked ComBat + native-q report of 165/173 drugs and 999 pairs.)*
