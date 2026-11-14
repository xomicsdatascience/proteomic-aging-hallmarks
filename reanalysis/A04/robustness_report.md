# A04 Reproducibility & Sensitivity Report

- **Sensitivity (Primary Observed Plate-Blocked vs. ComBat+Imputation)**:
  - Correlation between primary `camera_t` and ComBat+imputation `NES`: $r = 0.816$ (`84.9%` direction agreement).
  - Correlation between primary `nes_gsea` and ComBat+imputation `NES`: $r = 0.836$.
- **Null Calibration (`DMSO_A` vs `DMSO_B` across 11 plates)**:
  - `0` / `7872` proteins significant at BH FDR $q < 0.05$.
  - `0` / `14` aging programs significant at BH FDR $q < 0.05$ (`0` at nominal $p < 0.05$).
- **Cross-Plate Replication (Leave-One-Plate-Out Pearson $r$ of $\log_2$FC)**:
  - Computed for `DFO` and key candidates (`replication_metrics.tsv`).
