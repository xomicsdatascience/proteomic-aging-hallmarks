# Reanalysis Results Summary (A01-A08)

- **Title**: Proteomic profiling maps pharmacological responses across aging associated cellular programs
- **Primary Estimable Universe**: 161 contrasts (160 library drugs + DFO) across 12 plates (`Blue1-6`, `Red1-6`); 12 `RedMix` drugs excluded from primary inference due to 0 within-plate `DMSO` controls.
- **Recomputed Significance (Global BH FDR q < 0.05, 2,254 tests)**:
  - Correlation-aware CAMERA: 218 / 2,254 pairs across 142 / 161 compounds.
  - Pre-ranked GSEA (2,000 permutations): 826 / 2,254 pairs across 161 / 161 compounds.
- **Null Calibration**: 0 / 7,872 proteins and 0 / 14 aging programs significant at BH q < 0.05 in within-plate `DMSO` split.
- **Exploratory Candidates**: `Binimetinib` (`PROG_OXPHOS_resid_z = +1.88`, `PROG_INFLAMMATORY_resid_z = -4.15`) and `Baicalein` (`PROG_OXPHOS_resid_z = +1.37`, `PROG_INFLAMMATORY_resid_z = -4.28`) are nominated for orthogonal functional validation.
