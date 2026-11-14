# Reanalysis Results Summary (A01-A08)

- **Title**: Proteomic profiling maps pharmacological responses across aging associated cellular programs
- **Primary Estimable Universe**: 161 contrasts (160 library compounds + `DFO`) across 12 balanced assay plates (`Blue1`-`Blue6`, `Red1`-`Red6`, `7,871` proteins retained at `<=50%` missingness); 12 `RedMix` compounds excluded from primary inference due to `0` within-plate `DMSO` controls.
- **Recomputed Pathway Significance (Global BH FDR across 2,254 tests)**:
  - Standard regularized correlation-aware `CAMERA` (`rho = 0.01`, Wu & Smyth 2012): **504 / 2,254** pairs (`159 / 161` contrasts) at `BH q < 0.05` (`736` at `q < 0.10`; `1,071` at `q < 0.20`).
  - Pre-Ranked GSEA (`2,000` permutations): **826 / 2,254** pairs (`161 / 161` contrasts) at `BH q < 0.05`.
  - Unregularized empirical-VIF `CAMERA` (`VIF` up to `24.63` on `PROG_OXPHOS`): **218 / 2,254** pairs at `BH q < 0.05`.
- **14-Pathway Composite Geroprotective Score (`CGS`) Top Nominations (`reanalysis/A06/candidate_specificity.tsv`)**:
  1. **Binimetinib** (MEK1/2 inhibitor): `CGS_unadj = +1.91`, `CGS_adj = +3.29` (#1 in screen), **`CGS_combined = +3.68` (`BH q = 0.0028`)**; `PROG_OXPHOS` `CAMERA t = +7.39` (`q = 6.1e-11`), `PROG_INFLAMMATORY` `CAMERA t = -4.02` (`q = 8.2e-4`).
  2. **Cyclothiazide** (AMPA modulator / benzothiadiazine): `CGS_unadj = +2.58` (`q = 0.0273`, #1 unadjusted, 9/14 programs beneficial vs. `DMSO`), `CGS_adj = +1.53`, **`CGS_combined = +2.91` (`BH q = 0.0185`)**.
  3. **Zanubrutinib** (covalent BTK inhibitor): `CGS_unadj = +2.02` (8/14 beneficial), `CGS_adj = +1.38` (10/14 beneficial), **`CGS_combined = +2.40` (`BH q = 0.0593`)** (concordant with second BTK inhibitor `Orelabrutinib`, 9/14 adjusted programs beneficial).
  4. **Baicalein** (flavonoid): `CGS_unadj = +1.11`, `CGS_adj = +2.03`, **`CGS_combined = +2.22` (`BH q = 0.0836`)**.
  5. **Almonertinib (hydrochloride)** (3rd-gen WT-sparing EGFR inhibitor): `CGS_unadj = +0.85`, `CGS_adj = +2.26` (9/14 beneficial), **`CGS_combined = +2.20` (`BH q = 0.0870`)** (concordant with `Osimertinib` `+0.85`, opposing 1st-gen WT-EGFR inhibitors `Icotinib` `-3.93`, `Lapatinib` `-2.65`, `Erlotinib` `-1.75`).
  6. **Ferulic acid** (phenolic antioxidant): `CGS_unadj = +2.06` (10/14 beneficial), `CGS_adj = +0.95`, **`CGS_combined = +2.13` (`BH q = 0.0981`)**.
  7. **alpha-Lipoic Acid** (mitochondrial cofactor): `CGS_unadj = +1.31` (8/14 beneficial), `CGS_adj = +1.53` (10/14 beneficial), **`CGS_combined = +2.01` (`BH q = 0.1207`)** (concordant with second library entry `Lipoic acid`).
- **Pro-Aging / Stress Negative Controls Verified at Opposite Tail (`BH q <= 0.001`)**:
  - `DFO` (`CGS_combined = -6.15`, `q = 1.3e-7`, #161/161), `Amiodarone (hydrochloride)` (`-5.53`, `q = 2.5e-6`, #160/161), and anticholinergic H1-antihistamines (`Desloratadine` `-4.94`, `Benztropine` `-4.76`, `Procyclidine` `-4.74`, `Chloropyramine` `-4.67`, `Antazoline` `-4.52`).
- **Null Calibration Sanity Check (Methods/Supplement)**: `0 / 7,871` proteins and `0 / 14` aging programs significant at `BH q < 0.05` in within-plate `DMSO_A` vs `DMSO_B` split (`CGS_unadj = +0.016`, `p = 0.987`).
