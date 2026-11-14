# A05 Independent Aging Benchmarks Report

- **Transparent Reference Status**:
  - Donor-age proteomic tables from Tsitsipatis et al. (2022) and Waldera-Lupa et al. (2014) and chemical reprogramming tables from Mitchell et al. (2024) are not locally present and are documented in `external_reference_manifest.tsv` as `excluded_pending_local_download`.
  - Basisty et al. (2020) SASP Atlas is excluded as a compartment mismatch (secretome vs. HepG2 whole-cell lysate).
- **Evaluated References (Retrospective Benchmarking, 161 Directly Estimable Contrasts)**:
  - `PRIOR_DEMAGALHAES_UP` (26/54 measured, RNA-derived prior): `1` / `161` compounds significant at global BH FDR $q < 0.05$ under correlation-aware CAMERA.
  - `PRIOR_DEMAGALHAES_DN` (8/16 measured, RNA-derived prior): `5` / `161` compounds significant at global BH FDR $q < 0.05$.
  - `REF_REACTOME_OX_SENESCENCE` (47/94 measured): `0` / `161` compounds significant at BH FDR $q < 0.05$.
