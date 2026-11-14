# A06 Shared Stress, Growth & Relative Composition Specificity Report

- **Missing Measured Growth Data**: Cell counts, viability, and total protein per cell are not present in `raw (1).h5ad` or `drug_details (1).xlsx`. Cytostasis or toxicity cannot be ruled out; molecular proxies (`NUIS_E2F_PROLIF`, `NUIS_G2M_PROLIF`, `NUIS_APOPTOSIS`) are used for sensitivity only.
- **Correlation with Proliferation Proxy (`NUIS_E2F_PROLIF`)**:
  - `PROG_OXPHOS` vs `NUIS_E2F_PROLIF`: $r = -0.095$ ($p = 2.33e-01$)
  - `PROG_INFLAMMATORY` vs `NUIS_E2F_PROLIF`: $r = -0.272$ ($p = 4.90e-04$)
- **Candidate Specificity Residuals** (`candidate_specificity.tsv`, Z-scored after adjusting for `NUIS_E2F_PROLIF` and `mean_abs_t`):
  - `Binimetinib`: `PROG_OXPHOS_resid_z = +1.88`, `PROG_INFLAMMATORY_resid_z = -4.18`
  - `Baicalein`: `PROG_OXPHOS_resid_z = +1.37`, `PROG_INFLAMMATORY_resid_z = -4.28`
- **Relative-Composition Check**: Sample-median normalization measures relative proteome composition. `composition_sensitivity.tsv` reports `oxphos_mean_log2fc` alongside `cytosolic_ribo_mean_log2fc`.
