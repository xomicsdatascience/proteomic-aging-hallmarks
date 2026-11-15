# Proteomic responses across aging associated cellular programs

This repository reanalyzes the published HepG2 perturbation proteome from [Jiang et al.](https://doi.org/10.64898/2026.05.04.722718). It compares relative protein responses across **14 selected aging-associated programs** in **160 library entries plus the DFO assay-response control**. It does not directly measure all organismal hallmarks, biological age, or geroprotection.

Start with the [current manuscript](reanalysis/A08/Aging_Proteomics_Manuscript_Revised.md), [Word manuscript](reanalysis/A08/Aging_Proteomics_Manuscript_Revised.docx), [supplementary figures S1-S4](reanalysis/A08/Aging_Proteomics_Supplementary_Figures.docx), and [four-axis figure S5](reanalysis/A08/Aging_Proteomics_Supplementary_Figure_S5.docx). [All pathway results](reanalysis/A08/source_tables/Supplementary_Table_S1_all_compounds.tsv) and a [readable paginated full atlas](reanalysis/A08/figures/figS2_atlas_paginated.pdf) include the complete primary universe.

The [four-axis vignette](docs/FOUR_AXIS_VIGNETTE.md) prioritizes the requested higher OXPHOS, autophagy and proteasome with lower senescence-associated abundance profile. [The five shortlisted profiles](reanalysis/A08/source_tables/Supplementary_Table_S2_four_axis_candidates.tsv), [complete exploratory rankings](reanalysis/A08/source_tables/Supplementary_Table_S3_exploratory_indices.tsv), and [exact score definitions](reanalysis/A06/prioritization_definitions.tsv) expose uncertainties, target hypotheses and fixed positions. Full-screen primary component estimates remain in the existing [Supplementary Table S1](reanalysis/A08/source_tables/Supplementary_Table_S1_all_compounds.tsv). Redundant detailed component exports are reproducible locally under the ignored `reanalysis/_cache/` directory.

## Analysis scope

- Exclude the unlabeled well **before** filtering proteins. The <=50% zero-value filter retains 7,871 features from the 8,630-feature export.
- Fit observed log2 abundances with compound and plate effects using 1,004 assay observations. Do not impute missing abundances. Protein tests require at least three observed values in each contrasted group.
- Use matched plates as the replicate unit for primary pathway inference: program mean log2 fold change minus the measured background mean on each plate, followed by a one-sample t test across plates. Correct the 161 x 14 primary tests together with BH. Proteins within a plate are correlated measurements, not independent biological replicates.
- Report fixed-correlation CAMERA-PR-style and pairwise empirical-correlation analyses as secondary sensitivity checks. GSEA NES is descriptive. Figure 3 marks primary q < 0.05 effects with purple diamonds, with green for binimetinib and orange for baicalein retaining the same significance shapes.
- Restore the 14-program combined response index as exploratory prioritization using its earlier target-sign hypotheses and an equal-weight mean of SD-scaled primary effect estimates. The raw index and its E2F/G2M residual remain separate. There are no score p values, q values, significance boundaries, null covariance estimates from gene-set overlap, or inferred universally favorable directions.
- Require all four primary point-estimate signs for the requested-profile shortlist before ranking its four-axis mean. This sign filter is not a joint significance criterion: 68 library entries pass the signs, but none has all four requested-direction effects at global q < 0.05. The leading three entries share Blue column 12, so position confounding is directly relevant.
- Report the non-null vehicle-position diagnostics, fixed-coordinate context, reference sensitivity, missingness sensitivity, and uncertainty. Binimetinib and baicalein occupy adjacent Red H3/H2 wells, and compound-position effects are not separable in the deposited design. They do not establish strict full-screen FDR calibration.
- Keep clinical organ-age and functional validation separate from the results. [Follow-up options](docs/VALIDATION_OPTIONS.md) cover human liver proteomics and a possible Cedars-Sinai CT collaboration. No Medicare analysis or imaging result is claimed.

## Reproduce

Tested with Python 3.12. Input files and gene sets are committed; no gene-set download is needed during the analysis.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reanalysis/run_all.py
python -m unittest discover -s tests -v
```

On Windows, activate the virtual environment with its Scripts activation command. The workflow runs A01/A03, A02, A04-A07, exploratory prioritization, figures, generated manuscript prose, Word exports, and output validation in dependency order, stopping on errors. Variance fitting uses one BLAS thread for reproducibility. GSEA runs 2,000 permutations with seed 42; allow several minutes.

To regenerate only figures, prose, and Word exports from existing inference tables:

```bash
python reanalysis/run_all.py --figures-only
```

Figure generation verifies hashes of all A01-A07 files and cannot rewrite the statistical evidence. Word layout should be visually reviewed after export; the included build checks do not replace that review.

To regenerate the derived exploratory indices after deliberately updating their definitions, run `python reanalysis/run_prioritization.py` before `--figures-only`. The figure-only workflow consumes the saved A06 indices and leaves all A01-A07 files unchanged.

## Outputs and provenance

| Directory | Contents |
|---|---|
| `data/` | Deposited intensity matrix and drug metadata |
| `reanalysis/A01/` | Sample/compound mapping, eligibility, explicit filtering order, input hashes |
| `reanalysis/A03/` | Frozen GMT, annotation provenance, membership, coverage and limits |
| `reanalysis/A02/` | Protein models, primary plate-level pathway comparisons, secondary correlation sensitivity, direct kinase contrasts |
| `reanalysis/A04/` | Technical QC, vehicle-position diagnostics, all reference comparisons, missingness sensitivity and plate-bootstrap summaries |
| `reanalysis/A05/` | RNA-prior context and an explicit inventory of proposed, unevaluated external references |
| `reanalysis/A06/` | Growth/composition sensitivity, fixed positions, exploratory indices, exact definitions and rankings |
| `reanalysis/A07/` | Complete example protein modules, observation counts, intervals and interpretation |
| `reanalysis/A08/` | Four main figures, supplementary figures, source tables, generated prose, claim ledger, synchronized Word exports and build provenance |

The complete protein table is stored in [indexed gzip partitions](reanalysis/A02/protein_effects/manifest.tsv), avoiding a single oversized generated file. `analysis_common.read_table('reanalysis/A02/protein_effects.tsv.gz')` reads the local full build output when available or reconstructs it from the committed partitions. Never treat an unavailable statistic as zero.

SVG figures retain full color. Compact PNG and Word copies use indexed color at their native pixel dimensions; no image is downsampled. Word image extents and all numerical source values are unchanged by this export encoding.

The [audit resolution note](docs/AUDIT_RESOLUTION.md) records substantive changes, including a correction to the earlier audit's null-protein calculation. Remaining limitations include the transformed cell model, one concentration/time point, missing viability/cell-count metadata, incomplete culture/acquisition provenance, and the absence of independent functional or organismal aging validation.

`archive/pre-audit/` contains superseded exploratory outputs retained for provenance. They are not the current analysis. The historical claim that a single HepG2 lysate screen directly captures all twelve hallmarks is not supported. The root Word report, A08 Final Word file, and results Word report are compatibility copies of the same canonical manuscript.
