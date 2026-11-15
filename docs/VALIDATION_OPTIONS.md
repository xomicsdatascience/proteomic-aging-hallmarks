# Follow-up options for liver aging validation

## Recommendation for the present manuscript

The current paper should conclude that deep proteomics resolves drug-induced changes in selected aging-associated protein programs in HepG2 cells. It does not measure all organismal hallmarks, establish a cellular age, or demonstrate rejuvenation. Remove the Medicare proposal and keep clinical organ-age work separate until a feasible cohort and a validated endpoint exist.

The first validation should be a repositioned, randomized replication of the selected treatments with distributed vehicle controls: the current binimetinib and baicalein conditions occupy adjacent fixed Red-plate wells H3 and H2. Their related response cannot be separated from a position effect in this export. The most direct next computational analysis is a comparison with an independent human liver aging proteome. The most attractive clinical follow-up is paired routine CT with a Cedars-Sinai imaging collaborator. These answer different questions and should not be combined into an unvalidated composite score.

## What imaging studies actually exist

| Study | Relevant resource | What it does not establish |
|---|---|---|
| Le Goallec and colleagues, Nature Communications 2022, [DOI 10.1038/s41467-022-29525-9](https://doi.org/10.1038/s41467-022-29525-9) | Liver-specific and abdominal age estimates from UK Biobank liver/pancreas MRI. Code is linked by the authors. | A trained MRI model cannot simply be applied to routine oncology CT, and age prediction accuracy does not validate treatment-induced rejuvenation. |
| Ecker and colleagues, npj Aging 2026, [DOI 10.1038/s41514-026-00377-7](https://doi.org/10.1038/s41514-026-00377-7) | Organ-specific imaging age models, including liver MRI, with disease/mortality associations. Public [training and preprocessing code](https://github.com/lab-midas/biological_age) is linked. | Pretrained liver-weight availability needs confirmation. The study does not demonstrate that binimetinib makes livers younger. UK Biobank access and local-model calibration would still be needed. |
| Deng and colleagues, medRxiv 2026, [DOI 10.64898/2026.05.19.26353206](https://doi.org/10.64898/2026.05.19.26353206), [PubMed record](https://pubmed.ncbi.nlm.nih.gov/42245027/) | A Cedars-Sinai-led CT study using 68,682 scans from 32,882 people, including liver-specific radiomic age estimates. The author list includes Debiao Li, Nicolas Musi, Sara Espinoza, Damini Dey, and Chen Yuan. This is the closest local methodological lead. | This is a preprint and a risk-prediction study, not a drug-response validation study. The local clinical data are not an openly downloadable treatment cohort. Model availability, permission for reuse, and repeat-scan performance need confirmation. |

The literature search did not identify a published comparison showing younger liver imaging age after binimetinib versus another treatment. The opportunity is to design that analysis, not cite an existing demonstration.

## A concrete CT pilot

The initial question should be: **Do patients receiving a binimetinib-containing regimen have a different change in liver imaging age gap than clinically comparable patients receiving another regimen?** Binimetinib is targeted therapy; regimen identity matters more than the generic label chemotherapy.

Start with a feasibility count in the local imaging/EHR resource: patients with the relevant cancer genotype and treatment line, a pretreatment abdominal CT, a comparable follow-up CT, documented regimen dates, and usable non-tumor liver parenchyma. Do not infer that the published CT cohort already contains enough binimetinib exposures.

Use paired scans and a fixed follow-up window rather than comparing cross-sectional post-treatment ages. The endpoint is change in an age-bias-corrected liver predicted-age gap, with baseline age gap and elapsed chronological time accounted for. Compare changes between groups, and report intervals and individual trajectories. Choose follow-up windows from clinical scanning practice and test-retest performance before examining drug associations.

In BRAF-mutant melanoma, a contemporary BRAF/MEK combination such as dabrafenib plus trametinib is a more credible clinical comparator than an unmatched historical monotherapy group. The [FDA binimetinib label](https://www.accessdata.fda.gov/drugsatfda_docs/label/2024/210498s010lbl.pdf) and [FDA trametinib label](https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/204114Orig1s036%2C217513Orig1s010correctedlbl.pdf) establish the relevant combination indications. However, both regimens inhibit MEK, so this tests regimen-specific differences rather than MEK inhibition versus its absence. Genotype, line, prior immunotherapy, dose, concurrent therapy, and cancer response need matching or adjustment. A randomized trial imaging archive would offer a stronger design if accessible; access is not established here.

Before using the endpoint, establish repeatability in stable patients and compatibility of scanner, reconstruction, contrast phase, and contrast timing. An MRI model needs separate CT validation. Exclude or mask hepatic tumors/metastases, resection, and focal treatment effects. Evaluate sensitivity to steatosis, liver injury, edema, weight loss/cachexia, and tumor response. Otherwise a model may report a younger-looking texture because treatment changed pathology or acquisition rather than aging.

A successful observational result would support a treatment-associated imaging change. It would not, by itself, establish restored liver function, reversed cellular aging, or a survival benefit. The imaging and proteomic endpoints are biologically complementary, but they do not directly validate one another without an explicit connection.

## A nearer-term alternative using existing tissue proteomics

Ding and colleagues profiled 516 samples from 13 human tissues across five decades and developed tissue-specific proteomic aging signatures: [Cell 2025, DOI 10.1016/j.cell.2025.06.047](https://doi.org/10.1016/j.cell.2025.06.047). This offers a liver-relevant reference to investigate, conditional on obtaining the protein-by-sample matrix and usable age/sex metadata.

A defined analysis would:

1. Estimate a signed liver protein-age association in the external data, accounting for available sex, batch, and tissue-composition variables. Lock that signature before testing drug responses.
2. Map the external proteins to the HepG2 export and report the shared measured universe. Keep the donor unit independent and audit normal tissue versus disease or transplant-related selection.
3. Test whether each drug's observed protein changes oppose that external age-associated direction. Use plate-level replication and a declared testing family across all eligible drugs, rather than choosing favorable pathways afterward.
4. Repeat after excluding proliferation-related proteins and after removing genes used in the primary program panel. Distinguish a whole-proteome association from a mitochondrial-only association.
5. Present concordant and contradictory drugs, intervals, and coverage. Do not convert an acute HepG2 response into years of liver age or call transcript-based priors independent proteomic validation.

This is a proposed analysis, not a result included in the revised manuscript. It is a closer biological bridge to the current measurements than a large clinical imaging project. A focused follow-up in primary hepatocytes or aged-donor cells, with viability and cell-normalized respiration, could then test the most reproducible response combination experimentally.
