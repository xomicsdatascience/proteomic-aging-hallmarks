# Analyses and figures for the aging focused proteomics manuscript

Prepared for Jesse G. Meyer on 13 November 2026.

## Scientific objective

Use the existing HepG2 drug perturbation proteome to determine which pharmacological responses affect aging-associated cellular programs, which response combinations distinguish compounds, and whether those combinations agree with independent aging biology after accounting for shared drug responses.

Working title: **Proteomic profiling maps pharmacological responses across aging associated cellular programs**.

The underlying experimental platform and screen are already reported in Jiang et al., *Robotic perturbation proteomics and AI agents enable scalable drug mechanism discovery*, bioRxiv, DOI https://doi.org/10.64898/2026.05.04.722718. The new paper must identify its new analyses and scientific results explicitly. Do not claim the dataset, robotic throughput, generic GSEA, or use of proteomics in aging as new contributions. Do not claim to measure all twelve hallmarks, to predict biological age, to demonstrate rejuvenation, or to replace functional assays.

The reframed DOCX is a working manuscript, not a completed reanalysis. Its numerical examples and six retained source graphics come from the supplied report. Do not reuse their significance calls without verification. The planned final figures below differ from the source figure numbers.

## Sources and required inputs

Author-deposited dataset: https://huggingface.co/datasets/jgmeyer/mississippi

Parent study code: https://github.com/xomicsdatascience/AI-robotics-proteomics

Source working files named in the report: `raw (1).h5ad` and `drug_details (1).xlsx`. The public dataset describes `raw.h5ad` and `drug_details.xlsx`. Establish file identity with checksums and metadata; do not assume that similarly named files are identical.

Required local inputs are the least-processed protein intensity matrix, protein accession mapping, complete sample and plate manifests, QC composition and injection records, biological culture identifiers, compound identities and concentrations, exposure time, exclusion history, and available cell count or viability measurements. Retain accessions and provenance for all external reference data. If a required input is unavailable, document the missing item and its consequence; do not fabricate metadata or substitute a molecular proxy without labeling it.

The parent study reports HepG2 cells, 172 compounds, six biological replicates, 1,232 proteomes, and 8,703 quantified proteins. The author-deposited dataset description specifies 24 hours at 10 µM. The supplied report instead starts with 8,630 protein features and 173 annotated entries. Its retained sample totals are 1,029 drug-treated, 155 QC, 24 DFO, and 23 DMSO profiles after one unlabeled sample is removed. These differences must be reconciled in A01.

## Work order

| Task | Priority | Depends on | Primary deliverable |
|---|---|---|---|
| A01 Sample provenance and estimability | Required first | Original records | Frozen manifest and a comparison eligibility table |
| A02 Differential and pathway inference | Required | A01 and frozen A03 definitions | Recomputed effect, uncertainty, and multiplicity tables |
| A03 Program definitions and evidence | Required | Identifier mapping from A01 | Versioned program manifest and evidence matrix |
| A04 Reproducibility and sensitivity | Required | A01, A02, A03 | Robustness and null-calibration results |
| A05 Independent aging benchmarks | Central biological analysis | A01, A02, A03 | External reference comparisons and uncertainty |
| A06 Shared stress and growth responses | Central specificity analysis | A01, A02, A03 | Unadjusted and sensitivity-adjusted candidate profiles |
| A07 Candidate mechanisms and pharmacology | Needed for interpretation | A02 through A06 | Auditable shortlist and complete candidate panels |
| A08 Figure and manuscript assembly | Required last | Frozen outputs of A01 through A07 | Five figure packages and claim-to-result ledger |

A03 literature curation and A05 reference discovery can proceed while A01 is running. A02 must wait for the eligibility and definition freezes. A04, A05 scoring, and A06 can then run independently. Candidate selection follows these checks. This ordering prevents different agents from analyzing different sample sets or changing program definitions to favor a candidate.

## A01 Sample provenance and estimability

**Question:** What are the independent experimental units, what was actually measured, and which drug-versus-control contrasts are identifiable?

1. Reconcile library compounds, DFO controls, DMSO controls, pooled QC samples, annotations, and unique chemical entities. Normalize aliases and salts without merging distinct exposures silently. Determine whether DFO is included in the 172 or 173 counts. Verify Blue5_F12 and every exclusion against records.
2. Reconcile 8,703 versus 8,630 versus 7,871 features. Determine whether rows represent proteins, protein groups, or mapped genes. Track mapping losses, duplicates, contaminants, and filtering separately. Check whether zero means missingness.
3. Create one row per observation with `sample_id`, `compound_id`, `condition_role`, `dose`, `duration`, `culture_id`, `culture_day`, `plate`, `well`, `preparation_batch`, `instrument_batch`, `injection_order`, `injection_id`, and exclusion fields. Mark unknown fields explicitly.
4. Distinguish independently treated cultures from multiple wells from the same culture preparation and repeated injections. The parent study calls the replicates biological; verify what was independently repeated. Use the actual independent unit for inference.
5. Produce batch-by-condition and compound-by-batch tables. Test the design matrix rank and contrast estimability. Determine whether each candidate has within-batch DMSO controls. Pay particular attention to RedMix, reported to have no DMSO profiles.
6. Inspect what the pooled QCs are: analytical peptide pools, process controls, or independently prepared biological references. QCs can anchor certain technical shifts; they are not biological DMSO replicates. If reference anchoring is needed, state its assumptions and validate them on batches with genuine controls.
7. Produce an eligibility table classifying contrasts as directly estimable, estimable under an explicit reference-anchoring assumption, or not estimable. Main claims should use supported contrasts. Analyze the reference-anchored subset separately unless equivalence is demonstrated.
8. Compare the planned aging manuscript with the parent manuscript and list reused data, reused methods, reused results, and genuinely new analyses. The main-text novelty statement must match that inventory.

**Return:** `sample_manifest.tsv`, `compound_manifest.tsv`, `protein_mapping.tsv`, `sample_accounting.tsv`, `contrast_eligibility.tsv`, `design_audit.md`, and `parent_study_overlap.tsv`. Include a concise report explaining any unresolved discrepancy. Do not estimate independent replication from injection counts.

## A02 Differential and pathway inference

**Question:** Which compound responses survive a statistical analysis that respects the actual experiment?

1. Freeze the sample universe, contrasts, protein filtering, gene mapping, and gene-set release before testing. Preserve observed missingness. Do not use every absent measurement as evidence of low abundance.
2. Fit an appropriate model to observed log2 intensities with batch and other identifiable design terms. A moderated linear model can be appropriate, but select the model from A01. Aggregate technical injections or model their dependence. If treatment and batch are unidentifiable, do not present a fitted coefficient as a rescued comparison.
3. State the normalization and estimand explicitly. Equal peptide loading or median normalization measures relative proteome composition. It cannot establish increased absolute protein abundance per cell, mitochondrial mass, or protein production without suitable cell-number or spike-in information.
4. Use an observed-data or missingness-aware primary workflow appropriate to the data, with sensitivity to filtering and censoring assumptions. Retain the original left-censored imputation and ComBat workflow as a named sensitivity analysis unless independently justified as primary. Never apply ComBat blindly where biological treatment is confounded with batch.
5. Return protein effect sizes, standard errors or confidence intervals, sample counts, missingness by group, and nominal p values. State whether protein-level BH correction is per contrast or another family. Distinguish this from the pathway-level family.
6. For pathways, use model-based signed statistics and document test type, null, background, set-size rules, and software. Assess inter-protein correlation using an appropriate method, such as a correlation-aware competitive test, alongside preranked enrichment. A gene-permutation GSEA p value and an experimental-label permutation p value do not test the same null.
7. Keep `p_nominal`, software-native `q_native`, and explicitly calculated `q_BH_global` separate. Apply the primary pathway BH correction to all eligible compound-by-program hypotheses, not only hits or highlighted drugs. Define any secondary confirmatory family in advance. Do not call native GSEA q values BH-adjusted unless they were calculated that way.
8. Use sufficiently resolved p values for the intended testing family. If using permutations, state their number and use a nonzero finite-sampling convention. For very small p values, an appropriate adaptive method may be needed. Use treatment-label permutations only within genuinely exchangeable experimental blocks and at the independent unit.
9. Report nonsignificant and non-estimable results as well as positive results. Recompute the source claims of 165 responsive treatments and 999 pairs only after the compound denominator and correction family are resolved.

**Return:** `protein_effects.tsv`, `program_effects.tsv`, `contrast_definitions.json`, `analysis_config.yaml`, `session_info.txt`, scripts or notebooks, and a short `inference_report.md`. Every result row needs a pipeline identifier, compound identifier, contrast identifier, program identifier, effect direction, counts, and test provenance.

## A03 Program definitions and evidence

**Question:** What does each protein program actually measure, and how directly is it connected to aging biology?

1. Replace one-generic-pathway-per-hallmark labeling with a versioned, evidence-graded program collection. Multiple complementary modules per process may be warranted, but keep the primary collection compact and freeze it before rerunning inference. Treat additional modules as exploratory.
2. Record `program_id`, exact source set ID, release, retrieval date, source DOI or URL, species, identifiers, membership checksum, biological description, aging evidence, directness of the readout, and whether a beneficial direction is justified in this context. Distinguish pathway activity hypotheses from protein abundance programs.
3. Include measured-set size, total set size, measured coverage, mapping loss, and gene overlap. Prespecify a minimum measured-set size and evaluate sensitivity; a set with five measured proteins should not be called a comprehensive hallmark assay. Use a measured-proteome background where appropriate.
4. Separate DNA repair proteins from DNA damage burden; chromatin programs from epigenetic state; telomere maintenance from telomere length; autophagy machinery from flux; respiratory-complex abundance from respiratory capacity; and senescence-associated expression from stable senescence.
5. Call TLR cascades host innate sensing. Do not use this as a dysbiosis measurement. Do not claim stem cell exhaustion from HepG2 stemness annotations. Explicitly identify the hallmarks that this cellular model cannot assess.
6. Distinguish the MSigDB HALLMARK collection from the hallmarks of aging. Treat DEMAGALHAES_AGING_UP and DEMAGALHAES_AGING_DN as transcript-derived aging priors. Keep both directions, their overlap, species mapping, and any cell-context mismatch visible.
7. Quantify pairwise membership overlap and overlap among leading-edge proteins for important candidate responses. Determine whether several named programs are driven by the same small protein subset. Repeat candidate interpretation after removing shared contributors or using complementary nonredundant modules, where scientifically justified.

**Return:** `program_manifest.tsv`, `program_membership.tsv`, `programs.gmt`, `program_evidence_matrix.tsv`, `program_overlap.tsv`, `leading_edge_overlap.tsv`, and `annotation_decisions.md`. The evidence matrix must state what functional assay would be needed for a stronger interpretation.

## A04 Reproducibility and sensitivity

**Question:** Do the candidate patterns persist under defensible analysis choices and independent units?

1. Assess QC stability separately from biological consistency. Report protein missingness, QC variability, drift versus injection order, DMSO behavior, and repeated DFO profiles. A mixed PCA plot is not evidence that all batch effects are eliminated.
2. Evaluate treatment-effect consistency at the level supported by the design. Use repeated conditions across batches as genuine cross-batch checks. For compounds confined to one batch, a leave-one-batch-out check cannot establish drug-level replication. Clearly distinguish leave-one-culture-out from leave-one-injection-out.
3. Compare supported batch-aware primary results with the source pipeline, alternative missingness handling, exclusion of RedMix, alternative normalization, and reasonable gene-set definitions. State why each sensitivity is appropriate; do not search a large grid for the most favorable pipeline.
4. Use permissible control-versus-control splits or blocked label permutations to assess null behavior. Treat reused split controls as a diagnostic, not independent validation. If controls are too few for meaningful calibration, report that limitation.
5. Report effect-direction agreement and rank stability, with uncertainty, for all candidates. Use bootstrap resampling at the independent experimental unit. Treat model-based intervals from one experimental run as within-design uncertainty, not independent reproducibility.

**Return:** `qc_metrics.tsv`, `replication_metrics.tsv`, `sensitivity_effects.tsv`, `null_calibration.tsv`, `candidate_stability.tsv`, and `robustness_report.md`. Identify which original conclusions change.

## A05 Independent aging benchmarks

**Question:** Do drug-induced protein responses oppose independently observed cellular aging changes?

1. Build a reference inventory with priority given to age-associated protein effects in primary cultured cells. Start with Tsitsipatis et al. 2022 and Waldera-Lupa et al. 2014. Use donor-level age associations where data permit, accounting for available covariates. Their observational donor-age associations are not validated cellular aging rates.
2. Add induced-senescence proteomes and chemical-reprogramming proteomes as separate reference classes. Basisty et al. secretomes represent a different compartment from whole-cell lysates; do not combine them into one interchangeable reference. Mitchell et al. reprogramming effects are intervention responses, not a generic aging ground truth.
3. Record accession, species, cell type, compartment, aging or senescence definition, induction method, comparison, independent n, protein coverage, and processing. State why each reference is biologically relevant to HepG2 and where it is not.
4. Compare the signed drug effect with the signed age-associated effect on measured shared proteins. Use continuous-effect agreement or opposition and signed up/down enrichment where appropriate. Report shared-feature counts, confidence intervals, corrected significance, and sensitivity to protein overlap. Avoid judging opposition solely by a negative correlation p value across correlated proteins.
5. Use an appropriate empirical null or correlation-aware inference to account for protein dependence and selection. Clearly distinguish uncertainty in the drug effect from uncertainty in the external age effect. Show reference-specific results before considering a combined score; combine only compatible references and state weights before scoring candidates.
6. Evaluate transcript-derived signatures separately as supportive priors. Do not call a projected RNA signature an independently measured proteomic aging clock.
7. Freeze reference inclusion and scoring before final candidate selection. Because the current candidates have already been inspected, characterize this as retrospective external benchmarking. Do not claim a prospective blinded discovery.
8. Where positive and negative intervention comparators exist, test separation with uncertainty. AUC is appropriate only with sufficient independently labeled compounds and an honest split. Hold out entire compounds or related chemotypes when training any score. Do not train an age clock on this dataset, which has no age labels.

**Return:** `external_reference_manifest.tsv`, reference-specific signed effect tables, `drug_reference_comparisons.tsv`, `reference_coverage.tsv`, code, and `aging_benchmark_report.md`. Report lack of agreement as a result. A reference that does not map adequately should be excluded with a reason rather than forced into a score.

## A06 Shared stress and growth responses

**Question:** Do proposed candidates carry information beyond a generic change in proteome composition or proliferation?

1. Recover measured cell counts, viability, cell size, total protein per cell, or growth data when they exist. Analyze these separately from protein-derived proxies. If they do not exist, do not state that cytostasis or toxicity has been ruled out.
2. Quantify proliferation and cell-cycle programs, apoptosis, stress responses, and global perturbation magnitude using frozen definitions. Report their relationships to aging-associated programs and candidate rankings. The same proteins must not define both a nuisance score and its supposedly independent target signal without a documented sensitivity analysis.
3. Describe dominant shared-response axes across compounds and determine whether candidates are ordinary or unusual relative to matched perturbation strength and growth effects. A screen-relative outlier is not automatically beneficial.
4. Present unadjusted results alongside sensitivity analyses using measured covariates or clearly labeled molecular proxies. Growth and stress changes can mediate a real drug effect; automatic residualization can remove biology. State the causal interpretation and avoid presenting adjusted associations as causal estimates.
5. Evaluate relative-composition effects. Apparent increases in mitochondrial proteins may reflect decreases in other abundant protein classes under equal-loading normalization. Check mitochondrial versus cytosolic protein patterns and any available per-cell measurements. Without absolute measurements, retain a relative-abundance interpretation.
6. Use negative comparators with similar broad perturbation or antiproliferative effects where possible. Test whether binimetinib or baicalein remains distinctive without selectively choosing the comparator after viewing the result.

**Return:** `shared_response_scores.tsv`, `measured_growth_metadata.tsv` if available, `candidate_specificity.tsv`, `composition_sensitivity.tsv`, and `specificity_report.md`. State the strongest conclusion possible if growth data are unavailable.

## A07 Candidate mechanisms and pharmacology

**Question:** Which candidates merit follow-up, and what specific molecular hypotheses do they support?

1. Freeze nomination criteria using evidence from A02 through A06. Keep discovery and validation criteria distinct. Do not select only compounds with favorable anecdotes in the longevity literature.
2. Revisit binimetinib and baicalein as exploratory candidates, not predetermined winners. Include NDGA, lonafarnib, DFO, alpha-lipoic acid, and statins where relevant to the audited data. DFO is an assay response control, not automatically a positive geroprotector control.
3. Compare binimetinib with erlotinib and mobocertinib using complete prespecified respiratory-complex and inflammatory modules. Do not select only the most significant proteins. Identify which proteins are actually leading-edge members and which are supporting pathway members.
4. Use observed log2 abundance, effect size, confidence interval, missingness, and batch identifiers. Do not use DMSO-standardized values as the sole evidence. Provide raw observations and independent-unit summaries with n defined in the legend.
5. Verify drug identities, intended targets, likely off-target effects at the tested concentration, and whether a drug class is represented by enough compounds for class inference. Evidence that one MEK inhibitor extends lifespan in one organism does not establish that binimetinib does so, or that the HepG2 phenotype has the same mechanism.
6. Evaluate whether the program representation adds information over generic proteome response. Compare compound similarity in the complete proteome and selected program space. If compound classes are large enough, assess within-class similarity using plate-aware permutations or other appropriate controls. Describe classes with one member descriptively.
7. Return a shortlist with favorable evidence, contradictory evidence, unresolved confounders, and one falsifiable next experiment per candidate. If no candidates survive the checks, report that rather than lowering thresholds after inspection.

**Return:** `candidate_evidence.tsv`, `candidate_protein_effects.tsv`, `candidate_selection_rules.yaml`, `drug_similarity.tsv`, a pharmacology evidence table with primary citations, and `candidate_report.md`.

## A08 Final figures and manuscript assembly

Build from frozen result tables. Do not redraw a source chart by tracing its image or preserve erroneous hallmark labels for visual continuity. Show complete tested universes in supplemental data and clearly state missing or non-estimable results.

| Figure | Main question | Required panels and source tasks |
|---|---|---|
| **1 Dataset provenance and performance** | Are the comparisons supported by the experiment? | A: Existing screen and new reanalysis workflow. B: Sample accounting and independent units. C: Batch-by-condition coverage, highlighting RedMix. D: QC and control performance plus permitted replication checks. A01, A04. |
| **2 Evidence and coverage** | Which aging-associated programs are measured and how directly? | A: Evidence-graded program-to-process matrix, marking unmeasured states. B: Measured versus total membership. C: Gene-set and leading-edge overlap. A03. |
| **3 Compound response atlas** | What distinct combinations do drugs produce? | A: Every eligible compound by program with signed effects and audited significance. B: Selected interpretable contrasts with intervals. C: Shared response axes and candidate specificity. Put the full readable atlas in supplement if the main figure is too dense. A02, A04, A06. |
| **4 Independent aging benchmarks** | Are candidate effects related to aging beyond annotation? | A: External reference inventory and coverage. B: Drug-by-reference agreement or opposition. C: Candidate uncertainty and context sensitivity. Separate proteomic and transcript-derived references. A05, A06. |
| **5 Candidate mechanisms** | What testable biology emerges? | A: Prespecified compound comparison. B: Complete respiratory-complex and inflammatory-module effects. C: Observed values for prespecified sentinel proteins. D: Independent replication or functional validation if generated; otherwise show robustness and label the figure as nomination. A04, A07. |

Required supplementary figures: complete compound atlas including null results; detailed QC and drift; missingness and preprocessing sensitivity; control-null calibration; program overlap and leading edges; external reference context and coverage; and complete candidate protein panels. Merge related panels when this improves readability rather than creating one figure per diagnostic.

Every figure package must include editable SVG or PDF, a 300-dpi PNG preview, the exact source table, generation code, and a complete legend. Define n, biological versus technical replication, contrast direction, error bars, normalization, testing method, and correction family. Use consistent compound and program names and accessible colors. No em dash character in manuscript text or legends.

Create `claim_result_ledger.tsv` with claim text, result table, row identifiers, statistical family, figure panel, limitations, and final manuscript location. Create `results_summary.md` with verified counts and exact estimates for replacing the bracketed manuscript sections. Retain all changed or negative findings. No new figure or abstract claim should lack a ledger entry.

## Experiments that analysis agents cannot supply

These are separate from finishing the computational reanalysis. Do not invent them, imply that they are completed, or make their results a prerequisite for reporting a well-qualified resource analysis.

| Experiment | What it would resolve | Claim enabled if successful |
|---|---|---|
| **Independent repeat of a small prespecified candidate and comparator panel** on new culture days, randomized across plates with DMSO and DFO | Whether nominated signatures reproduce beyond the original run | Reproducible pharmacological protein programs |
| **Dose response with measured growth, viability, and cell number** | Whether signatures persist at tolerated exposures and whether relative loading or growth explains them | Better-supported candidate specificity |
| **Respiration assay with suitable normalization and mechanism controls** for a mitochondrial candidate | Whether respiratory-complex abundance corresponds to respiratory function | Modulation of measured mitochondrial function |
| **Autophagic flux assay** if autophagy is a central claim | Whether abundance changes correspond to altered throughput | Modulation of measured autophagic flux |
| **Relevant primary-cell aging or induced-senescence model** with vehicle, young or nonsenescent reference, senescent or old control, and appropriate intervention comparators | Whether candidates improve a defined aging-related phenotype and whether apparent benefit reflects selective cell killing | Functional rescue or senomorphic/senolytic activity within that model |
| **A second biologically relevant cell model** | Whether the response is restricted to transformed HepG2 cells | Broader cellular generality, within the tested models |

For any senescence claim, use a panel of complementary measurements, such as growth arrest and recovery, DNA damage, SA-beta-gal, and secreted factors, with cell counts and viability. Distinguish senolysis from phenotypic normalization. Do not infer organismal geroprotection from one cell assay. Each experiment needs an explicit primary endpoint and a sample-size justification based on the expected effect and variability, not an arbitrary universal replicate count.

## Decision rules for the manuscript

- If only the current pathway annotations remain after reanalysis, present a limited secondary molecular analysis. Evaluate whether it belongs within the parent study rather than claiming an independent aging platform.
- If robust compound-specific patterns also agree with independent aging references and exceed shared response explanations, the stronger contribution is an aging-focused pharmacological map with experimentally grounded candidate nomination.
- If prespecified candidate patterns reproduce and functional endpoints improve in aging-relevant cells, the manuscript can make those specific functional claims. A general geroprotective screening platform would still require demonstrated predictive performance across multiple interventions.

## Literature framing to retain

| Precedent | What it already establishes | What the present study must add |
|---|---|---|
| Waldera-Lupa et al. 2014, https://doi.org/10.18632/aging.100698 | Donor-derived cellular aging proteomics linked to several aging processes | Comparative pharmacological responses tied to independent aging evidence |
| Tsitsipatis et al. 2022, https://doi.org/10.1111/acel.13609 | Deep proteomic profiling of primary fibroblasts across adult donor ages | Drug response combinations and credible benchmark agreement |
| Basisty et al. 2020, https://doi.org/10.1371/journal.pbio.3000599 | DIA characterization of heterogeneous senescence-associated secretomes | Explicit compartment-aware comparisons and broader intracellular programs |
| Mitchell et al. 2024 version of record, https://doi.org/10.7554/eLife.90579 | Proteomic and phosphoproteomic responses to chemical reprogramming with functional follow-up | Comparison across a larger perturbation collection and specificity analysis |
| Yoo et al. 2026 preprint, https://doi.org/10.64898/2026.02.26.708335 | Spatial protein remodeling explicitly connected to aging hallmarks in yeast | Mammalian pharmacological comparisons with transparent limits |
| Gross et al. 2026, https://doi.org/10.1038/s43587-026-01161-8 | Hallmark-based drug prioritization using networks and transcriptional signatures | Experimental protein response evidence and appropriate external validation |
| Jiang et al. 2026 preprint, https://doi.org/10.64898/2026.05.04.722718 | The actual screen, its scale, robotic workflow, and general mechanism-discovery platform | A new biological question, new analyses, and new evidence specific to aging |

## Common instructions for all analysis agents

Use the frozen manifests and shared identifiers. Do not change another task's input or output files silently. Return code, parameters, software versions, input checksums, exclusions, and the requested tables even for null results. Distinguish performed analyses from proposals. Do not manufacture results from the preliminary report. Use the report as a list of hypotheses to check. Report any input that blocks the analysis and complete the unaffected work. Do not optimize thresholds, reference selection, or program definitions to retain binimetinib or baicalein.

One useful handoff message is: "Complete task Axx in this brief using the frozen upstream outputs. Return the specified files, a concise findings report, and explicit limitations. Keep all numerical claims traceable to generated result tables and do not present planned validation as completed."
