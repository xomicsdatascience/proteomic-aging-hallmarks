"""Build manuscript prose, legends and claim ledger from the current result tables."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from analysis_common import ROOT, PRIMARY_IDS, PROGRAM_NAMES, HIGHLIGHTS, read_table, write_table, json_write, short_drug
from prioritization import ASSAY_DIRECTIONS


def number(v,digits=3):
    if not np.isfinite(v):return 'not estimable'
    return f'{v:.{digits}g}'


def main():
    config=json.loads((ROOT/'reanalysis/A02/analysis_config.json').read_text())
    diagnostic=json.loads((ROOT/'reanalysis/A04/diagnostic_summary.json').read_text())
    inputs=json.loads((ROOT/'reanalysis/A01/input_summary.json').read_text())
    p=read_table('reanalysis/A02/program_effects.tsv');primary=p[p.program_id.isin(PRIMARY_IDS)]
    indices=read_table('reanalysis/A06/exploratory_prioritization.tsv').set_index('compound_id')
    eligible=indices[indices.four_axis_eligible_rank.notna()].sort_values('four_axis_eligible_rank')
    four=read_table('reanalysis/A08/source_tables/Supplementary_Table_S2_four_axis_candidates.tsv')
    four_primary=primary[primary.program_id.isin(ASSAY_DIRECTIONS)]
    all_four_q=int(four_primary.assign(hit=(four_primary.plate_competitive_log2fc*four_primary.program_id.map(ASSAY_DIRECTIONS)).gt(0)&four_primary.q_BH_global_plate.lt(.05)).groupby('compound_id').hit.all().drop(index='DFO',errors='ignore').sum())
    top=eligible.head(5);first=top.index[0]
    man=read_table('reanalysis/A03/program_manifest.tsv').set_index('program_id')
    plate=read_table('reanalysis/A02/plate_program_observations.tsv')
    proteins=read_table('reanalysis/A07/candidate_protein_effects.tsv')
    growth=json.loads((ROOT/'reanalysis/A06/growth_correlations.json').read_text())
    reference=read_table('reanalysis/A04/vehicle_reference_sensitivity.tsv')
    pair=read_table('reanalysis/A02/direct_comparator_program_effects.tsv')
    key=p.set_index(['compound_id','program_id']);b=key.loc['Binimetinib'];ba=key.loc['Baicalein']
    prot=proteins.set_index(['compound_id','feature_id'])
    ntests=len(primary);n_sig=int((primary.q_BH_global_plate<.05).sum());n_drugs=primary.loc[primary.q_BH_global_plate<.05,'compound_id'].nunique()
    cam_sig=int((primary.q_BH_global_camera<.05).sum());emp_sig=int((primary.q_BH_global_camera_emp<.05).sum())
    highlight=[d for d in HIGHLIGHTS if d in set(primary.compound_id)]
    clauses=[]
    for drug in ['Binimetinib','Baicalein']:
        for pid in ['PROG_OXPHOS','PROG_INFLAMMATORY']:
            r=key.loc[(drug,pid)]
            clauses.append(f"{drug} {PROGRAM_NAMES[pid]}: competitive mean log2FC {r.plate_competitive_log2fc:+.3f} (95% CI {r.plate_ci95_low:+.3f} to {r.plate_ci95_high:+.3f}), q = {number(r.q_BH_global_plate)}")
    selected=primary[primary.compound_id.isin(['Binimetinib','Baicalein'])][['compound_id','program_id','plate_competitive_log2fc','plate_ci95_low','plate_ci95_high','plate_t','q_BH_global_plate','n_effective_plates']]
    write_table(selected,ROOT/'reanalysis/A08/source_tables/selected_response_estimates.tsv')
    sensitivities=reference[reference.compound_id.eq('Binimetinib')&reference.program_id.isin(['PROG_OXPHOS','PROG_INFLAMMATORY'])]
    robust=all((row.plate_competitive_log2fc>0 if row.program_id=='PROG_OXPHOS' else row.plate_competitive_log2fc<0) for row in sensitivities.itertuples())
    sensitivity_text=('The mitochondrial-positive and inflammatory-negative directions were retained when the first listed vehicle well, second listed vehicle well, or their mean was used on the same paired-control plates.' if robust else 'Some directions were sensitive to which vehicle well was used; the reference-specific results are reported without selecting the most favorable control.')
    direct=pair[pair.compound_id.eq('Binimetinib_vs_Erlotinib (Hydrochloride)')&pair.program_id.eq('PROG_INFLAMMATORY')].iloc[0]
    direct_text=(f"On {int(direct.n_effective_plates)} shared plates, the direct binimetinib-versus-erlotinib inflammatory-program contrast was {direct.plate_competitive_log2fc:+.3f} log2 units (95% CI {direct.plate_ci95_low:+.3f} to {direct.plate_ci95_high:+.3f}; q = {number(direct.q_BH_global_plate)} within the separate comparator family)." if direct.n_effective_plates>=3 else 'The displayed kinase inhibitors do not all share assay plates. Direct additive-model comparisons are supplied separately; a plate-paired estimate is reported only when at least three shared plates are available.')
    index_results=(f"The exploratory 14-program combined response index ranked binimetinib first ({indices.loc['Binimetinib','combined_response_index']:+.3f}) and baicalein second ({indices.loc['Baicalein','combined_response_index']:+.3f}) among the 160 library entries (Figure 2B; Supplementary Table S3). They followed only {int(indices.loc['Binimetinib','combined_response_index_n_requested_directions'])} and {int(indices.loc['Baicalein','combined_response_index_n_requested_directions'])} of the 14 target signs, respectively. Their large changes on some axes outweighed opposing changes on others. This ranking therefore does not mean coordinated improvement of all programs. The raw index and E2F/G2M residual are displayed separately, with no index-level significance or rank-certainty claim. Different target signs, scales, weights or program redundancy can change this ordering.")
    top_text=', '.join(f"{short_drug(drug)} ({row.four_axis_index:.3f})" for drug,row in top.iterrows())
    first_effects=four[four.compound_id.eq(first)].set_index('program_id')
    first_text='; '.join(f"{PROGRAM_NAMES[pid]} {first_effects.loc[pid,'plate_competitive_log2fc']:+.3f} log2 units, q = {number(first_effects.loc[pid,'q_BH_global_plate'],6)}" for pid in ASSAY_DIRECTIONS)
    positions=', '.join(f"{short_drug(drug)} at {row.assay_group} {row.well}" for drug,row in top.head(3).iterrows())
    vignette=f"""### A four-axis profile for cellular follow-up

Restricting the target to higher OXPHOS, macroautophagy machinery and ubiquitination/proteasome responses with lower senescence-associated machinery produced a different shortlist (Supplementary Figure S5; Supplementary Table S2). Of 160 library entries, {len(eligible)} had primary point estimates in all four requested directions. Among those entries, the equal-weight four-axis index ranked {top_text}. These are descriptive candidates for the specified abundance profile, not validated functional responses. The unrestricted mean and the sign-filtered rank are both exported, because a strong response on one axis can otherwise conceal an opposite response on another.

For the leading entry, {short_drug(first)}, the component effects were {first_text}. The senescence estimate's unadjusted 95% interval excluded zero, but its global q value remained above 0.05. Across the screen, {all_four_q} library entries met q < 0.05 on all four axes in the requested directions. Thus the {len(eligible)}-compound sign filter is not evidence for a jointly significant four-axis response. The per-axis intervals, q values, effective plate counts, weakest signed components and complete rankings are provided rather than treating a rank as confidence in superiority.

The first three filtered candidates also share a potentially consequential layout feature: {positions}. All occupy column 12 of the Blue plates. Compound and coordinate effects cannot be separated in this design, so repositioned randomized replication is needed before attributing this pattern to treatment. Binimetinib and baicalein did not pass the four-axis sign filter because their primary proteasome-associated responses were negative. Their mitochondrial/inflammatory profile addresses a different hypothesis from the four-axis assay profile.

Cellular follow-up should test respiration, autophagic flux, proteasome activity and senescence phenotypes alongside cell counts and viability, using dose and time series with distributed vehicle wells. The present measurements concern relative protein abundance and do not establish improved function on any of these axes.
"""
    table_lines=['| Rank | Library entry | Four-axis index | Fixed position | OXPHOS effect; q | Autophagy effect; q | Proteasome effect; q | Senescence effect; q |',
                 '|---:|---|---:|---|---|---|---|---|']
    for drug,row in top.iterrows():
        values=four[four.compound_id.eq(drug)].set_index('program_id')
        cells=[f"{values.loc[pid,'plate_competitive_log2fc']:+.4f}; {number(values.loc[pid,'q_BH_global_plate'],6)}" for pid in ASSAY_DIRECTIONS]
        table_lines.append('| '+ ' | '.join([str(int(row.four_axis_eligible_rank)),drug,f'{row.four_axis_index:.4f}',f'{row.assay_group} {row.well}']+cells)+' |')
    (ROOT/'docs/FOUR_AXIS_VIGNETTE.md').write_text('# Four-axis cellular follow-up candidates\n\n'+vignette.split('\n',2)[2]+'\n'+ '\n'.join(table_lines)+'\n\nEffects are primary matched-plate competitive mean log2 fold changes. q values use the existing global 2,254-test primary family. The ranking has no p value or q value. [Shortlisted axis estimates and 95% intervals](../reanalysis/A08/source_tables/Supplementary_Table_S2_four_axis_candidates.tsv) and [index definitions](../reanalysis/A06/prioritization_definitions.tsv) accompany [Supplementary Figure S5](../reanalysis/A08/figures/figS5_four_axis_candidates.png).\n')
    abstract=f"""Deep perturbation proteomics can compare coordinated protein responses across cellular programs related to aging, but abundance changes do not directly establish pathway function or rejuvenation. We reanalyzed a published HepG2 screen of 160 library compounds and a deferoxamine assay control, retaining {inputs['retained_features_50pct']:,} protein features. Fourteen selected programs showed {man.loc[PRIMARY_IDS,'coverage_pct_retained'].min():.0f}% to {man.loc[PRIMARY_IDS,'coverage_pct_retained'].max():.0f}% measured membership. Primary pathway tests used competitive mean log2 fold changes across matched plates, with proteins treated as correlated measurements within each plate. Across {ntests:,} compound-program hypotheses, {n_sig:,} pairs met global Benjamini-Hochberg q < 0.05 under this model. Binimetinib combined increased relative mitochondrial protein abundance with a negative inflammatory-annotation response, whereas comparator kinase inhibitors showed differing combinations. Baicalein showed a related response. These examples also altered proliferation-associated proteins and produced mixed changes in RNA-derived aging priors. Vehicle-well differences limit calibration claims, and fixed treatment coordinates confound compound and position effects. This analysis provides a reproducible resource for comparing protein-program responses and selecting functional follow-up experiments. It does not establish reduced biological age, organismal geroprotection, or clinical benefit."""
    text=f"""# Proteomic profiling distinguishes drug responses across aging associated cellular programs

Jesse G. Meyer

## Abstract

{abstract}

## Introduction

Aging involves interacting molecular and cellular processes, including genome maintenance, proteostasis, nutrient sensing, mitochondrial remodeling, and inflammatory signaling (Lopez-Otin et al., 2023). These processes are often studied using individual functional assays. Deep protein profiling offers a complementary view of how many components respond to a perturbation in the same experiment. Its interpretation requires separating the abundance of pathway-associated proteins from the function or biological state those proteins participate in. For example, more autophagy machinery does not establish increased autophagic flux, and a DNA-repair annotation does not measure accumulated DNA damage.

Previous work has connected proteomic changes with donor age in human fibroblasts (Waldera-Lupa et al., 2014; Tsitsipatis et al., 2022), defined heterogeneous senescence-associated secretomes (Basisty et al., 2020), and characterized chemical reprogramming with proteomic and functional measurements (Mitchell et al., 2024). Spatial protein profiling has linked remodeling to aging hallmarks in yeast (Yoo et al., 2026). Network and perturbational expression approaches have also prioritized hallmark-related drug candidates (Gross et al., 2026). Consequently, neither the association of proteins with aging pathways nor hallmark-based drug prioritization is itself a new premise.

Jiang and colleagues generated a robotic DIA-MS perturbation dataset in HepG2 cells, with 172 library compounds and 1,232 proteomic observations (Jiang et al., 2026). The parent study established the experimental platform, drug-response organization, and selected validation examples. Here we reuse the deposited matrix to compare responses across a defined set of aging-associated protein programs. We ask which parts of those annotations are observed, how drug response combinations differ, and whether a mitochondrial/inflammatory response example can be distinguished from generic proliferation-related changes. The present analysis measures relative protein responses in a transformed cell line; it does not estimate biological age or rank compounds by demonstrated geroprotection.

## Methods

### Data provenance and experimental scope

The analyzed export contains {inputs['raw_samples']:,} observations and {inputs['raw_features']:,} gene-symbol-indexed protein features. This feature count differs from the 8,703 reported in the parent manuscript; accession-level export provenance is not provided in the deposited file. All counts here refer to the analyzed export. One unlabeled well was excluded before feature filtering. Twelve plates contained repeated library treatments, DMSO vehicle controls, and deferoxamine (DFO) assay-response controls. Primary comparisons comprise 160 library entries plus DFO. A plate without within-plate vehicle controls was excluded from primary inference. The primary model contains {config['n_primary_samples']:,} assay observations: 957 library-treated wells, 24 DFO wells, and 23 DMSO wells. Analytical pooled QC injections were excluded from biological inference and evaluated separately (Supplementary Figure S1).

Library treatments were 10 micromolar for 24 hours, as described in the parent study. Drug-treated wells were replicated across five or six assay plates; DFO was present on 12 plates. The parent study describes biological replicates. Independent culture days and injection order cannot be verified from the export, so plate/well identifiers are not interpreted as separately documented culture batches. Original compound names and canonical aliases are retained in a mapping table.

### Preprocessing and program definitions

After excluding the unlabeled well, features with more than 50% zero values across the remaining observations were removed, retaining {inputs['retained_features_50pct']:,} features. Zero values were treated as missing. Observed intensities were log2-transformed, and each sample was median-normalized to the dataset reference median. No missing abundance was imputed. This normalization measures relative proteome composition, not absolute protein abundance per cell. A stricter 20% missingness universe was used for sensitivity analysis.

The committed gene-set file defines 14 primary programs: 12 mechanistic protein annotations and two RNA-derived mammalian aging priors from de Magalhaes and colleagues (2009). Their sources are MSigDB Hallmark and Reactome libraries, with membership, retrieval provenance, and checksums supplied. Three additional modules describe E2F proliferation, G2M checkpoint, and apoptosis-associated responses. The annotation panel was fixed for this reproducible reanalysis, rather than prospectively registered before the original exploratory work.

The selected programs do not directly assay all organismal hallmarks. Telomere-maintenance proteins are not telomere length, chromatin-remodeling proteins are not an epigenetic age, whole-cell inflammatory proteins are not secreted SASP, and host Toll-like receptor proteins are not microbiome dysbiosis. No universal beneficial direction is inferred from these annotations. Target signs used for exploratory prioritization are explicit hypotheses, not evidence of geroprotection.

### Protein models and observation eligibility

Observed log2 abundances were fitted to an additive compound-plus-plate linear model with DMSO as the compound reference. The design has rank {config['design_rank']}. Missing observations were omitted per feature and non-estimable coefficients were left unavailable. Residual variances were moderated across features using a parametric empirical Bayes prior, following Smyth (2004), with prior degrees of freedom {config['d0_prior']:.3f} and prior variance {config['s0_2_prior']:.5f}. Protein-level tests require at least three observed treated and three observed reference values. Model estimates, standard errors, 95% intervals, available observation counts, eligibility, and per-contrast BH-adjusted q values are exported. Untested or unavailable protein results are not replaced with zero.

### Primary pathway tests using matched plates

For each compound and program, eligible protein fold changes were first calculated relative to the vehicle mean on the same plate. Within each plate, the mean log2 fold change of available program members was contrasted with the mean of available eligible proteins outside that program. At least five program members and 100 background proteins were required for a plate score. A two-sided one-sample t test of these competitive mean differences across matched plates was then used, with degrees of freedom equal to the number of effective plates minus one and at least three effective plates required. Thus proteins within a plate are not counted as independent biological replicates. A positive result denotes a larger relative response in the program than its measured background, not necessarily increased pathway activity.

BH correction was applied across the full {ntests:,}-test primary family. The three growth/stress annotations and the five direct kinase comparisons were separate families. Unavailable tests remain explicitly unavailable and are conservatively padded when accounting for the declared family. All primary compound-program results, plate-level observations, intervals, and available feature counts are provided in Supplementary Table S1. Parametric inference assumes independent, approximately normally distributed plate responses; the small number of plates and incomplete culture provenance limit this assumption.

### Secondary tests and sensitivity analyses

A Python implementation of the non-rank CAMERA-PR competitive formula was applied to precomputed moderated protein t statistics (Wu and Smyth, 2012). This is distinct from calling the expression-matrix limma camera function. Fixed inter-protein correlation 0.01 and pairwise empirical residual-correlation versions were evaluated separately. Empirical correlations require at least 20 common residual observations per protein pair; inadequate estimates remain unavailable rather than defaulting to zero. Each contrast uses its actual finite protein background. These analyses weight protein responses differently from the primary mean-fold-change test and have different assumptions; more discoveries do not demonstrate better calibration.

Pre-ranked GSEA used 2,000 gene-set permutations and seed 42 to provide normalized enrichment scores (NES) and leading-edge membership (Subramanian et al., 2005). NES is descriptive here. Native permutation p values and GSEA FDR values are not used as confirmatory significance estimates. The stricter missingness analysis re-estimated variance moderation within its feature universe. Within-plate vehicle splits, exhaustive paired vehicle sign assignments, and first-versus-second vehicle reference analyses were diagnostic checks. Bootstrap intervals from 1,000 plate resamples describe program mean fold-change stability, rather than a refitted geroprotection rank.

Across-drug residualization on E2F response and mean absolute protein t statistic was a descriptive sensitivity analysis. Those post-treatment molecular responses are not measured baseline confounders. Residuals do not establish causal independence from growth. Comparisons of mitochondrial and cytosolic ribosomal protein changes investigate one relative-composition explanation and cannot establish absolute organelle abundance.

### Exploratory combined response and four-axis indices

The 14-program combined response index retains the earlier analysis's target-sign hypotheses while replacing its inferential score. Positive targets are OXPHOS, aging-down RNA prior, macroautophagy machinery, unfolded protein response, ubiquitination/proteasome, DNA repair, telomere maintenance and chromatin organization. Negative targets are inflammatory response, host innate TLR sensing, senescence-associated machinery, mTORC1 response, aging-up RNA prior and cell-cell communication. These are exploratory choices, not universal beneficial directions. For each program, the primary competitive mean log2 fold change was divided by its sample SD across the 160 library entries, multiplied by the target sign, and averaged with equal weight across all 14 components. No screen mean was subtracted, so zero retains the treatment-versus-vehicle competitive-effect anchor. Scaling equalizes observed program variation and is not a null-variance estimate. Related or overlapping programs can still give the same biology repeated influence.

The four-axis index uses the same procedure, with equal weights only for positive OXPHOS, macroautophagy and proteasome responses and negative senescence-associated machinery responses. A requested-profile ranking additionally requires all four primary point estimates to have those signs. This sign filter does not require or imply joint significance. Missing components leave an index unavailable. DFO is displayed as an assay control but excluded from SD estimation, regression fits and ranks. Directions, weights, scales and every component estimate are exported.

For each index, an intercept-plus-E2F/G2M linear regression fitted across library entries supplies a separate post-treatment-proxy residual and rank. The raw index and residual are never combined. Neither index is assigned a normal-theory p value, BH q value, significance boundary or null covariance from gene-set overlap. The indices rank observed abundance profiles within this screen; they do not measure biological age, geroprotection, causal independence from growth or certainty in rank ordering.

## Results

### Protein coverage across the selected programs

Measured membership spanned {man.loc[PRIMARY_IDS,'coverage_pct_retained'].min():.0f}% to {man.loc[PRIMARY_IDS,'coverage_pct_retained'].max():.0f}% across the 14 programs (Figure 1). Coverage included {int(man.loc['PROG_OXPHOS','measured_retained_7871'])}/{int(man.loc['PROG_OXPHOS','total_genes'])} oxidative-phosphorylation genes, {int(man.loc['PROG_TELOMERE_MAINT','measured_retained_7871'])}/{int(man.loc['PROG_TELOMERE_MAINT','total_genes'])} telomere-maintenance genes, and {int(man.loc['PROG_DNA_REPAIR','measured_retained_7871'])}/{int(man.loc['PROG_DNA_REPAIR','total_genes'])} DNA-repair genes. Coverage was lower for the inflammatory response ({int(man.loc['PROG_INFLAMMATORY','measured_retained_7871'])}/{int(man.loc['PROG_INFLAMMATORY','total_genes'])}) and cell-cell communication ({int(man.loc['PROG_CELL_CELL_COMM','measured_retained_7871'])}/{int(man.loc['PROG_CELL_CELL_COMM','total_genes'])}), consistent with whole-cell lysate sampling. These are annotation-coverage counts across the retained matrix; the number of eligible observed members in a particular contrast can be smaller.

![Figure 1](figures/fig1_program_coverage.png)

**Figure 1. Measured protein coverage of selected aging-associated programs.** Grey markers show annotated gene-set size and blue markers show measured membership after the specified feature filter. Labels give measured/annotated counts and percentages. Coverage does not imply a functional measurement of the corresponding biological state. RNA-derived priors are identified separately by name.

### Drug responses form distinct program combinations

Across {ntests:,} primary hypotheses, the matched-plate analysis identified {n_sig:,} pairs in {n_drugs} contrasts at global BH q < 0.05. The selected heatmap illustrates {len(highlight)} retrospectively chosen response examples (Figure 2), with the complete 161-contrast atlas in Supplementary Figure S2 and all results in Supplementary Table S1. The examples include the mitochondrial/inflammatory response combination, comparator kinase mechanisms, literature-context compounds, and the DFO assay control. They are not ordered by a validated aging outcome.

Binimetinib had a positive competitive mitochondrial response ({b.loc['PROG_OXPHOS','plate_competitive_log2fc']:+.3f} log2 units; 95% CI {b.loc['PROG_OXPHOS','plate_ci95_low']:+.3f} to {b.loc['PROG_OXPHOS','plate_ci95_high']:+.3f}; q = {number(b.loc['PROG_OXPHOS','q_BH_global_plate'])}) and a negative inflammatory-annotation response ({b.loc['PROG_INFLAMMATORY','plate_competitive_log2fc']:+.3f}; 95% CI {b.loc['PROG_INFLAMMATORY','plate_ci95_low']:+.3f} to {b.loc['PROG_INFLAMMATORY','plate_ci95_high']:+.3f}; q = {number(b.loc['PROG_INFLAMMATORY','q_BH_global_plate'])}). Baicalein showed the same response directions, with mitochondrial and inflammatory q values of {number(ba.loc['PROG_OXPHOS','q_BH_global_plate'])} and {number(ba.loc['PROG_INFLAMMATORY','q_BH_global_plate'])}, respectively. Other kinase inhibitors produced different combinations, so mitochondrial protein elevation alone did not define this example. Binimetinib and baicalein occupy adjacent fixed wells H3 and H2 on all six Red plates (Supplementary Figure S4). Their related response is therefore also an adjacent-position observation; these data do not independently vary compound and coordinate.

The response was mixed across the remaining programs. Binimetinib reduced several replication-related annotations, including telomere maintenance, while increasing both aging-up and aging-down RNA-derived priors. Those priors therefore do not provide a simple, consistent reversal-of-aging interpretation. DFO remains an iron-chelation/hypoxia-mimetic assay-response control, rather than an established positive or negative control for geroprotection.

![Figure 2](figures/fig2_highlighted_heatmap.png)

**Figure 2. Selected program responses and the exploratory combined response index.** (A) Heatmap cells show primary matched-plate competitive t statistics. Asterisks denote BH q < 0.05 (*) or q < 0.01 (**) across all {ntests:,} primary tests. Color is clipped at +/-10 for readability; full statistics and intervals are in Supplementary Table S1. (B) Black circles show the equal-weight mean of 14 signed, SD-scaled primary effects. Open blue diamonds show the separate E2F/G2M regression residual. These are descriptive values without score-level p or q values. Grey segments connect the two summaries and are not confidence intervals. DFO is excluded from scaling, regression fitting and ranks. Displayed compounds are retrospective examples in the stated order; target directions are hypotheses, not universal benefits.

{index_results}

{vignette}

### Screen-wide distributions and sensitivity to growth and preprocessing

NES distributions showed different response tendencies across programs rather than a single common drug effect (Figure 3). Positive and negative scores summarize enrichment within each compound's eligible ranked protein response and do not measure activity or years of age. Binimetinib and baicalein are identified within the complete distributions rather than compared only with favorable selected examples.

DNA-repair and telomere-maintenance responses covaried with E2F response across drugs (secondary competitive-statistic r = {growth['PROG_DNA_REPAIR']:.3f} and {growth['PROG_TELOMERE_MAINT']:.3f}, respectively; Supplementary Figure S3). These associations support examining proliferation-related context, but cannot fully explain or causally remove it. Residual sensitivity and the complete protein panels retain contradictory observations. The ribosomal comparison is compatible with mitochondrial remodeling beyond a simple ribosome-depletion explanation, while retaining a relative-composition interpretation.

The fixed-correlation secondary test identified {cam_sig:,} significant pairs and its empirical-correlation sensitivity identified {emp_sig:,}, compared with {n_sig:,} under the matched-plate analysis. Such differences reflect weighting, protein correlation, missing observations, and statistical assumptions; they are not a basis for selecting whichever method gives the most significant candidates.

The vehicle-control experiment revealed {diagnostic['vehicle_split_significant_proteins']} significant protein differences and {diagnostic['vehicle_split_significant_programs']} primary plate-level program differences between the first and second listed vehicle wells at BH q < 0.05. These position-associated differences limit calibration claims and motivate the reference sensitivity analysis. {sensitivity_text} The complete reference-specific estimates and reduced plate counts are supplied. No claim of strict type I error control follows from this small control experiment.

![Figure 3](figures/fig3_nes_distributions.png)

**Figure 3. Screen-wide distributions of normalized enrichment scores.** Boxplots show the median and interquartile range of descriptive NES values for each program, with whiskers extending to 1.5 times the interquartile range. Purple diamonds identify primary matched-plate effects with global BH q < 0.05; other or untested effects are grey circles. Binimetinib is green and baicalein orange, retaining diamond versus circle shape for the same q-value distinction. Point positions show descriptive GSEA NES, while shape and significance color refer to the primary competitive-mean test. These are different summaries and significance is not a test of the plotted NES. The distributions include 160 library entries and DFO where estimable. No GSEA permutation-significance thresholds or anti-aging directions are assigned.

### Protein observations across the kinase inhibitor panel

Observed values for six mitochondrial proteins and six annotation/signaling sentinels are displayed for DMSO and six kinase inhibitors (Figure 4). The panel spans MEK, BTK, and EGFR/HER2 mechanisms. It is not a within-MEK-class experiment. Matched-plate vehicle differences were pooled and then centered and scaled to unit sample variance per protein across displayed observations. Individual values and plate identifiers are available in the figure source table; complete module and sentinel estimates for all displayed treatments are also supplied.

Binimetinib increased IDH2 and ATP5F1A in the additive protein model (log2FC {prot.loc[('Binimetinib','IDH2'),'log2fc']:+.3f} and {prot.loc[('Binimetinib','ATP5F1A'),'log2fc']:+.3f}) and reduced ICAM1 and LDLR ({prot.loc[('Binimetinib','ICAM1'),'log2fc']:+.3f} and {prot.loc[('Binimetinib','LDLR'),'log2fc']:+.3f}). The pattern was not uniform: COX4I1, NFKB1, RELA, and STAT3 did not all shift in the same direction or meet the same significance threshold. Reduced LDLR is also not intrinsically a favorable hepatic phenotype. Thus these sentinels support a specific protein-response description, rather than demonstrating generalized respiratory activation or anti-inflammatory benefit.

{direct_text} Direct comparator results are a separate testing family and include explicitly unavailable matched-plate comparisons where drugs do not share plates. The complete contrasts are provided alongside the additive-model estimates rather than inferred from differences between separate treatment-versus-vehicle significance calls.

![Figure 4](figures/fig4_kinase_inhibitor_panel.png)

**Figure 4. Scaled protein observations across kinase mechanisms.** (A) Mitochondrial proteins. (B) Inflammatory-annotation and signaling proteins. Within each plate, the DMSO mean log2 abundance was subtracted before centering and unit-variance scaling across displayed observations for each protein. Points are individual observed wells; boxes show median and interquartile range, with whiskers to 1.5 times that range. Treatment groups contain five or six plate-replicated wells and DMSO contains 23 wells across 12 plates; feature-specific missing values are omitted. Scaling is for visualization and does not supply the inferential test. STAT3 is shown as a signaling sentinel and is not a member of the supplied inflammatory annotation.

## Discussion

This reanalysis uses an existing perturbation proteome to quantify coverage and distinguish drug-induced responses across selected aging-associated protein programs. The main contribution is the comparative protein-response map and its transparent inference and sensitivity analyses. The experimental platform, dataset scale, and general mechanism-discovery workflow were established in the parent study. The present analysis does not newly demonstrate those platform capabilities or independently validate prior hallmark-based drug predictions.

The combined mitochondrial and inflammatory-annotation response in binimetinib is a testable biological observation, with baicalein providing a related example. It coexists with proliferation-associated changes and mixed aging-prior responses. A single MEK inhibitor, a single 10-micromolar concentration, and a 24-hour exposure in transformed HepG2 cells cannot establish a class mechanism, clinical exposure response, slowed aging, or rejuvenation. Off-target effects and growth-related changes remain possible. Direct cell counts, viability, protein mass per cell, and independent culture-day records are unavailable. Whole-cell abundance does not measure respiration, autophagic flux, telomere state, secreted SASP, or organismal physiology.

The pathway methods provide different views of the data. Plate-level mean responses retain the plate as the replicate unit and avoid treating correlated pathway proteins as biological replicates, but have only five or six library treatment plates and depend on parametric and missingness assumptions. The observed vehicle-position differences require caution. All 160 primary library entries have one fixed coordinate across replicate plates. Compound identity and position are therefore confounded, and plate replication cannot resolve a systematic coordinate effect. The layout diagnostic does not prove that the response of any particular compound is an artifact; it identifies a limitation that cannot be removed by row/column residualization of this dataset. Secondary competitive tests depend strongly on their correlation model; NES is a descriptive enrichment summary. None of these tests converts a protein response into a universal favorable direction for aging.

An independent signed human liver proteome reference would more directly test the aging relevance of the observed responses. A human tissue proteome atlas across five decades provides a candidate resource for such a comparison (Ding et al., 2025), subject to appropriate donor metadata, tissue-composition checks, and protein coverage. This comparison has not been performed here. It should be distinguished from reuse of RNA priors already included in the current panel. The first experimental follow-up should reposition and randomize selected treatments, including binimetinib and baicalein, with vehicle wells distributed across rows and columns. Replication in primary hepatocytes or aged-donor cells could then assess dose response, viability, cell-normalized respiration, and compartment-appropriate inflammatory outputs. Until such evidence exists, the supported conclusion is that deep proteomics measures distinguishable protein responses in these selected cellular programs and can guide specific follow-up experiments.

## Data and code availability

The deposited intensity matrix, metadata, frozen gene sets, rebuilt source tables, and reproducible analysis code are available at https://github.com/xomicsdatascience/proteomic-aging-hallmarks. The analysis configuration states filtering order, observation eligibility, contrasts, testing families, software versions, and input hashes. Figures and prose are generated from the final tables. Supplementary Table S1 contains all primary pathway comparisons; the complete protein model table includes counts and intervals. A separate research note describes potential external liver-proteome and imaging validation designs as future work, without presenting them as results.

## References

- Basisty, N., et al. (2020). A proteomic atlas of senescence-associated secretomes for aging biomarker development. PLoS Biology 18, e3000599. https://doi.org/10.1371/journal.pbio.3000599
- de Magalhaes, J. P., Curado, J., and Church, G. M. (2009). Meta-analysis of age-related gene expression profiles identifies common signatures of aging. Bioinformatics 25, 875-881. https://doi.org/10.1093/bioinformatics/btp073
- Ding, Y., et al. (2025). Comprehensive human proteome profiles across a 50-year lifespan reveal aging trajectories and signatures. Cell 188, 5763-5784.e26. https://doi.org/10.1016/j.cell.2025.06.047
- Gross, B., et al. (2026). Network-driven discovery of repurposable drugs targeting hallmarks of aging. Nature Aging 6, 1516-1531. https://doi.org/10.1038/s43587-026-01161-8
- Jiang, Y., et al. (2026). Robotic perturbation proteomics and AI agents enable scalable drug mechanism discovery. bioRxiv. https://doi.org/10.64898/2026.05.04.722718
- Lopez-Otin, C., Blasco, M. A., Partridge, L., Serrano, M., and Kroemer, G. (2023). Hallmarks of aging: An expanding universe. Cell 186, 243-278. https://doi.org/10.1016/j.cell.2022.11.001
- Mitchell, W., et al. (2024). Multi-omics characterization of partial chemical reprogramming reveals evidence of cell rejuvenation. eLife 12, RP90579. https://doi.org/10.7554/eLife.90579
- Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3, Article 3. https://doi.org/10.2202/1544-6115.1027
- Subramanian, A., et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences 102, 15545-15550. https://doi.org/10.1073/pnas.0506580102
- Tsitsipatis, D., et al. (2022). Proteomes of primary skin fibroblasts from healthy individuals reveal altered cell responses across the life span. Aging Cell 21, e13609. https://doi.org/10.1111/acel.13609
- Waldera-Lupa, D. M., et al. (2014). Proteome-wide analysis reveals an age-associated cellular phenotype of in situ aged human fibroblasts. Aging 6, 856-878. https://doi.org/10.18632/aging.100698
- Wu, D., and Smyth, G. K. (2012). Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Research 40, e133. https://doi.org/10.1093/nar/gks461
- Yoo, S., et al. (2026). Single-Cell Spatial Proteomics Uncovers Molecular Interconnectivity among Hallmarks of Aging. bioRxiv. https://doi.org/10.64898/2026.02.26.708335

## Supplementary figures and tables

**Supplementary Figure S1. Design, technical precision, replication, and vehicle diagnostics.** (A) Observation counts by role on the 12 primary plates. (B) Median per-protein log2 standard deviation among analytical pooled QC injections on each plate (sample SD). (C) Leave-one-plate-out correlations of within-plate protein fold changes for selected examples; segments show the observed range and points the mean. (D) Primary plate-level competitive statistics for the first-versus-second listed DMSO wells on the 11 plates containing both. Plot range and hit counts derive from the actual results. These diagnostics show vehicle-position differences and are not proof of calibrated whole-screen error control.

**Supplementary Figure S2. Complete primary response atlas.** All 160 library entries and the DFO assay control, ordered alphabetically, across the 14 programs. Color and asterisks use the same definitions as Figure 2. A full vector atlas and {len(range(0, primary.compound_id.nunique(), 41))} readable continuation pages are provided. Unavailable tests are left blank.

**Supplementary Figure S3. Annotation context and sensitivity analyses.** (A) The two RNA-derived aging-prior scores are components of the primary panel, not independent validation. (B) Growth-proxy versus inflammatory response, colored by mitochondrial response. (C) Comparison of secondary competitive statistics under the primary and stricter missingness filters. (D) Mean mitochondrial versus cytosolic ribosomal protein log2 fold changes. These panels describe sensitivity and relative composition rather than an aging outcome or causal adjustment.

**Supplementary Figure S4. Fixed-coordinate response context.** Competitive mean log2 fold changes for oxidative phosphorylation and inflammatory response plotted at each library entry's well coordinate in the Blue and Red assay groups. Colors are clipped at +/-0.5 and blank positions contain no displayed library entry. Bini marks binimetinib (Red H3) and Baic marks baicalein (Red H2). All library entries retain a fixed position across replicate plates, so this display cannot separate drug and coordinate effects.

**Supplementary Figure S5. Component effects for the four-axis shortlist.** The five highest four-axis indices among library entries with all four primary point estimates in the requested directions, shown in that rank order. Signs orient all four panels so positive values follow the requested profile, including reversal of the senescence effect and interval. Points show primary competitive mean log2 fold changes; segments show unadjusted 95% t intervals across effective plates. Purple diamonds denote global BH q < 0.05 within the existing {ntests:,}-test family; grey circles denote other or unavailable primary results. An interval excluding zero need not meet the global q threshold. No library entry met q < 0.05 on all four axes in the requested directions. The top three entries occupy Blue column 12, so fixed-position confounding remains a direct concern. The ranking has no composite significance test and does not establish improved pathway function.

**Supplementary Table S1.** All {ntests:,} primary compound-program comparisons, including effective plate counts, competitive mean fold changes, intervals, primary q values, secondary correlation sensitivity, descriptive NES, and eligible feature counts. Direct kinase contrasts, complete model protein estimates, every displayed observation, plate-level observations, vehicle-reference comparisons, and missingness sensitivity are supplied as separate source tables.

**Supplementary Table S2.** The 20 primary axis effects for the five displayed four-axis candidates, including 95% intervals, existing primary global q values, effective plate counts, target signs, index values, sign-filtered ranks, weakest signed components and fixed well positions. Full-screen primary effects remain in Supplementary Table S1, and complete exploratory rankings and E2F/G2M sensitivity are in Supplementary Table S3. The {len(eligible)} sign-eligible library entries are not a set of jointly significant hits.

**Supplementary Table S3.** Descriptive 14-program and four-axis indices for all library entries and DFO, with separate growth-proxy residuals. DFO has no library rank. Exact directions, equal weights and SD scales are supplied in the A06 definition table; underlying primary component effects are in Supplementary Table S1. No index-level p or q values are computed.
"""
    if '\u2014' in text:raise ValueError('Em dash found')
    if '$' in text:raise ValueError('Unrendered math marker found')
    path=ROOT/'reanalysis/A08/Aging_Proteomics_Manuscript_Revised.md';path.write_text(text,encoding='utf-8')
    ledger=[{'claim':'Primary modeled observations','value':config['n_primary_samples'],'source':'A01/sample_manifest.tsv; A02/analysis_config.json'},
            {'claim':'Retained protein features','value':inputs['retained_features_50pct'],'source':'A01/protein_mapping.tsv'},
            {'claim':'Primary pathway hypotheses','value':ntests,'source':'A02/program_effects.tsv'},
            {'claim':'Primary significant pairs q<0.05','value':n_sig,'source':'A02/program_effects.tsv:q_BH_global_plate'},
            {'claim':'Vehicle split significant proteins','value':diagnostic['vehicle_split_significant_proteins'],'source':'A04/diagnostic_summary.json'},
            {'claim':'Vehicle split significant plate programs','value':diagnostic['vehicle_split_significant_programs'],'source':'A04/null_calibration.tsv'},
            {'claim':'Binimetinib directions stable across paired vehicle references','value':str(robust),'source':'A04/vehicle_reference_sensitivity.tsv'}]
    ledger.extend([{'claim':'Four-axis point-estimate sign-eligible library entries','value':len(eligible),'source':'A06/exploratory_prioritization.tsv:four_axis_eligible_rank'},
                   {'claim':'Library entries with four requested-direction global q<0.05 effects','value':all_four_q,'source':'A02/program_effects.tsv; A06/prioritization_definitions.tsv'},
                   {'claim':'First four-axis sign-filtered candidate','value':first,'source':'A06/exploratory_prioritization.tsv:four_axis_eligible_rank'}])
    for drug in ['Binimetinib','Baicalein']:
        for pid in ['PROG_OXPHOS','PROG_INFLAMMATORY']:
            r=key.loc[(drug,pid)];ledger.append({'claim':drug+' '+pid+' competitive mean log2FC','value':r.plate_competitive_log2fc,'source':'A02/program_effects.tsv'})
    write_table(pd.DataFrame(ledger),ROOT/'reanalysis/A08/claim_result_ledger.tsv')
    (ROOT/'reanalysis/A08/results_summary.md').write_text('# Current results\n\n'+f'Primary plate-level analysis: {n_sig}/{ntests} pairs at q < 0.05 across {n_drugs} contrasts.\n\n'+'\n'.join('- '+c for c in clauses)+'\n\n'+sensitivity_text+'\n\n'+index_results+'\n\n'+vignette+'\nNo index-level significance, age reduction, clinical survival, or independent aging validation is claimed.\n')
    provenance={'input_summary':inputs,'analysis_config':config,'abstract_words':len(abstract.split()),'source_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((ROOT/'reanalysis').glob('*.py'))},'claim_ledger':'reanalysis/A08/claim_result_ledger.tsv'}
    json_write(provenance,ROOT/'reanalysis/A08/build_provenance.json')
    print(f'Manuscript generated: {len(text.split())} words; abstract {len(abstract.split())} words',flush=True)

if __name__=='__main__':main()
