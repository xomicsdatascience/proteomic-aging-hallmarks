# Audit resolution and remaining limits

The revision removes composite geroprotection scores, inferred universally favorable program directions, clinical drug rankings, and the unexecuted Medicare section. It regenerates provenance, statistics, figures, source tables, and manuscript text from one workflow.

The filtering order is now explicit: exclude the unlabeled well before the feature missingness filter. This retains 7,871 protein features. The primary design has 1,004 observations and rank 173. Compound aliases are normalized once and preserved in an original-to-canonical map. Protein-level tests require at least three observed values in each contrasted group.

Primary pathway inference uses plate-level competitive mean log2 fold changes. Proteins within a plate are not independent biological replicates. Fixed-rho CAMERA-PR-style and pairwise empirical-correlation tests remain separately labeled sensitivity analyses. GSEA NES is descriptive; its native gene-set permutation p values are not presented as confirmatory significance. The three analyses do not have identical null hypotheses or weighting.

The earlier audit's statement that raw DMSO controls produced no significant proteins was incorrect. Missing p values were not handled correctly in that audit calculation. With explicit missing-test handling, the former 7,872-feature vehicle split gives 53 significant proteins; the corrected 7,871-feature split gives 52 at BH q < 0.05. Excluding low-observation extreme t statistics also exposes a non-null pathway diagnostic that was masked in the old competitive analysis. These controls are now reported without a claim of strict calibration. The additional corruption caused by repeatedly rescaling the saved null table was a separate problem and has been removed.

Figure generation is read-only with respect to A01-A07 outputs and checks their hashes. Figure axes and hit labels derive from the actual control data. Quantities in the manuscript derive from the final tables. Complete displayed-protein observations, matched-plate pathway observations, direct kinase comparisons, and missingness sensitivity are exported.

All primary library entries retain one coordinate across replicate plates; binimetinib and baicalein are adjacent at Red H3 and H2. A new layout diagnostic makes this compound-position confounding explicit. It is not corrected by regression or by plate replication, and randomized repositioned validation is the first experimental priority.

The clinical imaging and external liver-proteome options in VALIDATION_OPTIONS.md are research proposals, not completed validation. The current conclusion is limited to measurable protein-program responses. No functional rejuvenation, organ age reduction, longevity, or clinical treatment recommendation is claimed.

## Verification of the revision

All inference stages were rebuilt from the deposited matrix. The final figure/manuscript workflow and output validator pass for 2,254 primary pathway comparisons and 1,267,231 protein-model rows. Seven regression tests pass, including missing-p-value handling, observed-data estimability, finite test backgrounds, and preservation of the plate replicate unit. A clean-checkout recovery check reconstructs the full protein table exactly from its 41 committed gzip partitions. Repeated figure generation preserves all A01-A07 hashes, including nested protein partitions.

The four compatibility Word exports are byte-identical. The 11-page canonical manuscript and 7-page supplementary figure document were rendered and visually reviewed, including the four full-atlas continuation pages. These checks verify the implementation and exports; they do not supply independent biological validation or resolve the fixed-position design limitation.
