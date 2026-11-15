# Design and provenance

The deposited matrix contains 1,232 observations and 8,630 features. One unlabeled well is excluded before filtering. The <=50% zero-value filter retains 7,871 features; the <=20% filter retains 7,296. This order excludes the boundary feature ZNF552 that was retained when filtering before the unlabeled-well exclusion.

The primary analysis compares 160 library entries and DFO with within-plate DMSO controls on 12 plates. RedMix is excluded from primary inference. Culture and acquisition metadata not in the deposited file remain marked unavailable. Plate/well identifiers do not verify independent culture days. The parent study calls these biological replicates.

Input hashes and derived counts are in input_summary.json. Compound aliases are resolved once in compound_identifier_mapping.tsv.
