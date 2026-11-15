"""Export exploratory indices separately from inferential pathway tables."""
from analysis_common import ROOT,read_table,write_table,short_drug
from prioritization import indices
p=read_table('reanalysis/A02/program_effects.tsv')
r,s,c=indices(p)
layout=read_table('reanalysis/A06/layout_context.tsv')[['compound_id','well','assay_group']].drop_duplicates()
r=r.merge(layout,on='compound_id',how='left',validate='one_to_one')
write_table(r,ROOT/'reanalysis/A06/exploratory_prioritization.tsv')
write_table(s,ROOT/'reanalysis/A06/prioritization_definitions.tsv')
write_table(c,ROOT/'reanalysis/_cache/prioritization_components.tsv')
eligible=r[r.four_axis_eligible_rank.notna()].sort_values('four_axis_eligible_rank')
top=', '.join(short_drug(v) for v in eligible.head(5).compound_id)
(ROOT/'reanalysis/A06/prioritization_report.md').write_text(f'''# Exploratory response prioritization

The 14-program combined response index retains the earlier score's explicit target-direction hypotheses. It averages each signed primary matched-plate competitive mean log2 fold change divided by that program's sample SD across the 160 library entries. It does not subtract the screen mean. The four-axis index uses only positive OXPHOS, macroautophagy and proteasome responses and negative senescence-associated machinery responses. All directions, equal weights, scales and available library counts are exported in `prioritization_definitions.tsv`. Complete underlying primary effects, 95% intervals and global q values remain in the existing A02 table and Supplementary Table S1. The script regenerates a detailed signed-component table into the ignored local `_cache` directory; this redundant export is not published.

These indices are descriptive prioritization rules, not geroprotection measurements or normally distributed Z statistics. There are no index-level p values, q values, significance cutoffs or null-covariance estimates. Cross-program overlap is not used as a null correlation matrix. Program overlap and cross-program correlation can still give related biology repeated influence in the mean. The legacy 14-program signs are exploratory choices, not universally beneficial directions. The earlier combined raw/residual formula is retired. Each raw index and its E2F/G2M regression residual are reported separately; neither their average nor an independence claim is made. Residualization may remove real treatment biology. DFO is excluded from all scale estimates, regression fits and ranks. Missing components leave an index unavailable rather than reducing its denominator.

For the requested four-axis profile, {len(eligible)} library entries have primary point estimates in all four specified directions. This is a sign filter, not a joint significance criterion. The filtered ranking starts with {top}. The unrestricted four-axis rank is also exported, because an average can conceal an oppositely directed component. Weakest signed components and every per-axis interval and q value remain visible. Rank differences do not quantify confidence in an ordering. Fixed well and assay-group fields expose the unresolved compound-position confounding.

The protein abundance response does not establish respiration, autophagic flux, proteasome activity or reduced senescence. Randomized repositioned replication with viability and cell-number measurements should precede functional interpretation.
''')
