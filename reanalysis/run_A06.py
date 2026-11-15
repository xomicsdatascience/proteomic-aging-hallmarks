"""Descriptive growth/stress sensitivity, with no composite geroprotection score."""
import numpy as np
import pandas as pd
from scipy import stats
from analysis_common import ROOT, PRIMARY_IDS, read_table, write_table, json_write, load_gmt
p=read_table('reanalysis/A02/program_effects.tsv')
protein=read_table('reanalysis/A02/protein_effects.tsv.gz')
gmt=load_gmt();pivot=p.pivot(index='compound_id',columns='program_id',values='camera_t')
eligible=protein[protein.eligible_for_testing]
shared=eligible.groupby('compound_id').agg(mean_abs_t=('t_mod',lambda s:np.nanmean(np.abs(s))),mean_abs_log2fc=('log2fc',lambda s:np.nanmean(np.abs(s))),n_de_q05=('q_BH_contrast',lambda s:int((s<.05).sum()))).join(pivot).reset_index()
write_table(shared,ROOT/'reanalysis/A06/shared_response_scores.tsv')
specific=shared.copy();X=np.column_stack([np.ones(len(shared)),shared.NUIS_E2F_PROLIF,shared.mean_abs_t])
for target in ['PROG_OXPHOS','PROG_INFLAMMATORY','PROG_MTORC1_SIGNALING','PROG_MACROAUTOPHAGY']:
    y=shared[target].to_numpy();ok=np.isfinite(y)&np.isfinite(X).all(axis=1)
    z=np.full(len(y),np.nan);resid=y[ok]-X[ok]@np.linalg.lstsq(X[ok],y[ok],rcond=None)[0]
    z[ok]=(resid-resid.mean())/resid.std(ddof=1);specific[target+'_resid_z']=z
specific['interpretation']='Screen-relative sensitivity residual; not causal adjustment or geroprotection'
write_table(specific,ROOT/'reanalysis/A06/candidate_specificity.tsv')
composition=[]
for drug,group in eligible.groupby('compound_id'):
    ox=group[group.feature_id.isin(gmt['PROG_OXPHOS'])].log2fc.mean()
    ribo=group[group.feature_id.str.startswith(('RPL','RPS'))].log2fc.mean()
    composition.append({'compound_id':drug,'oxphos_mean_log2fc':ox,'cytosolic_ribo_mean_log2fc':ribo,'oxphos_minus_ribo_log2fc':ox-ribo})
write_table(pd.DataFrame(composition),ROOT/'reanalysis/A06/composition_sensitivity.tsv')
r={k:float(pivot[k].corr(pivot.NUIS_E2F_PROLIF)) for k in ['PROG_OXPHOS','PROG_INFLAMMATORY','PROG_DNA_REPAIR','PROG_TELOMERE_MAINT','NUIS_G2M_PROLIF']}
json_write(r,ROOT/'reanalysis/A06/growth_correlations.json')
(ROOT/'reanalysis/A06/specificity_report.md').write_text("""# Growth and composition sensitivity\n\nResidual Z scores describe deviations from relationships fitted across this compound collection. E2F response and overall perturbation magnitude are post-treatment molecular proxies, not measured baseline confounders. Residualization can remove real drug biology and does not demonstrate an effect independent of growth. No beneficial target directions, composite geroprotection scores, or composite significance tests are computed.\n\nThe ribosome comparison examines one compositional explanation. It cannot establish absolute protein or organelle abundance per cell. Direct cell counts, viability, and total protein per cell are unavailable.\n""")

# Fixed-coordinate context is descriptive; drug and position effects are not identifiable separately.
samples=read_table('reanalysis/A01/sample_manifest.tsv')
library=samples[samples.condition_role.eq('library_compound')&~samples.excluded&samples.plate.ne('RedMix')].copy()
layout=library.groupby('compound_id').agg(n_unique_wells=('well','nunique'),well=('well','first'),assay_group=('plate',lambda x:','.join(sorted(set(str(v).rstrip('0123456789') for v in x))))).reset_index()
layout['well_row']=layout.well.str[0];layout['well_column']=layout.well.str[1:].astype(int)
layout['interpretation']='Fixed compound-position association; position is confounded with compound identity'
layout=layout.merge(p[p.program_id.isin(PRIMARY_IDS)][['compound_id','program_id','plate_competitive_log2fc','q_BH_global_plate']],on='compound_id',validate='one_to_many')
write_table(layout,ROOT/'reanalysis/A06/layout_context.tsv')
(ROOT/'reanalysis/A06/layout_report.md').write_text("""# Fixed-position context\n\nAll 160 primary library entries occupy one fixed well coordinate across their replicate plates. Binimetinib occupies Red-plate H3 and baicalein Red-plate H2, so their related response is also an adjacent-position observation. Vehicle wells occupy E1 and F1. Repeated plates do not independently vary compound and position. The vehicle differences motivate concern about position effects, but do not establish that a particular drug response is an artifact. These data cannot identify separate compound and well effects. No row/column residualization is presented as a correction. Randomized or repositioned treatment replication with distributed vehicle controls is required to resolve this limitation.\n""")
