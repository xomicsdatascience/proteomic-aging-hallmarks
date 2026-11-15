"""Exploratory response examples and complete panel source data."""
import numpy as np
import pandas as pd
from analysis_common import (ROOT, PRIMARY_IDS, HIGHLIGHTS, KINASE_PANEL, SENTINELS,
    read_table, write_table, load_gmt)
p=read_table('reanalysis/A02/program_effects.tsv');protein=read_table('reanalysis/A02/protein_effects.tsv.gz')
spec=read_table('reanalysis/A06/candidate_specificity.tsv');metadata=read_table('reanalysis/A01/compound_manifest.tsv')
primary=p[p.program_id.isin(PRIMARY_IDS)]
counts=primary.groupby('compound_id').q_BH_global_camera.apply(lambda x:int((x<.05).sum())).rename('n_programs_q05')
evidence=spec.merge(counts,on='compound_id').merge(metadata[['compound_id','Category','Target','CAS_Number']],on='compound_id',how='left',validate='one_to_one')
evidence['displayed_example']=evidence.compound_id.isin(HIGHLIGHTS)
evidence['selection_scope']='Retrospective examples for response comparison; no validated geroprotector ranking'
write_table(evidence.sort_values('compound_id'),ROOT/'reanalysis/A07/candidate_evidence.tsv')
gmt=load_gmt();genes=gmt['PROG_OXPHOS']|gmt['PROG_INFLAMMATORY']|set(sum(SENTINELS.values(),[]))
panel=list(dict.fromkeys(KINASE_PANEL+['Baicalein','DFO','Nordihydroguaiaretic acid','Lonafarnib']))
module=protein[protein.compound_id.isin(panel)&protein.feature_id.isin(genes)].copy()
module['in_oxphos_annotation']=module.feature_id.isin(gmt['PROG_OXPHOS'])
module['in_inflammatory_annotation']=module.feature_id.isin(gmt['PROG_INFLAMMATORY'])
module['plotted_sentinel']=module.feature_id.isin(sum(SENTINELS.values(),[]))
write_table(module,ROOT/'reanalysis/A07/candidate_protein_effects.tsv')
# Preserve missing observations: pairwise-complete protein response similarities.
pt=protein.pivot(index='feature_id',columns='compound_id',values='t_mod')
pg=primary.pivot(index='program_id',columns='compound_id',values='camera_t')
cp=pt.corr(min_periods=100);cg=pg.corr(min_periods=10).reindex(index=cp.index,columns=cp.columns)
iu=np.triu_indices(len(cp),1);a,b=cp.to_numpy()[iu],cg.to_numpy()[iu];valid=np.isfinite(a)&np.isfinite(b)
r=float(np.corrcoef(a[valid],b[valid])[0,1])
write_table(pd.DataFrame([{'comparison':'descriptive correlation of pairwise drug similarities in protein and program spaces','value':r,'n_pairs':int(valid.sum()),'inference':'pairs share compounds; no independent-pair p value'}]),ROOT/'reanalysis/A07/drug_similarity.tsv')
(ROOT/'reanalysis/A07/candidate_selection_rules.yaml').write_text("""scope: retrospective illustrative comparison, not prospective candidate validation\nprimary_pathway_threshold: fixed-correlation CAMERA-PR-style global BH q < 0.05\ncompound_universe: 160 library entries plus DFO assay response control\nprimary_examples: [Binimetinib, Baicalein]\nselection_rationale: response combination of mitochondrial and inflammatory annotations, alongside growth sensitivity and contradictory programs\ncomposite_geroprotection_score: not used\nfunctional_or_organ_age_validation: not performed\n""")
(ROOT/'reanalysis/A07/candidate_report.md').write_text("""# Interpretation of response examples\n\nBinimetinib and baicalein are retrospective examples of a mitochondrial/inflammatory protein-response combination. They are not predetermined winners of a validated aging screen. Complete mitochondrial and inflammatory annotation members, every displayed sentinel, contradictory results, observation counts, and model intervals are available for all displayed kinase inhibitors. A single MEK inhibitor does not establish a MEK-class effect. Ten micromolar, 24-hour HepG2 responses do not establish clinical exposure or on-target specificity.\n\nDFO remains an assay response control. No compound is relabeled an established geroprotector, pro-aging control, or cytotoxic control from these data alone. A second MEK inhibitor or genetic perturbation, primary-cell dose/viability work, cell-normalized respiration, and compartment-appropriate inflammatory assays would test the proposed mechanism.\n""")
