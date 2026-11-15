"""Contextual aging annotations, explicitly separated from independent validation."""
from analysis_common import ROOT, read_table, write_table
p=read_table('reanalysis/A02/program_effects.tsv')
t=p.pivot(index='compound_id',columns='program_id',values='camera_t')
q=p.pivot(index='compound_id',columns='program_id',values='q_BH_global_camera')
context=t[['PRIOR_DEMAGALHAES_UP','PRIOR_DEMAGALHAES_DN']].copy()
context.columns=['demagalhaes_up_camera_t','demagalhaes_dn_camera_t']
context['demagalhaes_up_q_BH']=q.PRIOR_DEMAGALHAES_UP
context['demagalhaes_dn_q_BH']=q.PRIOR_DEMAGALHAES_DN
context['interpretation']='Panel components, not independent validation'
write_table(context.reset_index(),ROOT/'reanalysis/A05/drug_reference_comparisons.tsv')
manifest=[
 {'reference':'de Magalhaes 2009','doi':'10.1093/bioinformatics/btp073','compartment':'RNA-derived mammalian aging annotation','status':'included in primary panel; contextual only'},
 {'reference':'Tsitsipatis 2022','doi':'10.1111/acel.13609','compartment':'Human donor fibroblast whole-cell proteome','status':'prospective independent reference; not evaluated'},
 {'reference':'Waldera-Lupa 2014','doi':'10.18632/aging.100698','compartment':'Human donor fibroblast whole-cell proteome; PMC4247387','status':'prospective independent reference; not evaluated'},
 {'reference':'Mitchell 2024','doi':'10.7554/eLife.90579','compartment':'Chemical reprogramming proteome','status':'prospective intervention reference; not evaluated'},
 {'reference':'Ding 2025','doi':'10.1016/j.cell.2025.06.047','compartment':'Human tissue proteome aging atlas, including liver','status':'prospective liver-relevant signed reference; not evaluated'},
 {'reference':'Basisty 2020','doi':'10.1371/journal.pbio.3000599','compartment':'Secretome; different from the current whole-cell lysate','status':'contextual literature; not evaluated as matching validation'}]
write_table(__import__('pandas').DataFrame(manifest),ROOT/'reanalysis/A05/external_reference_manifest.tsv')
coverage=read_table('reanalysis/A03/program_manifest.tsv')
write_table(coverage[coverage.program_id.str.startswith('PRIOR_')][['program_id','total_genes','measured_retained_7871']],ROOT/'reanalysis/A05/reference_coverage.tsv')
(ROOT/'reanalysis/A05/aging_benchmark_report.md').write_text("""# Aging reference scope\n\nNo independent proteomic aging dataset has been analyzed in this build. The two de Magalhaes RNA-derived priors are already members of the primary panel and must not be described as independent validation. No hard-coded variance inflation is used for an additional canonical senescence pathway. The reference inventory distinguishes evaluated panel context from proposed external comparisons. A signed human liver aging proteome comparison is a proposed next analysis, not a completed result or an estimate of years rejuvenated.\n""")
