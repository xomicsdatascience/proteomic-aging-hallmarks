"""Fail on mismatched counts, invalid joins, stale score fields, or broken exports."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from docx import Document
from analysis_common import ROOT, PRIMARY_IDS, KINASE_PANEL, SENTINELS, read_table, load_raw, primary_mask, bh, sha256


def main():
    config=json.loads((ROOT/'reanalysis/A02/analysis_config.json').read_text())
    sm=read_table('reanalysis/A01/sample_manifest.tsv');mapping=read_table('reanalysis/A01/protein_mapping.tsv')
    X,obs,features=load_raw();keep=(X[obs.drug.ne('Drug')]==0).mean(axis=0)<=.5
    assert np.array_equal(keep,mapping.retained_50pct_filter)
    assert int(keep.sum())==config['n_features']
    assert int(primary_mask(sm).sum())==config['n_primary_samples']
    p=read_table('reanalysis/A02/program_effects.tsv');primary=p[p.program_id.isin(PRIMARY_IDS)]
    assert len(primary)==config['n_primary_compounds']*len(PRIMARY_IDS)==config['n_primary_program_tests']
    assert not p.duplicated(['compound_id','program_id']).any()
    assert not any(c.startswith('CGS') or c.startswith('n_beneficial') for c in p)
    for pc,qc in [('p_plate','q_BH_global_plate'),('p_camera','q_BH_global_camera'),('p_camera_emp','q_BH_global_camera_emp')]:
        np.testing.assert_allclose(bh(primary[pc]),primary[qc],rtol=2e-8,atol=1e-10,equal_nan=True)
    assert (primary.loc[primary.p_plate.notna(),'n_effective_plates']>=3).all()
    assert primary.loc[primary.n_effective_plates.lt(3),'p_plate'].isna().all()
    protein=read_table('reanalysis/A02/protein_effects.tsv.gz')
    test=protein.eligible_for_testing
    assert (protein.loc[test,'n_treated_obs']>=3).all()
    assert (protein.loc[test,'n_reference_obs']>=3).all()
    assert protein.loc[~test,'q_BH_contrast'].isna().all()
    partition_dir=ROOT/'reanalysis/A02/protein_effects';manifest=pd.read_csv(partition_dir/'manifest.tsv',sep='\t')
    assert manifest.compound_id.nunique()==config['n_primary_compounds']
    assert int(manifest.n_rows.sum())==len(protein)
    for name,group in manifest.groupby('file_name'):
        assert group.sha256.nunique()==1 and sha256(partition_dir/name)==group.sha256.iloc[0]
        part=pd.read_csv(partition_dir/name,sep='\t')
        expected=group.set_index('compound_id').n_rows.sort_index()
        pd.testing.assert_series_equal(part.groupby('compound_id').size().sort_index(),expected,check_names=False)
    assert not protein.duplicated(['compound_id','feature_id']).any()
    source=read_table('reanalysis/A08/source_tables/Supplementary_Table_S1_all_compounds.tsv')
    cols=['compound_id','program_id','plate_t','q_BH_global_plate','camera_t','q_BH_global_camera']
    a=primary[cols].sort_values(cols[:2]).reset_index(drop=True);b=source[cols].sort_values(cols[:2]).reset_index(drop=True)
    pd.testing.assert_frame_equal(a,b,check_dtype=False,rtol=1e-8)
    plotted=read_table('reanalysis/A08/source_tables/fig4_scaled_protein_observations.tsv')
    assert not plotted.duplicated(['sample_id','feature_id']).any()
    assert set(plotted.compound_id)==set(['DMSO']+KINASE_PANEL)
    for _,g in plotted.groupby('feature_id'):
        assert abs(g.scaled_abundance.mean())<1e-8
        assert abs(g.scaled_abundance.std(ddof=1)-1)<1e-8
    module=read_table('reanalysis/A07/candidate_protein_effects.tsv')
    for drug in KINASE_PANEL:
        assert set(sum(SENTINELS.values(),[]))<=set(module.loc[module.compound_id.eq(drug),'feature_id'])
    markdown=(ROOT/'reanalysis/A08/Aging_Proteomics_Manuscript_Revised.md').read_text()
    assert '\u2014' not in markdown and '$' not in markdown and 'CGS_' not in markdown
    assert 'Medicare' not in markdown
    documents=[ROOT/'Aging_Hallmarks_Proteomics_Report.docx',ROOT/'results/Aging_Hallmarks_Proteomics_Report.docx',ROOT/'reanalysis/A08/Aging_Proteomics_Manuscript_Final.docx',ROOT/'reanalysis/A08/Aging_Proteomics_Manuscript_Revised.docx']
    assert len({sha256(f) for f in documents})==1
    doc=Document(documents[-1]);text='\n'.join(p.text for p in doc.paragraphs)
    assert '\\frac' not in text and '$' not in text and '\u2014' not in text
    assert 'q < 0.01 (**)' in text
    assert len(doc.inline_shapes)==4
    for name in ['fig1_program_coverage','fig2_highlighted_heatmap','fig3_nes_distributions','fig4_kinase_inhibitor_panel','figS1_qc_and_null_calibration','figS2_full_160_compound_atlas','figS3_benchmarks_and_sensitivity','figS4_layout_context']:
        for extension in ['png','svg']:assert (ROOT/'reanalysis/A08/figures'/(name+'.'+extension)).exists()
    summary={'status':'passed','primary_samples':config['n_primary_samples'],'retained_features':config['n_features'],'primary_program_tests':len(primary),'protein_rows':len(protein),'protein_partitions':int(manifest.file_name.nunique()),'synchronized_word_exports':len(documents),'primary_plate_pairs_q05':int((primary.q_BH_global_plate<.05).sum()),'note':'Structural and numerical checks; Word layout is reviewed separately.'}
    (ROOT/'reanalysis/A08/validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary,flush=True)

if __name__=='__main__':main()
