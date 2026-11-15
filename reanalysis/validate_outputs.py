"""Fail on mismatched counts, invalid joins, stale score fields, or broken exports."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from docx import Document
from analysis_common import ROOT, PRIMARY_IDS, KINASE_PANEL, SENTINELS, read_table, load_raw, primary_mask, bh, sha256
from prioritization import indices, COMBINED_DIRECTIONS, ASSAY_DIRECTIONS


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
    # The indices are derived, descriptive outputs; primary significance is unchanged.
    ranking=read_table('reanalysis/A06/exploratory_prioritization.tsv')
    expected,definitions,components=indices(p)
    pd.testing.assert_frame_equal(ranking[expected.columns],expected,check_dtype=False,rtol=1e-8,atol=1e-9)
    pd.testing.assert_frame_equal(read_table('reanalysis/A06/prioritization_definitions.tsv'),definitions,check_dtype=False,rtol=1e-8)
    assert set(COMBINED_DIRECTIONS)==set(PRIMARY_IDS) and len(ASSAY_DIRECTIONS)==4
    assert not any(c.startswith(('p_','q_','CGS')) for c in ranking.columns)
    assert ranking.loc[ranking.compound_id.eq('DFO')].filter(regex='_rank$').isna().all(axis=None)
    assert ranking.loc[ranking.four_axis_eligible_rank.notna(),'all_four_requested_directions'].all()
    pd.testing.assert_frame_equal(read_table('reanalysis/A08/source_tables/Supplementary_Table_S3_exploratory_indices.tsv'),ranking,check_dtype=False)
    four=read_table('reanalysis/A08/source_tables/Supplementary_Table_S2_four_axis_candidates.tsv')
    selected_ids=set(ranking[ranking.four_axis_eligible_rank.le(5)].compound_id)
    assert set(four.compound_id)==selected_ids and len(four)==20
    assert four.groupby('compound_id').size().eq(4).all()
    selected_cols=['compound_id','program_id','plate_competitive_log2fc','plate_ci95_low','plate_ci95_high','q_BH_global_plate','n_effective_plates']
    a=primary[primary.program_id.isin(ASSAY_DIRECTIONS)&primary.compound_id.isin(selected_ids)][selected_cols].sort_values(cols[:2]).reset_index(drop=True)
    b=four[selected_cols].sort_values(cols[:2]).reset_index(drop=True)
    pd.testing.assert_frame_equal(a,b,check_dtype=False,rtol=1e-8)
    assert np.array_equal(four.target_direction,four.program_id.map(ASSAY_DIRECTIONS))
    np.testing.assert_allclose(four.signed_effect,four.plate_competitive_log2fc*four.target_direction,rtol=1e-8)
    all_four=components[components.four_axis_component]
    all_four_q=int(all_four.assign(hit=all_four.signed_scaled_effect.gt(0)&all_four.q_BH_global_plate.lt(.05)).groupby('compound_id').hit.all().drop(index='DFO').sum())
    fig3=read_table('reanalysis/A08/source_tables/fig3_nes_all_compounds.tsv')
    a=primary[['compound_id','program_id','nes_gsea','q_BH_global_plate']].sort_values(cols[:2]).reset_index(drop=True)
    b=fig3[a.columns].sort_values(cols[:2]).reset_index(drop=True)
    pd.testing.assert_frame_equal(a,b,check_dtype=False,rtol=1e-8)
    hit=fig3.q_BH_global_plate.lt(.05)
    assert np.array_equal(fig3.primary_q_lt_005,hit)
    assert np.array_equal(fig3.point_marker,np.where(hit,'D','o'))
    for drug,color in [('Binimetinib','#009E73'),('Baicalein','#D55E00')]:
        assert fig3.loc[fig3.compound_id.eq(drug),'point_color'].eq(color).all()
    others=~fig3.compound_id.isin(['Binimetinib','Baicalein'])
    assert np.array_equal(fig3.loc[others,'point_color'],np.where(hit[others],'#CC79A7','#999999'))
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
    supplement=Document(ROOT/'reanalysis/A08/Aging_Proteomics_Supplementary_Figures.docx')
    assert len(supplement.inline_shapes)==7
    vignette=Document(ROOT/'reanalysis/A08/Aging_Proteomics_Supplementary_Figure_S5.docx')
    assert len(vignette.inline_shapes)==1
    assert any('Supplementary Figure S5.' in para.text for para in vignette.paragraphs)
    for name in ['fig1_program_coverage','fig2_highlighted_heatmap','fig3_nes_distributions','fig4_kinase_inhibitor_panel','figS1_qc_and_null_calibration','figS2_full_160_compound_atlas','figS3_benchmarks_and_sensitivity','figS4_layout_context','figS5_four_axis_candidates']:
        for extension in ['png','svg']:assert (ROOT/'reanalysis/A08/figures'/(name+'.'+extension)).exists()
    summary={'status':'passed','primary_samples':config['n_primary_samples'],'retained_features':config['n_features'],'primary_program_tests':len(primary),'protein_rows':len(protein),'protein_partitions':int(manifest.file_name.nunique()),'synchronized_word_exports':len(documents),'primary_plate_pairs_q05':int((primary.q_BH_global_plate<.05).sum()),'four_axis_sign_eligible_library_entries':int(ranking.four_axis_eligible_rank.notna().sum()),'four_axis_all_global_q05_library_entries':all_four_q,'published_four_axis_candidate_effects':len(four),'figure3_primary_q_styles_checked':len(fig3),'note':'Structural and numerical checks; Word layout is reviewed separately.'}
    (ROOT/'reanalysis/A08/validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary,flush=True)

if __name__=='__main__':main()
