"""Observed-data protein model and CAMERA-PR-style pathway comparisons."""
from __future__ import annotations
import importlib.metadata
import platform
import time
import numpy as np
import pandas as pd
from scipy import stats
import gseapy as gp
from analysis_common import (ROOT, PRIMARY_IDS, PRIMARY_RHO, MIN_TREATED_OBS, MIN_CONTROL_OBS,
    KINASE_PANEL, bh, write_table, read_table, json_write, normalized_data, load_gmt,
    make_design, fit_observed_model, moderate_variances, competitive_test,
    residual_correlations, empirical_rho, primary_mask, plate_program_comparisons, write_protein_partitions)


def pathway_results(T, B, features, labels, gene_sets, correlations, pipeline, with_nes=False):
    rows=[]
    for ci,label in enumerate(labels):
        values=T[:,ci]
        available=np.isfinite(values)
        nes={}
        if with_nes:
            ranking=pd.Series(values[available],index=features[available]).sort_values(ascending=False,kind='stable')
            result=gp.prerank(rnk=ranking,gene_sets={k:sorted(v) for k,v in gene_sets.items()},min_size=5,max_size=500,permutation_num=2000,seed=42,threads=1,outdir=None,verbose=False)
            nes=result.results
        for pid,genes in gene_sets.items():
            member=np.isin(features,list(genes))
            tc,pc,m,vif=competitive_test(values,member,PRIMARY_RHO)
            idx,corr=correlations[pid]
            rho,pair_fraction=empirical_rho(corr,available[idx])
            te,pe,_,vife=competitive_test(values,member,rho)
            in_genes=member & available
            rec={'pipeline':pipeline,'compound_id':label,'contrast_id':label+'_vs_DMSO' if pipeline=='primary_plate_blocked_eb' else label,
                 'program_id':pid,'program_family':'nuisance_control' if pid.startswith('NUIS_') else 'primary_aging',
                 'n_genes_measured':m,'n_background_proteins':int(available.sum()),'rho_camera_fixed':PRIMARY_RHO,'vif_camera_fixed':vif,
                 'camera_t':tc,'p_camera':pc,'rho_inter_protein':rho,'vif_emp':vife,'empirical_pair_fraction':pair_fraction,
                 'camera_t_emp':te,'p_camera_emp':pe,'mean_log2fc':float(np.nanmean(B[in_genes,ci])) if m else np.nan,
                 'test_status':'estimable' if np.isfinite(pc) else 'insufficient_observed_features'}
            rec['nes_gsea']=nes.get(pid,{}).get('nes',np.nan)
            rec['lead_genes']=nes.get(pid,{}).get('lead_genes','')
            rows.append(rec)
        if with_nes and ((ci+1)%20==0 or ci+1==len(labels)):
            print(f'NES {ci+1}/{len(labels)}',flush=True)
    frame=pd.DataFrame(rows)
    for family in frame.program_family.unique():
        mask=frame.program_family.eq(family)
        for pcol,qcol in [('p_camera','q_BH_global_camera'),('p_camera_emp','q_BH_global_camera_emp')]:
            frame.loc[mask,qcol]=bh(frame.loc[mask,pcol])
    return frame


def write_inference_report(primary, config):
    p=primary[primary.program_id.isin(PRIMARY_IDS)]
    report=f"""# Protein and pathway inference

The model uses {config['n_primary_samples']:,} primary assay observations and {config['n_features']:,} retained features, with design rank {config['design_rank']}. Empirical Bayes prior: d0={config['d0_prior']:.6g}, s0 squared={config['s0_2_prior']:.6g}. Protein tests require at least three observed values in each contrasted group.

At BH q < 0.05 across {len(p):,} compound-program hypotheses, the primary matched-plate mean-response analysis identifies {int((p.q_BH_global_plate<.05).sum())} pairs. The secondary fixed-correlation CAMERA-PR-style analysis identifies {int((p.q_BH_global_camera<.05).sum())} pairs, and the pairwise empirical-correlation sensitivity identifies {int((p.q_BH_global_camera_emp<.05).sum())} pairs. These analyses have different assumptions and weights. NES is descriptive. Direct kinase comparisons have their own testing family.

Plate-response inference assumes independent plates. Fixed treatment coordinates, incomplete culture provenance, and observed vehicle-position differences limit interpretation; see A04 and A06. Neither discovery counts nor correlation sensitivity establish calibrated whole-screen error control or beneficial aging effects.
"""
    (ROOT/'reanalysis/A02/inference_report.md').write_text(report)
    return report


def main():
    start=time.monotonic()
    sm=read_table('reanalysis/A01/sample_manifest.tsv')
    eligibility=read_table('reanalysis/A01/contrast_eligibility.tsv')
    mapping=read_table('reanalysis/A01/protein_mapping.tsv')
    Y,obs,features=normalized_data()
    if not np.array_equal(obs.sample_id,sm.sample_id): raise ValueError('Sample order mismatch')
    mask=primary_mask(sm).to_numpy(); samples=sm.loc[mask].reset_index(drop=True); A=Y[mask]
    compounds=eligibility.loc[eligibility.eligibility_status.eq('directly_estimable'),'compound_id'].tolist()
    D,columns=make_design(samples,compounds)
    direct_pairs=[('Binimetinib',d) for d in KINASE_PANEL if d!='Binimetinib']
    C=np.zeros((len(columns),len(compounds)+len(direct_pairs)))
    for i,c in enumerate(compounds): C[columns.index(c),i]=1
    for i,(a,b) in enumerate(direct_pairs,start=len(compounds)):
        C[columns.index(a),i]=1; C[columns.index(b),i]=-1
    print(f'Fitting {len(samples)} samples, {len(features)} proteins, design rank {len(columns)}',flush=True)
    B,U,V,df,resid=fit_observed_model(A,D,C)
    post,d0,s0=moderate_variances(V,df)
    SE=np.sqrt(U*post[:,None]); T=B/SE
    labels=compounds+[a+'_vs_'+b for a,b in direct_pairs]
    treated_n=np.column_stack([np.isfinite(A[samples.compound_id.eq(c).to_numpy()]).sum(axis=0) for c in compounds])
    dmso_n=np.isfinite(A[samples.compound_id.eq('DMSO').to_numpy()]).sum(axis=0)
    valid=np.column_stack([(treated_n[:,i]>=MIN_TREATED_OBS)&(dmso_n>=MIN_CONTROL_OBS) for i in range(len(compounds))]+[
        (treated_n[:,compounds.index(a)]>=MIN_TREATED_OBS)&(treated_n[:,compounds.index(b)]>=MIN_TREATED_OBS) for a,b in direct_pairs])
    T=np.where(valid,T,np.nan)
    P=2*stats.t.sf(np.abs(T),df=df[:,None]+d0)
    protein_frames=[]
    for ci,label in enumerate(labels):
        if ci<len(compounds): nt=treated_n[:,ci]; nr=dmso_n; ref='DMSO'
        else:
            a,b=direct_pairs[ci-len(compounds)];nt=treated_n[:,compounds.index(a)];nr=treated_n[:,compounds.index(b)];ref=b
        critical=stats.t.ppf(.975,df+d0)
        protein_frames.append(pd.DataFrame({'compound_id':label,'feature_id':features,'reference':ref,
            'n_treated_obs':nt,'n_reference_obs':nr,'eligible_for_testing':valid[:,ci],
            'log2fc':B[:,ci],'se_mod':SE[:,ci],'ci95_low':B[:,ci]-critical*SE[:,ci],
            'ci95_high':B[:,ci]+critical*SE[:,ci],'t_mod':T[:,ci],'p_nominal':P[:,ci],
            'q_BH_contrast':bh(P[:,ci])}))
    protein=pd.concat(protein_frames[:len(compounds)],ignore_index=True)
    direct=pd.concat(protein_frames[len(compounds):],ignore_index=True)
    write_table(protein,ROOT/'reanalysis/A02/protein_effects.tsv.gz')
    write_protein_partitions(protein)
    write_table(direct,ROOT/'reanalysis/A02/direct_comparator_protein_effects.tsv.gz')
    gmt=load_gmt();correlations=residual_correlations(resid,features,gmt)
    primary=pathway_results(T[:,:len(compounds)],B[:,:len(compounds)],features,compounds,gmt,correlations,'primary_plate_blocked_eb',True)
    plate,units=plate_program_comparisons(A,samples,features,{c:np.isfinite(T[:,i]) for i,c in enumerate(compounds)},
        [(c,c,'DMSO') for c in compounds],gmt)
    primary=primary.merge(plate,on=['compound_id','program_id'],validate='one_to_one')
    write_table(units,ROOT/'reanalysis/A02/plate_program_observations.tsv')
    write_table(primary,ROOT/'reanalysis/A02/program_effects.tsv')
    pair_table=pathway_results(T[:,len(compounds):],B[:,len(compounds):],features,labels[len(compounds):],gmt,correlations,'direct_kinase_comparators',False)
    pair_plate,pair_units=plate_program_comparisons(A,samples,features,
        {labels[len(compounds)+j]:np.isfinite(T[:,len(compounds)+j]) for j in range(len(direct_pairs))},
        [(labels[len(compounds)+j],a,b) for j,(a,b) in enumerate(direct_pairs)],gmt)
    pair_table=pair_table.merge(pair_plate,on=['compound_id','program_id'],validate='one_to_one')
    write_table(pair_units,ROOT/'reanalysis/A02/direct_plate_program_observations.tsv')
    write_table(pair_table,ROOT/'reanalysis/A02/direct_comparator_program_effects.tsv')
    # Re-estimate the prior on the stricter feature universe; no new imputations.
    strict=mapping.loc[mapping.retained_50pct_filter,'retained_20pct_filter'].to_numpy(bool)
    post20,d020,s020=moderate_variances(V,df,strict)
    T20=np.where(valid,(B/np.sqrt(U*post20[:,None])),np.nan)[strict,:len(compounds)]
    corr20=residual_correlations(resid[:,strict],features[strict],gmt)
    sensitivity=pathway_results(T20,B[strict,:len(compounds)],features[strict],compounds,gmt,corr20,'max20pct_missing_min3',False)
    strict_plate,_=plate_program_comparisons(A[:,strict],samples,features[strict],
        {c:np.isfinite(T20[:,i]) for i,c in enumerate(compounds)},[(c,c,'DMSO') for c in compounds],gmt)
    sensitivity=sensitivity.merge(strict_plate,on=['compound_id','program_id'],validate='one_to_one')
    write_table(sensitivity,ROOT/'reanalysis/A04/sensitivity_effects.tsv')
    # Cache only intermediates. The tracked outputs and raw matrix suffice to rebuild.
    cache=ROOT/'reanalysis/_cache';cache.mkdir(exist_ok=True)
    np.savez_compressed(cache/'model.npz',effects=B,unscaled=U,variances=V,df_residual=df,residuals=resid,features=features.astype(str),contrast_names=np.array(labels),eligible=valid)
    diagnostics=[]
    for pid,(idx,corr) in correlations.items():
        rho,fraction=empirical_rho(corr,np.ones(len(idx),bool))
        diagnostics.append({'program_id':pid,'n_program_features':len(idx),'rho_pairwise_residual':rho,'finite_pair_fraction':fraction,
            'complete_residual_rows':int(np.isfinite(resid[:,idx]).all(axis=1).sum()),'min_pair_observations':20})
    write_table(pd.DataFrame(diagnostics),ROOT/'reanalysis/A02/correlation_diagnostics.tsv')
    config={'primary_model':'observed-data OLS with additive compound and plate effects; parametric empirical Bayes variance moderation',
        'primary_pathway_method':'one-sample t test of competitive mean log2 fold changes across independent matched plates',
        'secondary_pathway_method':'Python implementation of the non-rank CAMERA-PR formula on precomputed moderated t statistics',
        'inter_protein_correlation_fixed':PRIMARY_RHO,'empirical_correlation':'mean pairwise residual Pearson correlation; min 20 pair observations; estimates unavailable below 80% finite pairs',
        'n_primary_samples':len(samples),'n_features':len(features),'design_rank':int(np.linalg.matrix_rank(D)),
        'n_primary_compounds':len(compounds),'n_primary_program_tests':int(primary.program_family.eq('primary_aging').sum()),
        'd0_prior':d0,'s0_2_prior':s0,'strict20_d0_prior':d020,'strict20_s0_2_prior':s020,
        'normalization':'log2 sample-median normalization, relative composition','minimum_treated_observations':MIN_TREATED_OBS,'minimum_reference_observations':MIN_CONTROL_OBS,
        'gsea_permutations':2000,'gsea_seed':42,'gsea_inference':'NES and leading edge are descriptive only; no nominal or FDR significance claims from gene-set permutations',
        'composite_geroprotection_score':'removed; no universal beneficial program directions',
        'low_observation_protein_contrasts_excluded':int((~valid[:,:len(compounds)]).sum())}
    json_write(config,ROOT/'reanalysis/A02/analysis_config.json')
    import yaml
    (ROOT/'reanalysis/A02/analysis_config.yaml').write_text(yaml.safe_dump(config,sort_keys=False))
    json_write({'primary_contrasts':[c+'_vs_DMSO' for c in compounds],'direct_comparators':[{'treated':a,'reference':b} for a,b in direct_pairs]},ROOT/'reanalysis/A02/contrast_definitions.json')
    versions={'Python':platform.python_version()}
    for name in ['numpy','pandas','scipy','h5py','gseapy','seaborn','matplotlib','openpyxl','python-docx']:
        versions[name]=importlib.metadata.version(name)
    (ROOT/'reanalysis/A02/session_info.txt').write_text('\n'.join(f'{k}: {v}' for k,v in versions.items())+'\n')
    report=write_inference_report(primary,config)
    print(report,flush=True);print(f'A02 complete in {time.monotonic()-start:.1f} s',flush=True)


if __name__=='__main__': main()
