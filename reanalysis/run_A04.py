"""Vehicle diagnostics, replication, missingness and abundance stability."""
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from analysis_common import (ROOT, PRIMARY_IDS, HIGHLIGHTS, PRIMARY_RHO, normalized_data, read_table,
    write_table, json_write, load_gmt, primary_mask, competitive_test, bh, plate_program_comparisons)


def main():
    sm=read_table('reanalysis/A01/sample_manifest.tsv');Y,obs,features=normalized_data();gmt=load_gmt()
    prog=read_table('reanalysis/A02/program_effects.tsv')
    protein=read_table('reanalysis/A02/protein_effects.tsv.gz')
    qc=[]
    for plate in sorted(sm.plate.unique()):
        ix=sm.plate.eq(plate).to_numpy();q=ix&sm.condition_role.eq('analytical_qc').to_numpy()
        qc.append({'plate':plate,'n_qc':int(q.sum()),'n_dmso':int((ix&sm.compound_id.eq('DMSO').to_numpy()).sum()),
                   'qc_missing_rate':float(np.isnan(Y[q]).mean()),'qc_median_protein_sd':float(np.nanmedian(np.nanstd(Y[q],axis=0,ddof=1)))})
    write_table(pd.DataFrame(qc),ROOT/'reanalysis/A04/qc_metrics.tsv')
    reps=[];stability=[];rng=np.random.default_rng(42)
    for drug in [d for d in HIGHLIGHTS if d in set(sm.compound_id)]:
        fcs=[]
        for plate in sorted(sm.loc[sm.compound_id.eq(drug)&sm.plate.ne('RedMix'),'plate'].unique()):
            treated=sm.compound_id.eq(drug)&sm.plate.eq(plate);control=sm.compound_id.eq('DMSO')&sm.plate.eq(plate)
            if not control.any():continue
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                fcs.append(np.nanmean(Y[treated],axis=0)-np.nanmean(Y[control],axis=0))
        fcs=np.array(fcs)
        if len(fcs)<2:continue
        cor=[]
        for j in range(len(fcs)):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                rest=np.nanmean(np.delete(fcs,j,axis=0),axis=0)
            ok=np.isfinite(fcs[j])&np.isfinite(rest)
            cor.append(float(stats.pearsonr(fcs[j,ok],rest[ok]).statistic))
        reps.append({'compound_id':drug,'n_plates':len(fcs),'mean_lopo_pearson_r':np.mean(cor),'min_lopo_pearson_r':min(cor),'max_lopo_pearson_r':max(cor)})
        eligible=set(protein.loc[protein.compound_id.eq(drug)&protein.eligible_for_testing,'feature_id'])
        selections=rng.integers(0,len(fcs),size=(1000,len(fcs)))
        for pid in PRIMARY_IDS:
            member=np.isin(features,list(gmt[pid]&eligible))
            if member.sum()<5:continue
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                unit=np.nanmean(fcs[:,member],axis=1)
            if not np.isfinite(unit).all():continue
            boot=unit[selections].mean(axis=1);estimate=float(unit.mean())
            stability.append({'compound_id':drug,'program_id':pid,'statistic':'mean within-plate program log2 fold change','n_plates':len(fcs),'n_features':int(member.sum()),
                'mean_plate_log2fc':estimate,'boot_ci95_low':np.quantile(boot,.025),'boot_ci95_high':np.quantile(boot,.975),
                'sign_consistency_pct':100*np.mean(np.sign(boot)==np.sign(estimate))})
    write_table(pd.DataFrame(reps),ROOT/'reanalysis/A04/replication_metrics.tsv')
    write_table(pd.DataFrame(stability),ROOT/'reanalysis/A04/candidate_stability.tsv')
    missing=protein.groupby('compound_id').agg(n_features=('feature_id','size'),n_eligible_features=('eligible_for_testing','sum'),minimum_observed_treated=('n_treated_obs','min'),median_observed_treated=('n_treated_obs','median')).reset_index()
    missing['n_below_three_treated']=protein.assign(low=protein.n_treated_obs<3).groupby('compound_id').low.sum().reindex(missing.compound_id).to_numpy()
    write_table(missing,ROOT/'reanalysis/A04/contrast_missingness.tsv')
    # Compare the two vehicle positions on the same set of paired-control plates.
    compounds=sorted(protein.compound_id.unique())
    eligible_features={c:np.isin(features,protein.loc[protein.compound_id.eq(c)&protein.eligible_for_testing,'feature_id']) for c in compounds}
    paired_plates=[];first=[];second=[]
    for plate in sorted(sm.plate.unique()):
        ix=np.where(sm.plate.eq(plate)&sm.compound_id.eq('DMSO'))[0]
        if len(ix)==2:paired_plates.append(plate);first.append(ix[0]);second.append(ix[1])
    reference_sensitivity=[]
    for reference,indices in [('first_vehicle',first),('second_vehicle',second),('mean_paired_vehicles',first+second)]:
        keep=primary_mask(sm).to_numpy()&sm.plate.isin(paired_plates).to_numpy()
        keep&=~sm.compound_id.eq('DMSO').to_numpy()|np.isin(np.arange(len(sm)),indices)
        frame,_=plate_program_comparisons(Y[keep],sm.loc[keep].reset_index(drop=True),features,eligible_features,
            [(c,c,'DMSO') for c in compounds],{k:gmt[k] for k in PRIMARY_IDS})
        frame['vehicle_reference']=reference;reference_sensitivity.append(frame)
    write_table(pd.concat(reference_sensitivity,ignore_index=True),ROOT/'reanalysis/A04/vehicle_reference_sensitivity.tsv')
    # The paired vehicle null is always computed afresh, never rescaled from a saved statistic.
    a=[];b=[]
    for plate in sorted(sm.plate.unique()):
        ix=np.where(sm.plate.eq(plate)&sm.compound_id.eq('DMSO'))[0]
        if len(ix)==2:a.append(ix[0]);b.append(ix[1])
    diff=Y[a]-Y[b];count=np.isfinite(diff).sum(axis=0);good=count>=3
    diff0=np.nan_to_num(diff,nan=0)[:,good];n=count[good];ss=np.nansum(diff[:,good]**2,axis=0)
    signs=2*((np.arange(2**len(a))[:,None]>>np.arange(len(a)))&1)-1
    sums=signs@diff0
    with np.errstate(invalid='ignore',divide='ignore'):
        variances=(ss-sums**2/n)/(n-1)
        null_t=(sums/n)/np.sqrt(variances/n)
    null_t[~np.isfinite(null_t)]=np.nan
    diagnostic=[];scores=[]
    for pid in PRIMARY_IDS:
        member=np.isin(features[good],list(gmt[pid]));one=[]
        for ts in null_t:
            val,p,m,vif=competitive_test(ts,member,PRIMARY_RHO);one.append((val,p))
        one=np.array(one);scores.append(one)
        diagnostic.append({'program_id':pid,'null_sd':float(np.nanstd(one[:,0],ddof=1)),'nominal_p05_fraction':float(np.mean(one[:,1]<.05))})
    scores=np.stack(scores,axis=1);q=np.array([bh(v[:,1]) for v in scores]);actual=scores[-1]
    null=pd.DataFrame({'program_id':PRIMARY_IDS,'null_camera_t':actual[:,0],'null_p_camera':actual[:,1],'null_q_BH':q[-1],
                       'n_paired_plates':len(a),'rho_fixed':PRIMARY_RHO})
    pseudo=pd.DataFrame({'compound_id':['Vehicle_A']*len(a)+['Vehicle_B']*len(b),
                         'plate':sm.iloc[a].plate.tolist()+sm.iloc[b].plate.tolist()})
    paired,paired_units=plate_program_comparisons(np.vstack([Y[a],Y[b]]),pseudo,features,
        {'DMSO_split':good},[('DMSO_split','Vehicle_A','Vehicle_B')],{k:gmt[k] for k in PRIMARY_IDS})
    null=null.merge(paired[['program_id','plate_t','p_plate','q_BH_global_plate']],on='program_id',validate='one_to_one')
    write_table(paired_units,ROOT/'reanalysis/A04/null_plate_observations.tsv')
    positions=pd.concat([sm.iloc[a].assign(vehicle_position='first listed vehicle well'),sm.iloc[b].assign(vehicle_position='second listed vehicle well')])
    write_table(positions[['sample_id','plate','well','vehicle_position']],ROOT/'reanalysis/A04/vehicle_positions.tsv')
    write_table(null,ROOT/'reanalysis/A04/null_calibration.tsv')
    write_table(pd.DataFrame(diagnostic),ROOT/'reanalysis/A04/null_signflip_diagnostics.tsv')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        _,protein_null_p=stats.ttest_1samp(diff[:,good],0,axis=0,nan_policy='omit')
    prim=prog[prog.program_id.isin(PRIMARY_IDS)]
    strict=read_table('reanalysis/A04/sensitivity_effects.tsv')
    joined=prim.merge(strict[['compound_id','program_id','camera_t','q_BH_global_camera']],on=['compound_id','program_id'],suffixes=('_primary','_strict'))
    finite=joined.camera_t_primary.notna()&joined.camera_t_strict.notna()
    plate_units=paired_units.pivot(index='plate',columns='program_id',values='competitive_mean_log2fc').reindex(columns=PRIMARY_IDS).to_numpy()
    plate_flip_p=[]
    for sign in signs:
        plate_flip_p.append(stats.ttest_1samp(plate_units*sign[:,None],0,axis=0,nan_policy='omit').pvalue)
    plate_flip_q=np.array([bh(pv) for pv in plate_flip_p])
    summary={'vehicle_split_significant_programs':int((null.q_BH_global_plate<.05).sum()),
             'vehicle_split_fixed_rho_significant_programs':int((null.null_q_BH<.05).sum()),
             'vehicle_split_max_abs_t':float(null.plate_t.abs().max()),
             'vehicle_split_significant_proteins':int((bh(protein_null_p)<.05).sum()),'vehicle_pair_count':len(a),
             'signflip_assignments':len(signs),'signflip_fraction_any_program_q05':float(np.any(q<.05,axis=1).mean()),
             'plate_signflip_fraction_any_program_q05':float(np.any(plate_flip_q<.05,axis=1).mean()),
             'strict_missingness_camera_pearson_r':float(stats.pearsonr(joined.loc[finite,'camera_t_primary'],joined.loc[finite,'camera_t_strict']).statistic),
             'strict_missingness_sign_agreement':float((np.sign(joined.loc[finite,'camera_t_primary'])==np.sign(joined.loc[finite,'camera_t_strict'])).mean()),
             'bootstrap_repetitions':1000,'bootstrap_unit':'plate; descriptive program mean log2 fold changes, not a refit of the moderated model',
             'calibration_scope':'Diagnostic paired vehicle experiment; not a full-screen type I error or FDR estimate'}
    json_write(summary,ROOT/'reanalysis/A04/diagnostic_summary.json')
    (ROOT/'reanalysis/A04/robustness_report.md').write_text(f"""# Robustness and diagnostic scope\n\nThe freshly reconstructed vehicle split has {summary['vehicle_split_significant_programs']} significant plate-level pathways, {summary['vehicle_split_fixed_rho_significant_programs']} fixed-rho competitive pathways, and {summary['vehicle_split_significant_proteins']} significant proteins at BH q < 0.05. The controls therefore show position-associated differences and do not establish strict calibration. Across all {len(signs):,} paired sign assignments, {100*summary['plate_signflip_fraction_any_program_q05']:.2f}% have at least one significant plate-level program. This small paired-control design is diagnostic, not full-screen FDR validation.\n\nThe stricter <=20% missingness universe gives secondary CAMERA-PR-style pathway-statistic r={summary['strict_missingness_camera_pearson_r']:.3f} and {100*summary['strict_missingness_sign_agreement']:.1f}% sign agreement with the primary universe. Per-compound observed counts are reported separately.\n\nBootstrap intervals describe mean within-plate protein-program log2 fold changes after resampling independent plates (1,000 draws). They do not establish function, causal specificity, a geroprotection rank, or a second independent validation dataset.\n""")
    print(summary,flush=True)

if __name__=='__main__':main()
