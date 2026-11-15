"""Rebuild provenance and coverage using the committed, frozen gene sets."""
import hashlib
import numpy as np
import pandas as pd
from analysis_common import ROOT, PROGRAM_NAMES, canonical_compound, load_raw, load_gmt, sha256, write_table, json_write

X, obs, features = load_raw()
obs['compound_id'] = obs.drug
obs['plate'] = obs.batch
obs['well'] = obs.sample_id.str.split('_', n=1).str[1]
roles = {'QC_sample':'analytical_qc', 'DMSO':'vehicle_control', 'DFO':'positive_assay_control', 'Drug':'unlabeled_excluded'}
obs['condition_role'] = obs.drug.map(roles).fillna('library_compound')
obs['excluded'] = obs.drug.eq('Drug')
obs['exclusion_reason'] = np.where(obs.excluded, 'Unlabeled well', '')
obs['dose'] = np.where(obs.condition_role.eq('library_compound'), '10 uM', 'not available for control/QC')
obs['duration'] = np.where(obs.condition_role.ne('analytical_qc'), '24 h (parent study)', 'not applicable')
obs['culture_id'] = obs.sample_id
obs['culture_id_provenance'] = 'derived plate/well identifier; not an independently verified culture batch'
for key in ['culture_day','preparation_batch','instrument_batch','injection_order']:
    obs[key] = 'not available in deposited matrix'
obs['injection_id'] = obs.sample_id
cols = ['sample_id','original_compound_id','compound_id','condition_role','dose','duration','culture_id','culture_id_provenance','culture_day','plate','well','preparation_batch','instrument_batch','injection_order','injection_id','excluded','exclusion_reason']
write_table(obs[cols], ROOT/'reanalysis/A01/sample_manifest.tsv')
write_table(obs[['original_compound_id','compound_id']].drop_duplicates(), ROOT/'reanalysis/A01/compound_identifier_mapping.tsv')
metadata = pd.read_excel(ROOT/'data/drug_details (1).xlsx')
metadata['compound_id'] = metadata.ProductName.map(canonical_compound)
compounds = obs.groupby(['compound_id','condition_role'], sort=True).agg(n_samples=('sample_id','size'),n_plates=('plate','nunique'),plates=('plate',lambda s: ','.join(sorted(s.unique())))).reset_index()
compounds = compounds.merge(metadata[['compound_id','BaseName','CAS Number','Category','Target','PathWay','Research Area']], how='left', on='compound_id', validate='one_to_one').rename(columns={'CAS Number':'CAS_Number','Research Area':'Research_Area'})
write_table(compounds, ROOT/'reanalysis/A01/compound_manifest.tsv')
valid = ~obs.excluded.to_numpy()
zero_rate = (X[valid] == 0).mean(axis=0)
nonqc = valid & obs.condition_role.ne('analytical_qc').to_numpy()
mapping = pd.DataFrame({'feature_id':features,'gene_symbol':features,'zero_fraction_all':zero_rate,'zero_fraction_non_qc':(X[nonqc]==0).mean(axis=0),'retained_50pct_filter':zero_rate<=.5,'retained_20pct_filter':zero_rate<=.2})
write_table(mapping, ROOT/'reanalysis/A01/protein_mapping.tsv')
accounting = obs.groupby(['plate','condition_role']).size().unstack(fill_value=0).reset_index()
accounting['total'] = accounting.drop(columns='plate').sum(axis=1)
write_table(accounting, ROOT/'reanalysis/A01/sample_accounting.tsv')
dmso_plates = set(obs.loc[obs.compound_id.eq('DMSO'),'plate'])
eligible=[]
for row in compounds.itertuples():
    if row.condition_role not in ['library_compound','positive_assay_control']: continue
    plates=row.plates.split(','); status='directly_estimable' if len(plates)>1 and all(p in dmso_plates for p in plates) else 'excluded_no_within_plate_vehicle'
    eligible.append({'compound_id':row.compound_id,'contrast_id':row.compound_id+'_vs_DMSO','condition_role':row.condition_role,'n_treated':row.n_samples,'n_plates':row.n_plates,'plates':row.plates,'eligibility_status':status,'eligibility_reason':'Repeated treated wells on plates containing vehicle controls' if status=='directly_estimable' else 'No matched within-plate vehicle in deposited data'})
write_table(pd.DataFrame(eligible), ROOT/'reanalysis/A01/contrast_eligibility.tsv')
# Freeze membership from version-controlled GMT, without a network lookup.
gmt = load_gmt(); retained=set(features[zero_rate<=.5]); all_features=set(features)
pman=pd.read_csv(ROOT/'reanalysis/A03/program_manifest.tsv',sep='\t')
mem=[]; overlaps=[]
for idx,row in pman.iterrows():
    pid=row.program_id; genes=gmt[pid]
    pman.loc[idx,'program_name']=PROGRAM_NAMES[pid]
    pman.loc[idx,'total_genes']=len(genes)
    pman.loc[idx,'measured_all_8630']=len(genes & all_features)
    pman.loc[idx,'measured_retained_7871']=len(genes & retained)
    pman.loc[idx,'coverage_pct_retained']=100*len(genes & retained)/len(genes)
    pman.loc[idx,'membership_checksum']=hashlib.sha256(','.join(sorted(genes)).encode()).hexdigest()[:16]
    for g in sorted(genes): mem.append({'program_id':pid,'gene_symbol':g,'in_raw_8630':g in all_features,'in_retained_7871':g in retained})
for i,a in enumerate(gmt):
    for b in list(gmt)[i+1:]:
        ga,gb=gmt[a]&retained,gmt[b]&retained
        overlaps.append({'program_1':a,'program_2':b,'overlap_count':len(ga&gb),'jaccard':len(ga&gb)/len(ga|gb)})
write_table(pman, ROOT/'reanalysis/A03/program_manifest.tsv')
write_table(pd.DataFrame(mem), ROOT/'reanalysis/A03/program_membership.tsv')
write_table(pd.DataFrame(overlaps), ROOT/'reanalysis/A03/program_overlap.tsv')
summary={'raw_samples':len(obs),'unlabeled_excluded':int(obs.excluded.sum()),'normalization_samples':int(valid.sum()),'raw_features':len(features),'retained_features_50pct':int((zero_rate<=.5).sum()),'retained_features_20pct':int((zero_rate<=.2).sum()),'filter_order':'exclude unlabeled well, then compute feature missingness across the remaining observations','raw_h5ad_sha256':sha256(ROOT/'data/raw (1).h5ad'),'drug_metadata_sha256':sha256(ROOT/'data/drug_details (1).xlsx'),'gmt_sha256':sha256(ROOT/'reanalysis/A03/programs.gmt')}
json_write(summary,ROOT/'reanalysis/A01/input_summary.json')
(ROOT/'reanalysis/A01/design_audit.md').write_text(f"""# Design and provenance\n\nThe deposited matrix contains {len(obs):,} observations and {len(features):,} features. One unlabeled well is excluded before filtering. The <=50% zero-value filter retains {summary['retained_features_50pct']:,} features; the <=20% filter retains {summary['retained_features_20pct']:,}. This order excludes the boundary feature ZNF552 that was retained when filtering before the unlabeled-well exclusion.\n\nThe primary analysis compares 160 library entries and DFO with within-plate DMSO controls on 12 plates. RedMix is excluded from primary inference. Culture and acquisition metadata not in the deposited file remain marked unavailable. Plate/well identifiers do not verify independent culture days. The parent study calls these biological replicates.\n\nInput hashes and derived counts are in input_summary.json. Compound aliases are resolved once in compound_identifier_mapping.tsv.\n""")
(ROOT/'reanalysis/A03/annotation_decisions.md').write_text("""# Annotation scope\n\nThe committed GMT defines 14 selected aging-associated programs and three growth/stress response modules. Membership is fixed for reproducibility; it was not prospectively registered. Whole-cell protein abundance does not directly measure DNA damage, telomere length, chromatin state, autophagic flux, irreversible senescence, secreted SASP, stem-cell exhaustion, or microbiome dysbiosis. TLR proteins describe host innate sensing. Context-dependent directions remain context dependent; no universal anti-aging signs or composite geroprotection score are assigned. The RNA-derived aging priors are part of this panel and are not independent validation.\n""")
(ROOT/'reanalysis/A01/parent_study_overlap.tsv').write_text('item\tstatus\tnotes\nHepG2 screen and acquisition\treused\tJiang et al. 2026, DOI 10.64898/2026.05.04.722718\nSelected aging-associated program comparisons\tnew reanalysis\tDescriptive protein responses, observed-data plate model, correlation and missingness sensitivity\nClinical or organismal validation\tnot performed\tNo human organ-age or longevity claim\n')
print(summary,flush=True)
