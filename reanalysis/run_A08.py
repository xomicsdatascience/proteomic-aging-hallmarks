"""Draw current figures from immutable A01-A07 outputs. Never edit inference tables."""
import hashlib
from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from PIL import Image, PngImagePlugin
from analysis_common import (ROOT, PROGRAM_NAMES, PRIMARY_IDS, HIGHLIGHTS, KINASE_PANEL, SENTINELS,
    read_table, write_table, normalized_data, short_drug)
from prioritization import ASSAY_DIRECTIONS

FIG=ROOT/'reanalysis/A08/figures';SOURCE=ROOT/'reanalysis/A08/source_tables'
FIG.mkdir(parents=True,exist_ok=True);SOURCE.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':13,'axes.labelsize':12,
                     'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':11.5,
                     'svg.hashsalt':'proteomic-aging-hallmarks','savefig.dpi':300})


def inference_hashes():
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for directory in (ROOT/'reanalysis').glob('A0[1-7]')
            for p in directory.rglob('*') if p.is_file()}


def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=300,bbox_inches='tight',metadata={'Software':'proteomic-aging-hallmarks'})
    fig.savefig(FIG/(name+'.svg'),bbox_inches='tight',metadata={'Date':None})
    if name in ['fig2_highlighted_heatmap','fig3_nes_distributions']:
        # Keep native 300-dpi dimensions; full-color vectors remain available.
        with Image.open(FIG/(name+'.png')) as source:
            indexed=source.convert('RGB').quantize(colors=256,method=Image.Quantize.MEDIANCUT,dither=Image.Dither.NONE)
        metadata=PngImagePlugin.PngInfo();metadata.add_text('Software','proteomic-aging-hallmarks')
        indexed.save(FIG/(name+'.png'),optimize=True,dpi=(300,300),pnginfo=metadata)
    plt.close(fig)


def labels(ids):
    return [PROGRAM_NAMES[k] for k in ids]


def atlas(ax,t,q,title,ticksize=10):
    annotations=q.map(lambda v:'**' if v<.01 else ('*' if v<.05 else ''))
    sns.heatmap(t,cmap='RdBu_r',center=0,vmin=-10,vmax=10,annot=annotations,fmt='',xticklabels=True,yticklabels=True,
                annot_kws={'fontsize':9.5},linewidths=.25,linecolor='white',ax=ax,
                cbar_kws={'label':'Plate-level t (colors clipped at +/-10)','extend':'both'})
    ax.set_xticklabels(labels(t.columns),rotation=50,ha='right',fontsize=13)
    ax.set_yticklabels([short_drug(d) for d in t.index],fontsize=ticksize,rotation=0)
    ax.set_xlabel('');ax.set_ylabel('');ax.set_title(title,pad=12)


def main():
    before=inference_hashes()
    p=read_table('reanalysis/A02/program_effects.tsv');primary=p[p.program_id.isin(PRIMARY_IDS)]
    pman=read_table('reanalysis/A03/program_manifest.tsv').set_index('program_id').loc[PRIMARY_IDS]
    sm=read_table('reanalysis/A01/sample_manifest.tsv')
    specific=read_table('reanalysis/A06/candidate_specificity.tsv')
    composition=read_table('reanalysis/A06/composition_sensitivity.tsv')
    priorities=read_table('reanalysis/A06/exploratory_prioritization.tsv').set_index('compound_id')
    tp=primary.pivot(index='compound_id',columns='program_id',values='plate_t')[PRIMARY_IDS]
    qp=primary.pivot(index='compound_id',columns='program_id',values='q_BH_global_plate')[PRIMARY_IDS]
    # Main figure 1: coverage, with no assumed favorable directions.
    fig,ax=plt.subplots(figsize=(8.8,6.3));y=np.arange(len(pman))
    measured=pman.measured_retained_7871;total=pman.total_genes
    ax.hlines(y,measured,total,color='#b8b8b8',lw=2)
    ax.scatter(total,y,s=60,color='#777777',label='Annotated genes')
    ax.scatter(measured,y,s=65,color='#0072B2',label='Measured proteins')
    for i,(m,n) in enumerate(zip(measured,total)):
        ax.text(n+6,i,f'{int(m)}/{int(n)} ({100*m/n:.0f}%)',va='center',fontsize=10.5)
    ax.set_yticks(y,labels(PRIMARY_IDS));ax.invert_yaxis();ax.set_xlim(0,390)
    ax.set_xlabel('Number of genes or protein features');ax.legend(loc='lower right',frameon=False)
    fig.tight_layout();save(fig,'fig1_program_coverage')
    write_table(pman.reset_index(),SOURCE/'fig1_program_coverage.tsv')
    # Main figure 2: the combined index is exploratory and has no significance test.
    highlight=[d for d in HIGHLIGHTS if d in tp.index]
    fig,axes=plt.subplots(1,2,figsize=(12.0,9.0),gridspec_kw={'width_ratios':[4.7,1.5]})
    atlas(axes[0],tp.loc[highlight],qp.loc[highlight],'A. Illustrative program responses',12)
    scores=priorities.loc[highlight]
    y=np.arange(len(scores))+.5
    axes[1].hlines(y,scores.combined_response_index,scores.combined_response_index_growth_proxy_residual,color='#b8b8b8',lw=1)
    axes[1].scatter(scores.combined_response_index,y,color='#222222',s=28,label='Raw index',zorder=3)
    axes[1].scatter(scores.combined_response_index_growth_proxy_residual,y,facecolors='none',edgecolors='#0072B2',marker='D',s=28,label='E2F/G2M residual',zorder=3)
    axes[1].set_ylim(len(scores),0);axes[1].set_yticks([]);axes[1].axvline(0,color='#777777',ls='--',lw=.8)
    axes[1].set_title('B. 14-program index',pad=12)
    axes[1].set_xlabel('Mean signed effect / SD\n(descriptive)',fontsize=11)
    axes[1].legend(loc='upper center',bbox_to_anchor=(.5,-.10),frameon=False,fontsize=9.5)
    fig.tight_layout();save(fig,'fig2_highlighted_heatmap')
    write_table(primary[primary.compound_id.isin(highlight)],SOURCE/'fig2_highlighted_programs.tsv')
    write_table(scores.reset_index(),SOURCE/'fig2_exploratory_indices.tsv')
    write_table(primary.sort_values(['compound_id','program_id']),SOURCE/'Supplementary_Table_S1_all_compounds.tsv')
    # Main figure 3: descriptive NES distributions, with two response examples identified.
    fig,ax=plt.subplots(figsize=(8.8,6.8));rng=np.random.default_rng(42)
    data=[primary.loc[primary.program_id.eq(k),'nes_gsea'].dropna().to_numpy() for k in PRIMARY_IDS]
    bp=ax.boxplot(data,positions=np.arange(len(PRIMARY_IDS)),vert=False,widths=.55,showfliers=False,patch_artist=True)
    for box in bp['boxes']:box.set(facecolor='#eeeeee',edgecolor='#777777')
    fig3_sources=[]
    for i,k in enumerate(PRIMARY_IDS):
        group=primary[primary.program_id.eq(k)].copy()
        group['primary_q_lt_005']=group.q_BH_global_plate.lt(.05)
        group['point_y']=i+rng.uniform(-.18,.18,len(group))
        group['point_color']=np.where(group.primary_q_lt_005,'#CC79A7','#999999')
        group['point_marker']=np.where(group.primary_q_lt_005,'D','o')
        for d,c in [('Binimetinib','#009E73'),('Baicalein','#D55E00')]:
            group.loc[group.compound_id.eq(d),'point_color']=c
        for color in ['#999999','#CC79A7','#009E73','#D55E00']:
            for marker in ['o','D']:
                shown=group[group.point_color.eq(color)&group.point_marker.eq(marker)]
                candidate=color in ['#009E73','#D55E00']
                ax.scatter(shown.nes_gsea,shown.point_y,s=46 if candidate else 13,c=color,marker=marker,
                           alpha=1 if candidate else .65,edgecolors='white' if candidate else 'none',
                           linewidths=.4,zorder=4 if candidate else 2)
        fig3_sources.append(group)
    ax.axvline(0,color='black',lw=.8,ls='--');ax.set_yticks(range(len(PRIMARY_IDS)),labels(PRIMARY_IDS));ax.invert_yaxis()
    ax.set_xlabel('Normalized enrichment score (descriptive)')
    ax.legend(handles=[Line2D([],[],marker='D',ls='',color='#CC79A7',label='Plate-level q < 0.05'),Line2D([],[],marker='o',ls='',color='#999999',label='Other / unavailable')]+[Line2D([],[],marker='o',ls='',color=c,label=d) for d,c in [('Binimetinib','#009E73'),('Baicalein','#D55E00')]],loc='upper center',bbox_to_anchor=(.55,1.14),ncol=2,frameon=False,fontsize=10.5)
    fig.tight_layout();save(fig,'fig3_nes_distributions')
    write_table(pd.concat(fig3_sources)[['compound_id','program_id','nes_gsea','q_BH_global_plate','n_genes_measured','plate_test_status','primary_q_lt_005','point_color','point_marker','point_y']],SOURCE/'fig3_nes_all_compounds.tsv')
    # Main figure 4: within-plate vehicle subtraction, then pooled unit-variance scaling.
    Y,obs,features=normalized_data();panel=['DMSO']+KINASE_PANEL
    m=sm.compound_id.isin(panel)&~sm.excluded&sm.plate.ne('RedMix');observations=[]
    for gene in sum(SENTINELS.values(),[]):
        j=list(features).index(gene)
        values=[]
        for idx in np.where(m)[0]:
            dmso=sm.compound_id.eq('DMSO')&sm.plate.eq(sm.iloc[idx].plate)
            delta=Y[idx,j]-np.nanmean(Y[dmso,j])
            values.append({'sample_id':sm.iloc[idx].sample_id,'compound_id':sm.iloc[idx].compound_id,'plate':sm.iloc[idx].plate,
                           'feature_id':gene,'normalized_log2_abundance':Y[idx,j],'within_plate_log2fc':delta})
        frame=pd.DataFrame(values);frame['scaled_abundance']=(frame.within_plate_log2fc-frame.within_plate_log2fc.mean())/frame.within_plate_log2fc.std(ddof=1)
        observations.append(frame)
    values=pd.concat(observations,ignore_index=True);write_table(values,SOURCE/'fig4_scaled_protein_observations.tsv')
    colors=['#777777','#009E73','#56B4E9','#CC79A7','#D55E00','#0072B2','#E69F00'];fig,axes=plt.subplots(2,1,figsize=(10,8.0),sharey=False)
    rng=np.random.default_rng(42)
    for ax,(title,genes),letter in zip(axes,SENTINELS.items(),['A','B']):
        for j,gene in enumerate(genes):
            for di,drug in enumerate(panel):
                v=values.loc[values.feature_id.eq(gene)&values.compound_id.eq(drug),'scaled_abundance'].dropna().to_numpy()
                x=j+(di-(len(panel)-1)/2)*.115
                if len(v):
                    ax.boxplot([v],positions=[x],widths=.09,showfliers=False,manage_ticks=False,patch_artist=True,
                        boxprops={'facecolor':colors[di],'alpha':.18,'edgecolor':colors[di]},medianprops={'color':colors[di]},
                        whiskerprops={'color':colors[di]},capprops={'color':colors[di]})
                    ax.scatter(x+rng.uniform(-.025,.025,len(v)),v,s=14,c=colors[di],alpha=.8,zorder=3)
        ax.axhline(0,color='black',ls='--',lw=.8);ax.set_xticks(range(len(genes)),genes)
        ax.set_ylabel('Scaled abundance (Z score)');ax.set_title(letter+'. '+title,loc='left')
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color=c,label=short_drug(d)) for d,c in zip(panel,colors)],loc='lower center',ncol=3,frameon=False,fontsize=11)
    fig.tight_layout(rect=(0,.10,1,1));save(fig,'fig4_kinase_inhibitor_panel')
    # Supplemental figure S1: truthful diagnostic names, data-derived axes and counts.
    qc=read_table('reanalysis/A04/qc_metrics.tsv');qc=qc[qc.plate.ne('RedMix')]
    rep=read_table('reanalysis/A04/replication_metrics.tsv');rep=rep[rep.compound_id.isin(['DFO','Binimetinib','Baicalein','Nordihydroguaiaretic acid','Lonafarnib','Erlotinib (Hydrochloride)','Mobocertinib (succinate)'])]
    null=read_table('reanalysis/A04/null_calibration.tsv')
    accounting=read_table('reanalysis/A01/sample_accounting.tsv');accounting=accounting[accounting.plate.ne('RedMix')]
    fig,axes=plt.subplots(2,2,figsize=(10.8,8.2))
    for role,color in [('library_compound','#0072B2'),('analytical_qc','#009E73'),('vehicle_control','#E69F00'),('positive_assay_control','#CC79A7')]:
        axes[0,0].plot(accounting.plate,accounting[role],'o',label=role.replace('_',' '),c=color,ms=4)
    axes[0,0].tick_params(axis='x',rotation=60,labelsize=9);axes[0,0].set_ylabel('Observations');axes[0,0].set_title('A. Primary plate composition');axes[0,0].legend(fontsize=8,frameon=False)
    axes[0,1].plot(qc.plate,qc.qc_median_protein_sd,'o',c='#0072B2');axes[0,1].tick_params(axis='x',rotation=60,labelsize=9)
    axes[0,1].set_ylabel('Median protein SD (log2)');axes[0,1].set_title('B. Pooled analytical QC precision')
    rep=rep.sort_values('mean_lopo_pearson_r');y=np.arange(len(rep));axes[1,0].hlines(y,rep.min_lopo_pearson_r,rep.max_lopo_pearson_r,color='#999999')
    axes[1,0].scatter(rep.mean_lopo_pearson_r,y,c='#009E73');axes[1,0].set_yticks(y,[short_drug(d) for d in rep.compound_id],fontsize=9);axes[1,0].set_xlim(0,1)
    axes[1,0].set_xlabel('Leave-one-plate-out Pearson r');axes[1,0].set_title('C. Cross-plate response similarity')
    axes[1,1].scatter(null.plate_t,np.arange(len(null)),c='#777777');axes[1,1].axvline(0,c='black',lw=.8)
    axes[1,1].set_yticks(range(len(null)),labels(null.program_id),fontsize=8.5)
    limit=max(1,float(null.plate_t.abs().max())*1.2);axes[1,1].set_xlim(-limit,limit)
    axes[1,1].set_xlabel('Vehicle split plate-level t');axes[1,1].set_title(f'D. Vehicle diagnostic ({int((null.q_BH_global_plate<.05).sum())} BH hits)')
    fig.tight_layout();save(fig,'figS1_qc_and_null_calibration')
    # Supplemental figure S2: full vector atlas and readable continuation pages.
    ordered=sorted(tp.index,key=str.casefold)
    fig,ax=plt.subplots(figsize=(12,32));atlas(ax,tp.loc[ordered],qp.loc[ordered],f'Complete {len(ordered)}-contrast response atlas',8)
    fig.tight_layout();fig.savefig(FIG/'figS2_full_160_compound_atlas.pdf',bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
    save(fig,'figS2_full_160_compound_atlas')
    with PdfPages(FIG/'figS2_atlas_paginated.pdf',metadata={'CreationDate':None,'ModDate':None}) as pdf:
        for page,start in enumerate(range(0,len(ordered),41),1):
            ds=ordered[start:start+41];fig,ax=plt.subplots(figsize=(10.8,10.0))
            atlas(ax,tp.loc[ds],qp.loc[ds],f'Supplementary Figure S2, part {page}',12.5)
            fig.tight_layout();pdf.savefig(fig,bbox_inches='tight');save(fig,f'figS2_atlas_part_{page}')
    # Supplemental S3: annotation context, growth sensitivity, missingness and composition.
    context=read_table('reanalysis/A05/drug_reference_comparisons.tsv');strict=read_table('reanalysis/A04/sensitivity_effects.tsv')
    compare=primary.merge(strict[['compound_id','program_id','camera_t']],on=['compound_id','program_id'],suffixes=('_primary','_strict'))
    fig,axes=plt.subplots(2,2,figsize=(10.8,8.2))
    axes[0,0].scatter(context.demagalhaes_up_camera_t,context.demagalhaes_dn_camera_t,c='#777777',s=20,alpha=.65)
    axes[0,0].set_xlabel('Aging-up RNA prior t');axes[0,0].set_ylabel('Aging-down RNA prior t');axes[0,0].set_title('A. Panel annotation context')
    context_colors=axes[0,1].scatter(specific.NUIS_E2F_PROLIF,specific.PROG_INFLAMMATORY,c=specific.PROG_OXPHOS,cmap='RdBu_r',s=22)
    fig.colorbar(context_colors,ax=axes[0,1],label='OXPHOS secondary t',pad=.02)
    axes[0,1].set_xlabel('E2F response t');axes[0,1].set_ylabel('Inflammatory response t');axes[0,1].set_title('B. Growth and response context')
    for d,offset in [('Binimetinib',(-45,18)),('Baicalein',(8,8)),('DFO',(5,5))]:
        s=specific[specific.compound_id.eq(d)].iloc[0];axes[0,1].annotate(d,(s.NUIS_E2F_PROLIF,s.PROG_INFLAMMATORY),fontsize=8,xytext=offset,textcoords='offset points',arrowprops={'arrowstyle':'-','lw':.5,'color':'#555555'})
    axes[1,0].scatter(compare.camera_t_primary,compare.camera_t_strict,s=5,alpha=.25,c='#0072B2');axes[1,0].set_xlabel('Primary <=50% missingness t');axes[1,0].set_ylabel('Stricter <=20% missingness t');axes[1,0].set_title('C. Missingness sensitivity')
    axes[1,1].scatter(composition.cytosolic_ribo_mean_log2fc,composition.oxphos_mean_log2fc,c='#009E73',s=18,alpha=.6)
    axes[1,1].set_xlabel('Ribosomal mean log2FC');axes[1,1].set_ylabel('OXPHOS mean log2FC');axes[1,1].set_title('D. Relative composition diagnostic')
    for ax in axes.flat:ax.axhline(0,c='#999999',lw=.6);ax.axvline(0,c='#999999',lw=.6)
    fig.tight_layout();save(fig,'figS3_benchmarks_and_sensitivity')
    layout=read_table('reanalysis/A06/layout_context.tsv')
    write_table(layout,SOURCE/'figS4_layout_context.tsv')
    fig,axes=plt.subplots(2,2,figsize=(10.8,7.8))
    for i,group in enumerate(['Blue','Red']):
        for j,pid in enumerate(['PROG_OXPHOS','PROG_INFLAMMATORY']):
            subset=layout[layout.assay_group.eq(group)&layout.program_id.eq(pid)]
            matrix=subset.pivot(index='well_row',columns='well_column',values='plate_competitive_log2fc').reindex(index=list('ABCDEFGH'),columns=range(1,13))
            sns.heatmap(matrix,ax=axes[i,j],cmap='RdBu_r',center=0,vmin=-.5,vmax=.5,
                        linewidths=.4,linecolor='white',cbar_kws={'label':'Competitive mean log2FC','extend':'both'})
            axes[i,j].set_title(group+' plates: '+PROGRAM_NAMES[pid])
            axes[i,j].set_xlabel('Well column');axes[i,j].set_ylabel('Well row')
            if group=='Red':
                for column,label in [(2,'Baic'),(3,'Bini')]:axes[i,j].text(column-.5,7.5,label,ha='center',va='center',fontsize=8,color='black',bbox={'facecolor':'white','alpha':.75,'pad':.5,'edgecolor':'none'})
    fig.tight_layout();save(fig,'figS4_layout_context')
    # S5 shows each axis, including uncertainty, for the sign-filtered shortlist.
    eligible=priorities[priorities.four_axis_eligible_rank.notna()].sort_values('four_axis_eligible_rank')
    top=eligible.head(5);top_ids=list(top.index)
    four=primary[primary.program_id.isin(ASSAY_DIRECTIONS)].merge(
        priorities.reset_index(),on='compound_id',validate='many_to_one')
    four['target_direction']=four.program_id.map(ASSAY_DIRECTIONS)
    four['signed_effect']=four.plate_competitive_log2fc*four.target_direction
    four['signed_ci95_low']=np.where(four.target_direction.eq(1),four.plate_ci95_low,-four.plate_ci95_high)
    four['signed_ci95_high']=np.where(four.target_direction.eq(1),four.plate_ci95_high,-four.plate_ci95_low)
    write_table(four.sort_values(['compound_id','program_id']),ROOT/'reanalysis/_cache/four_axis_all_compounds.tsv')
    selected_columns=['compound_id','program_id','target_direction','plate_competitive_log2fc',
                      'plate_ci95_low','plate_ci95_high','q_BH_global_plate','n_effective_plates',
                      'four_axis_index','four_axis_eligible_rank','four_axis_weakest_signed_component',
                      'well','assay_group','signed_effect','signed_ci95_low','signed_ci95_high']
    write_table(four[four.compound_id.isin(top_ids)][selected_columns].sort_values(['four_axis_eligible_rank','program_id']),SOURCE/'Supplementary_Table_S2_four_axis_candidates.tsv')
    write_table(priorities.reset_index(),SOURCE/'Supplementary_Table_S3_exploratory_indices.tsv')
    fig,axes=plt.subplots(2,2,figsize=(10.8,7.2),sharex=True)
    titles=['A. OXPHOS (+)','B. Macroautophagy machinery (+)',
            'C. Ubiquitination / proteasome (+)','D. Senescence-associated machinery (-)']
    for ax,(pid,direction),title in zip(axes.flat,ASSAY_DIRECTIONS.items(),titles):
        shown=four[four.program_id.eq(pid)].set_index('compound_id').loc[top_ids]
        y=np.arange(len(shown))
        ax.errorbar(shown.signed_effect,y,xerr=np.vstack([shown.signed_effect-shown.signed_ci95_low,shown.signed_ci95_high-shown.signed_effect]),fmt='none',ecolor='#777777',capsize=3,lw=1.2,zorder=1)
        for hit,marker in [(False,'o'),(True,'D')]:
            m=shown.q_BH_global_plate.lt(.05).eq(hit)
            ax.scatter(shown.loc[m,'signed_effect'],y[m],marker=marker,c='#CC79A7' if hit else '#777777',s=45,zorder=3)
        ax.axvline(0,c='black',ls='--',lw=.8);ax.set_yticks(y,[short_drug(d) for d in top_ids],fontsize=10.5);ax.invert_yaxis()
        ax.set_title(title,loc='left',fontsize=12)
        ax.set_xlim(-.065,.15)
    fig.supxlabel('Signed competitive mean log2FC (positive follows requested direction)',fontsize=12)
    fig.legend(handles=[Line2D([],[],marker='D',ls='',c='#CC79A7',label='Primary global q < 0.05'),Line2D([],[],marker='o',ls='',c='#777777',label='Other / unavailable')],loc='upper center',ncol=2,frameon=False,fontsize=11)
    fig.tight_layout(rect=(0,.035,1,.94));save(fig,'figS5_four_axis_candidates')
    after=inference_hashes()
    if before!=after: raise RuntimeError('Figure generation changed A01-A07 inference outputs')
    print('Figures regenerated; every A01-A07 file is unchanged.',flush=True)

if __name__=='__main__':main()
