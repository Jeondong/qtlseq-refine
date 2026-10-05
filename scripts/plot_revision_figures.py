"""Render revised simulation panels from saved results and fixed first examples.

Publication exports are retained in figures/. These regenerated plots use the
same reported quantities; typography/layout can differ from author-edited files.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from descriptive_core import make_chrom,windows,WINDOWS

ROOT=Path(__file__).resolve().parents[1]
COLORS=['#28628a','#c18725','#b64343']

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/figures'));args=ap.parse_args()
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):
        for ext in ['png','svg']:fig.savefig(out/f'{name}.{ext}',dpi=200,bbox_inches='tight',facecolor='white')
        plt.close(fig)
    t=pd.read_csv(args.root/'data/supplementary/Supplementary_Data_S1_threshold_sim.tsv',sep='\t')
    p=pd.read_csv(args.root/'data/supplementary/Supplementary_Data_S2_peakshift_sim.tsv',sep='\t')
    fig,axes=plt.subplots(1,3,figsize=(10,3.4),sharey=True,layout='constrained');rng=np.random.RandomState(42)
    for _ in range(5):
        pos,val=make_chrom(rng)
        for j,pct in enumerate([0,95,99]):
            keep=np.ones(len(val),bool) if pct==0 else abs(val)>=np.percentile(abs(val),pct)
            _,mid,means,_=windows(pos[keep],val[keep],100000,10000,include_empty=True)
            axes[j].plot(mid/1e6,means,color=COLORS[j],alpha=.45,lw=.7)
    for j,ax in enumerate(axes):
        ax.axvline(15,c='black',ls='--',lw=.8);ax.axvspan(13.5,16.5,color='gray',alpha=.12)
        ax.set(xlim=(0,30),ylim=(-.2,1.05),xlabel='Position (Mb)',title=f'({chr(65+j)}) '+['All SNPs','Observed p95','Observed p99'][j])
    axes[0].set_ylabel('Mean ΔSNP-index');save(fig,'Figure_3_revised')
    fig,axes=plt.subplots(3,3,figsize=(9,8),layout='constrained')
    for r,(metric,label) in enumerate([('qtl_signal','Signal-region mean'),('gw_noise','SD of retained window means'),('discriminability','Mean / SD ratio')]):
        for c,w in enumerate(WINDOWS):
            ax=axes[r,c];data=[t[(t.window==w)&(t.threshold==th)][metric] for th in ['All SNPs','p95','p99']]
            bp=ax.boxplot(data,patch_artist=True,showfliers=False)
            for box,color in zip(bp['boxes'],COLORS):box.set_facecolor(color);box.set_alpha(.65)
            ax.set_xticks([1,2,3],['All','p95','p99']);ax.set_title(f'({chr(65+r*3+c)}) '+w.split('\n')[0]);ax.set_ylabel(label if c==0 else '')
    save(fig,'Figure_4_revised')
    fig,axes=plt.subplots(3,2,figsize=(10,9),layout='constrained');scales=[2,1,.5,.1,.01]
    pos,val=make_chrom(np.random.RandomState(2025));keep=abs(val)>=np.percentile(abs(val),99)
    for ax,mb,letter in zip(axes[0],[2,.1],['A','B']):
        width=int(mb*1e6);st,mid,means,_=windows(pos[keep],val[keep],width,max(int(width*.05),1000),include_empty=True)
        ax.plot(st/1e6,means,c=COLORS[2],label='Start');ax.plot(mid/1e6,means,c=COLORS[0],ls='--',label='Midpoint')
        peak=np.nanargmax(means);ax.axvline(st[peak]/1e6,c=COLORS[2],ls=':',lw=.8);ax.axvline(mid[peak]/1e6,c=COLORS[0],ls=':',lw=.8)
        ax.axvline(15,c='green',lw=.8);ax.set(xlim=(9,21),xlabel='Reported position (Mb)',ylabel='Mean ΔSNP-index',title=f'({letter}) First replicate, {mb:g} Mb');ax.legend(fontsize=8)
    for method,color,offset in [('start',COLORS[2],-.16),('center',COLORS[0],.16)]:
        data=[p[p.window_mb==w]['error_'+method] for w in scales]
        bp=axes[1,0].boxplot(data,positions=np.arange(5)+offset,widths=.27,patch_artist=True,showfliers=True)
        for box in bp['boxes']:box.set_facecolor(color);box.set_alpha(.6)
        med=[p[p.window_mb==w]['abs_err_'+method].median() for w in scales]
        axes[1,1].plot(np.arange(5),med,'o-',c=color,label=method.title())
    axes[1,0].axhline(0,c='gray',lw=.7);axes[1,0].set(title='(C) Signed peak error',ylabel='Error (Mb)')
    axes[1,1].set(title='(D) Median absolute peak error',ylabel='Absolute error (Mb)');axes[1,1].legend(fontsize=8)
    signed=[p[p.window_mb==w].error_start.median() for w in scales];axes[2,0].bar(np.arange(5),signed,color=COLORS[2])
    axes[2,0].set(title='(E) Median signed error, start assignment',ylabel='Signed error (Mb)')
    reduction=[abs(p[p.window_mb==w].error_start.median())-abs(p[p.window_mb==w].error_center.median()) for w in scales]
    axes[2,1].plot(np.arange(5),reduction,'o-',c='#695394');axes[2,1].set(title='(F) Reduction in absolute median signed error',ylabel='Directional-bias reduction (Mb)')
    for ax in axes[1:].flat:ax.set_xticks(np.arange(5),[str(w) for w in scales]);ax.set_xlabel('Window size (Mb)')
    save(fig,'Figure_8_white')
    summary=pd.read_csv(args.root/'data/supplementary/Supplementary_Data_S5/benchmark_verified_summary.tsv',sep='\t')
    protocol=json.loads((args.root/'data/supplementary/Supplementary_Data_S5/benchmark_protocol.json').read_text());conditions=list(protocol['configs'])
    fig,axes=plt.subplots(1,4,figsize=(12,6.5),sharey=True,layout='constrained')
    metrics=[('coverage','(A) True-position inclusion'),('median_peak_error_mb','(B) Peak error (Mb)'),('median_called_width_mb','(C) Called width (Mb)'),('false_positive_rate','(D) Null detection')]
    for j,(metric,title) in enumerate(metrics):
        ax=axes[j]
        for k,(meth,label) in enumerate([('ALL_100kb','All'),('p95_100kb','p95 null'),('p99_100kb','p99 null')]):
            ss=summary[summary.method==meth].set_index('scenario').loc[conditions];vals=ss[metric].to_numpy();y=np.arange(13)+(k-1)*.18
            ax.plot(vals,y,'o',ms=3.7,c=COLORS[k],label=label)
            if metric in ['coverage','false_positive_rate']:
                pre='coverage' if metric=='coverage' else 'null';ax.errorbar(vals,y,xerr=[vals-ss[pre+'_lower'],ss[pre+'_upper']-vals],fmt='none',ecolor=COLORS[k],lw=.6)
        ax.set_title(title,fontsize=9);ax.grid(axis='x',alpha=.15)
        if metric=='coverage':ax.set_xlim(-.035,1.035)
        if metric=='false_positive_rate':ax.axvline(.05,c='gray',ls='--',lw=.7)
    axes[0].set_yticks(np.arange(13),[c.replace('_',' ') for c in conditions]);axes[0].invert_yaxis();axes[-1].legend(fontsize=7,loc='lower right')
    save(fig,'Figure_S2_population_sensitivity')
    fig,axes=plt.subplots(1,3,figsize=(9,3.4),layout='constrained')
    for ax,sc,title in zip(axes,['baseline','weak_QTL','low_depth'],['Baseline','Weak effect','Low depth']):
        ss=summary[summary.scenario==sc].set_index('method').loc[['ALL_100kb','p99_100kb','p99_n1_100kb']]
        for offset,(metric,label,color) in zip([-.25,0,.25],[('detection_rate','QTL detected',COLORS[0]),('coverage','True position included',COLORS[2]),('no_finite_peak_rate','No finite peak','#888888')]):
            ax.bar(np.arange(3)+offset,ss[metric],width=.24,label=label,color=color)
        ax.set_xticks(np.arange(3),['All\nmin 10','p99\nmin 10','p99\nmin 1']);ax.set(ylim=(0,1.05),title=title,ylabel='Fraction of all QTL tests')
    axes[0].legend(fontsize=6.5);save(fig,'Figure_S1_marker_support')
    print(f'Rendered five simulation figures to {out}')

if __name__=='__main__':main()
