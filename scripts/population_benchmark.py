"""Seeded F2 benchmark with independent null calibration and held-out populations.
All comparator scores are implemented here; this is NOT execution of QTLseqr or PyBSASeq.
Intervals are chromosome-wise null-threshold exceedance sets, NOT location confidence intervals.
"""
from pathlib import Path
import sys, json, argparse, time, math
P=Path('outputs/f2_new')
import numpy as np
import pandas as pd
from scipy.stats import binom
from scipy.special import xlogy
L=30_000_000
BASE=dict(N=400,bulk=20,depth=40,density=100,rate=4.,h2=.3,map='uniform',error=.005,bias=0.)
SCENARIOS={
 'baseline':{},'weak_QTL':{'h2':.1},'strong_QTL':{'h2':.6},
 'small_bulk':{'bulk':10},'large_bulk':{'bulk':50},
 'low_depth':{'depth':15},'high_depth':{'depth':80},
 'low_recombination':{'rate':.5},'high_recombination':{'rate':8.},
 'heterogeneous_map':{'map':'heterogeneous'},'sparse_markers':{'density':10},
 'dense_markers':{'density':400},'differential_error':{'bias':.08}}
SCALE=[(2000000,100000),(500000,50000),(100000,10000)]

def exact_cutoffs(bulk,depth,error):
 """Exact single-marker absolute-difference null distribution, independent random F2 bulks."""
 k=np.arange(2*bulk+1);p=error+(1-2*error)*k/(2*bulk)
 pmf=(binom.pmf(np.arange(depth+1)[:,None],depth,p[None,:])*binom.pmf(k,2*bulk,.5)).sum(axis=1)
 diff=np.abs(np.arange(depth+1)[:,None]-np.arange(depth+1)[None,:])
 prob=np.bincount(diff.ravel(),weights=(pmf[:,None]*pmf[None,:]).ravel(),minlength=depth+1)
 return np.searchsorted(np.cumsum(prob),[.95,.99])/depth

def setup(c):
 m=int(30*c['density']);pos=(np.arange(m)+.5)*L/m
 rate=np.full(m,c['rate'])
 if c['map']=='heterogeneous':
  rate[(pos>=8e6)&(pos<18e6)]*=.05
  rate[(pos>=23e6)&(pos<25e6)]*=5
 # Genetic coordinates in Morgans. Piecewise constant physical map.
 g=np.cumsum(rate*(L/m)/1e8);g-=rate*(L/m)/2e8
 total=float(np.sum(rate*(L/m)/1e8))
 cuts=exact_cutoffs(c['bulk'],c['depth'],c['error'])
 grids={}
 for w,step in SCALE:
  s=np.arange(0,L-w+1,step)
  grids[w]=(s,s+w/2,np.searchsorted(pos,s),np.searchsorted(pos,s+w))
 # Explicit tricube G smoothing with half-width 1 Mb, not locfit/package p-values.
 radius=1e6;off=np.arange(-math.ceil(radius/(L/m)),math.ceil(radius/(L/m))+1)*L/m
 kernel=np.maximum(1-(abs(off)/radius)**3,0)**3
 return dict(pos=pos,g=g,total=total,cuts=cuts,grids=grids,kernel=kernel)

def simulate(c,s,rng,alt,diagnostics=False):
 """Poisson F1 crossovers -> F2 genotypes -> phenotype extremes -> pooled reads."""
 n=c['N'];m=len(s['pos']);nh=2*n
 ne=rng.poisson(s['total'],nh)
 h=np.zeros((nh,m),np.uint8)
 rows=np.repeat(np.arange(nh),ne);events=rng.uniform(0,s['total'],len(rows))
 ix=np.searchsorted(s['g'],events);keep=ix<m
 np.add.at(h,(rows[keep],ix[keep]),1)
 h[:,0]+=rng.integers(0,2,nh,dtype=np.uint8)
 np.cumsum(h,axis=1,dtype=np.uint8,out=h);h&=1
 geno=h.reshape(n,2,m).sum(axis=1,dtype=np.uint8)
 qi=int(rng.integers(np.searchsorted(s['pos'],3e6),np.searchsorted(s['pos'],27e6)))
 truth=float(s['pos'][qi])
 # Additive genetic variance = .5 for unselected F2; h2 means expected single-QTL PVE.
 if alt: phenotype=geno[:,qi]+rng.normal(0,np.sqrt(.5*(1-c['h2'])/c['h2']),n)
 else: phenotype=rng.normal(size=n)
 order=np.argsort(phenotype);b=c['bulk'];lo=order[:b];hi=order[-b:]
 p1=geno[lo].mean(axis=0)/2;p2=geno[hi].mean(axis=0)/2
 p1=c['error']+(1-2*c['error'])*p1;p2=c['error']+(1-2*c['error'])*p2
 if c['bias']:
  # Bulk-specific alt-read contamination in a fixed non-target-defined region.
  mask=(s['pos']>=3e6)&(s['pos']<6e6);p2[mask]=(1-c['bias'])*p2[mask]+c['bias']
 a1=rng.binomial(c['depth'],p1);a2=rng.binomial(c['depth'],p2)
 if diagnostics:
  return a1,a2,truth,geno
 return a1,a2,truth

def window_mean(v,keep,l,r,minn):
 count=np.r_[0,np.cumsum(keep.astype(int))];cnt=count[r]-count[l]
 cs=np.r_[0,np.cumsum(np.where(keep,v,0.))]
 val=np.divide(cs[r]-cs[l],cnt,out=np.full(len(l),np.nan),where=cnt>=minn)
 return val,cnt

def score(a1,a2,c,s):
 depth=c['depth'];delta=(a2-a1)/depth
 result={}
 for w,step in SCALE:
  starts,mid,l,r=s['grids'][w]
  for filt,keep,minn in [('ALL',np.ones(len(delta),bool),10),('p95',abs(delta)>s['cuts'][0],10),('p99',abs(delta)>s['cuts'][1],10),('p99_n1',abs(delta)>s['cuts'][1],1)]:
   val,cnt=window_mean(delta,keep,l,r,minn)
   result[f'{filt}_{w//1000}kb']=(np.abs(val),starts,starts+w,cnt)
 # G statistic handles zero cell counts by the limit 0 log 0 = 0.
 obs=np.array([depth-a1,depth-a2,a1,a2],float)
 pooled=(a1+a2)/2
 expected=np.array([depth-pooled,depth-pooled,pooled,pooled])
 g=2*np.sum(xlogy(obs,np.divide(obs,expected,out=np.ones_like(obs),where=expected>0)),axis=0)
 kernel=s['kernel'];gm=np.convolve(g,kernel,'same')/np.convolve(np.ones(len(g)),kernel,'same')
 starts,mid,l,r=s['grids'][2000000];pick=np.searchsorted(s['pos'],mid)
 result['G_tricube_h1Mb']=(gm[pick],starts,starts+2000000,r-l)
 # Significant-marker fraction using the SAME bulk-aware p99 cutoff; not Fisher/PyBSASeq.
 val,cnt=window_mean((abs(delta)>s['cuts'][1]).astype(float),np.ones(len(delta),bool),l,r,10)
 result['p99_fraction_2000kb']=(val,starts,starts+2000000,cnt)
 return result

def maxscore(entry):
 v=entry[0];return float(np.nanmax(v)) if np.any(np.isfinite(v)) else -np.inf

def evaluate(entry,cut,truth):
 v,st,en,cnt=entry;valid=np.isfinite(v)
 if not valid.any():return dict(detected=False,covered=False,total_width_mb=0.,top_width_mb=np.nan,peak_error_mb=np.nan,n_intervals=0,peak_n=0,tied_peaks=0)
 mx=np.nanmax(v);ties=np.flatnonzero(valid&np.isclose(v,mx,rtol=0,atol=1e-10));k=int(ties[0])
 peak=(st[k]+en[k])/2
 hit=np.flatnonzero(valid&(v>cut));intervals=[]
 for j in hit:
  if intervals and st[j]<=intervals[-1][1]:intervals[-1][1]=max(intervals[-1][1],float(en[j]))
  else:intervals.append([float(st[j]),float(en[j])])
 widths=[b-a for a,b in intervals];top=[b-a for a,b in intervals if a<=peak<=b]
 return dict(detected=bool(intervals),covered=any(a<=truth<=b for a,b in intervals),total_width_mb=sum(widths)/1e6,top_width_mb=(top[0]/1e6 if top else np.nan),peak_error_mb=abs(peak-truth)/1e6,n_intervals=len(intervals),peak_n=int(cnt[k]),tied_peaks=len(ties))

def wilson(k,n):
 if n==0:return (np.nan,np.nan)
 z=1.95996398454;p=k/n;den=1+z*z/n;mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return mid-half,mid+half

def main():
 global P
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=1000);ap.add_argument('--n-cal',type=int,default=2000);ap.add_argument('--scenario',choices=['all',*SCENARIOS],default='all');ap.add_argument('--output-dir',type=Path,default=P);ap.add_argument('--skip-plots',action='store_true');args=ap.parse_args()
 if args.n < 1 or args.n_cal < 19:ap.error('--n must be positive and --n-cal must be at least 19')
 P=args.output_dir.resolve()
 if (P/'benchmark_protocol.json').exists():ap.error('Choose a new output directory; existing benchmark results will not be overwritten.')
 for folder in [P,P/'tables',P/'figures',P/'qa']:folder.mkdir(parents=True,exist_ok=True)
 configs={k:{**BASE,**v} for k,v in SCENARIOS.items() if args.scenario in ['all',k]}
 config=dict(seed=20260914,n_calibration=args.n_cal,n_null_test=args.n,n_qtl_test=args.n,alpha=.05,length_bp=L,minimum_markers=10,configs=configs,interval='union of physical support of windows above held-out chromosome-max null threshold',cutoff='ceil((n_calibration+1)*0.95) order statistic; strict exceedance',limitations=['One chromosome, not genome-wide FWER','F2 additive single QTL only','Equal constant depth within each condition','Stylized genetic maps; no crossover interference','No read alignment or structural-variant simulation','Differential-error scenario assumes its artifact is represented in null calibration','Comparator scores are reimplementations, not official package executions'])
 (P/'benchmark_protocol.json').write_text(json.dumps(config,indent=2))
 t0=time.time();allrecords=[];cutrecords=[]
 for si,(name,c) in enumerate(configs.items()):
  s=setup(c);seeds=np.random.SeedSequence([20260914,list(SCENARIOS).index(name)]).spawn(3)
  cal_rng=np.random.default_rng(seeds[0]);null_rng=np.random.default_rng(seeds[1]);test_rng=np.random.default_rng(seeds[2]);maxima={}
  for i in range(args.n_cal):
   a1,a2,truth=simulate(c,s,cal_rng,False)
   for method,entry in score(a1,a2,c,s).items():maxima.setdefault(method,[]).append(maxscore(entry))
  rank=min(args.n_cal,math.ceil((args.n_cal+1)*.95));cuts={k:float(np.sort(v)[rank-1]) for k,v in maxima.items()}
  for method,cut in cuts.items():cutrecords.append(dict(scenario=name,method=method,threshold=cut,single_marker_p95=s['cuts'][0],single_marker_p99=s['cuts'][1]))
  for alt,rng,stage in [(False,null_rng,'null'),(True,test_rng,'qtl')]:
   for i in range(args.n):
    a1,a2,truth=simulate(c,s,rng,alt)
    for method,entry in score(a1,a2,c,s).items():
     allrecords.append(dict(scenario=name,stage=stage,replicate=i,method=method,truth_bp=truth if alt else np.nan,**evaluate(entry,cuts[method],truth)))
  print(f'{name}: completed {args.n_cal+2*args.n} populations; elapsed {time.time()-t0:.1f}s',flush=True)
  pd.DataFrame(allrecords).to_csv(P/'tables/benchmark_replicates.tsv.gz',sep='\t',index=False,compression='gzip')
 pd.DataFrame(cutrecords).to_csv(P/'tables/benchmark_thresholds.tsv',sep='\t',index=False)
 d=pd.DataFrame(allrecords);rows=[]
 for (sc,m),g in d.groupby(['scenario','method']):
  nu=g[g.stage=='null'];q=g[g.stage=='qtl'];n=len(q);fpr=nu.detected.mean();det=q.detected.mean();cov=q.covered.mean()
  fl,fu=wilson(int(nu.detected.sum()),len(nu));cl,cu=wilson(int(q.covered.sum()),n)
  rows.append(dict(scenario=sc,method=m,n=n,false_positive_rate=fpr,fpr_lower=fl,fpr_upper=fu,detection_rate=det,coverage=cov,coverage_lower=cl,coverage_upper=cu,coverage_given_detection=q.loc[q.detected,'covered'].mean(),median_peak_error_mb=q.peak_error_mb.median(),median_called_width_mb=q.loc[q.detected,'total_width_mb'].median(),median_total_width_all_mb=q.total_width_mb.median(),mean_total_width_all_mb=q.total_width_mb.mean(),median_peak_n=q.peak_n.median(),no_finite_peak_rate=q.peak_error_mb.isna().mean()))
 summary=pd.DataFrame(rows);summary.to_csv(P/'tables/benchmark_summary.tsv',sep='\t',index=False)
 pairs=[]
 for sc,g in d[d.stage=='qtl'].groupby('scenario'):
  for scale in [2000,500,100]:
   a=g[g.method==f'ALL_{scale}kb'].set_index('replicate');b=g[g.method==f'p99_{scale}kb'].set_index('replicate')
   delta=b.covered.astype(float)-a.covered.astype(float);se=delta.std(ddof=1)/np.sqrt(len(delta))
   err=b.peak_error_mb-a.peak_error_mb
   pairs.append(dict(scenario=sc,scale_kb=scale,coverage_difference=float(delta.mean()),paired_normal_lower=float(delta.mean()-1.96*se),paired_normal_upper=float(delta.mean()+1.96*se),median_paired_peak_error_difference_mb=float(err.median()),n_paired_finite_errors=int(err.notna().sum())))
 pd.DataFrame(pairs).to_csv(P/'tables/benchmark_paired_comparisons.tsv',sep='\t',index=False)
 if args.skip_plots:return
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'svg.fonttype':'none','font.family':'DejaVu Sans','font.size':9})
 methods=['ALL_2000kb','ALL_100kb','p99_2000kb','p99_100kb','G_tricube_h1Mb','p99_fraction_2000kb']
 colors=['#27658d','#73a4c3','#cf6239','#e4aa8b','#665191','#48845c']
 fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
 base=summary[summary.scenario=='baseline'].set_index('method')
 for ax,metric,ylabel in zip(axs,['false_positive_rate','coverage','median_called_width_mb'],['Null false-positive rate','True-position inclusion rate','Median called width (Mb)']):
  for j,m in enumerate(methods):
   value=base.loc[m,metric];ax.bar(j,value,color=colors[j])
   if metric in ['false_positive_rate','coverage']:
    prefix='fpr' if metric=='false_positive_rate' else 'coverage';lo=base.loc[m,prefix+'_lower'];hi=base.loc[m,prefix+'_upper'];ax.errorbar(j,value,yerr=[[value-lo],[hi-value]],color='black',capsize=2,lw=.8)
  ax.set_xticks(range(len(methods)),['All 2 Mb','All 0.1 Mb','p99 2 Mb','p99 0.1 Mb',"G tricube",'p99 fraction'],rotation=35,ha='right');ax.set_ylabel(ylabel)
  if metric=='false_positive_rate':ax.axhline(.05,ls='--',c='black',lw=.8)
 fig.savefig(P/'figures/Figure_4_independent_benchmark.svg');fig.savefig(P/'qa/Figure_4.png',dpi=140);plt.close(fig)
 fig,axs=plt.subplots(2,1,figsize=(12,7),layout='constrained')
 scs=list(configs)
 for ax,metric,ylabel in zip(axs,['coverage','median_peak_error_mb'],['True-position inclusion rate','Median peak error (Mb)']):
  for m,color in zip(methods,colors):ax.plot(range(len(scs)),[summary[(summary.scenario==sc)&(summary.method==m)][metric].iloc[0] for sc in scs],'o-',ms=3,label=m,color=color)
  ax.set_xticks(range(len(scs)),[s.replace('_',' ') for s in scs],rotation=30,ha='right');ax.set_ylabel(ylabel)
 axs[0].legend(ncol=3,fontsize=8)
 fig.savefig(P/'figures/Figure_5_sensitivity_benchmark.svg');fig.savefig(P/'qa/Figure_5.png',dpi=140);plt.close(fig)
 print(summary[summary.scenario=='baseline'].to_string(index=False),flush=True)

if __name__=='__main__':main()
