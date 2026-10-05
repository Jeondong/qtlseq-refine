"""Reproduce 600 descriptive observed-percentile experiments (seed 42).

gw_noise = SD of retained window means; discriminability = mean/SD.
Neither is an independently measured background-noise or localization metric.
"""
from pathlib import Path
import argparse, gzip
import numpy as np
import pandas as pd
from descriptive_core import make_chrom, windows, WINDOWS, QTL_MB

def simulate(replicates=600, seed=42, raw_path=None):
    rng=np.random.RandomState(seed); rows=[]
    raw=gzip.open(raw_path,'wt',encoding='utf-8',newline='') if raw_path else None
    try:
        for i in range(replicates):
            pos,val=make_chrom(rng)
            if raw:
                pd.DataFrame({'sim_i':i,'pos_bp':np.round(pos,1),'delta':np.round(val,6)}).to_csv(raw,sep='\t',index=False,header=i==0)
            for label,(width,step) in WINDOWS.items():
                for name,percentile in [('All SNPs',0),('p95',95),('p99',99)]:
                    cut=0 if percentile==0 else np.percentile(abs(val),percentile)
                    keep=abs(val)>=cut
                    _,mid,means,_=windows(pos[keep],val[keep],width,step)
                    if len(means)<5:continue
                    region=abs(mid/1e6-QTL_MB)<=2.5
                    signal=means[region].mean() if region.sum()>=2 else np.nan
                    sd=means.std();pk=int(np.argmax(means));error=abs(mid[pk]/1e6-QTL_MB)
                    above=np.flatnonzero(means>=means[pk]/2)
                    span=(mid[above[-1]]-mid[above[0]])/1e6 if len(above)>=2 else np.nan
                    rows.append(dict(window=label,threshold=name,n_snps=int(keep.sum()),qtl_signal=signal,
                                     gw_noise=sd,discriminability=signal/sd if sd>0 else np.nan,
                                     fwhm=span,peak_err=error,detected=int(error<=2.5)))
    finally:
        if raw:raw.close()
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--replicates',type=int,default=600);ap.add_argument('--seed',type=int,default=42)
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/descriptive'))
    ap.add_argument('--raw-snps',action='store_true',help='Also write the large raw SNP draws (Data S3).')
    args=ap.parse_args()
    if args.replicates<1:ap.error('--replicates must be positive')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    raw=args.output_dir/'Supplementary_Data_S3_snp_data.tsv.gz' if args.raw_snps else None
    df=simulate(args.replicates,args.seed,raw)
    target=args.output_dir/'Supplementary_Data_S1_threshold_sim.tsv';df.to_csv(target,sep='\t',index=False)
    print(f'Wrote {len(df):,} descriptive result rows to {target}')

if __name__=='__main__':main()
