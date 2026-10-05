"""Reproduce 500 coordinate-assignment experiments, preserving seed 2025.

A half-window coordinate shift changes signed error, not its dispersion.
The first replicate supplies both illustrative scales in Figure 8.
"""
from pathlib import Path
import argparse, gzip
import numpy as np
import pandas as pd
from descriptive_core import make_chrom, windows, QTL_MB

SCALES=[2.,1.,.5,.1,.01]

def simulate(replicates=500,seed=2025,raw_path=None):
    rng=np.random.RandomState(seed);rows=[]
    raw=gzip.open(raw_path,'wt',encoding='utf-8',newline='') if raw_path else None
    try:
        for i in range(replicates):
            pos,val=make_chrom(rng)
            if raw:pd.DataFrame({'sim_i':i,'pos_bp':np.round(pos,1),'delta':np.round(val,6)}).to_csv(raw,sep='\t',index=False,header=i==0)
            keep=abs(val)>=np.percentile(abs(val),99)
            for mb in SCALES:
                width=int(mb*1e6);step=max(int(width*.05),1000)
                start,mid,means,_=windows(pos[keep],val[keep],width,step)
                if len(means)<3:continue
                k=int(np.argmax(means));es=start[k]/1e6-QTL_MB;ec=mid[k]/1e6-QTL_MB
                rows.append(dict(window_mb=mb,window_bp=width,error_start=es,error_center=ec,abs_err_start=abs(es),abs_err_center=abs(ec)))
    finally:
        if raw:raw.close()
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--replicates',type=int,default=500);ap.add_argument('--seed',type=int,default=2025)
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/descriptive'))
    ap.add_argument('--raw-snps',action='store_true',help='Also write the large raw SNP draws (Data S4).')
    args=ap.parse_args()
    if args.replicates<1:ap.error('--replicates must be positive')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    raw=args.output_dir/'Supplementary_Data_S4_snp_data.tsv.gz' if args.raw_snps else None
    df=simulate(args.replicates,args.seed,raw)
    target=args.output_dir/'Supplementary_Data_S2_peakshift_sim.tsv';df.to_csv(target,sep='\t',index=False)
    print(f'Wrote {len(df):,} coordinate result rows to {target}')

if __name__=='__main__':main()
