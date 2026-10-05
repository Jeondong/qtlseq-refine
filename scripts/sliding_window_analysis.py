"""Independent window scans of a supplied p99-filtered empirical input.

Selection is performed upstream. This script does not test marker significance.
Each window is [start, end); unsupported windows are omitted.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from descriptive_core import windows

CONFIGS=[('2Mb_100Kb',2000000,100000),('1Mb_50Kb',1000000,50000),
         ('0.5Mb_10Kb',500000,10000),('0.1Mb_5Kb',100000,5000),('0.01Mb_5Kb',10000,5000)]

def scan(input_file,chrom,minimum=1):
    df=pd.read_excel(input_file).iloc[:,:2].dropna()
    data=df.to_numpy(dtype=float)
    if not len(data) or not np.isfinite(data).all():raise ValueError('Input must contain finite positions and values.')
    if (data[:,0]<0).any():raise ValueError('Positions must be nonnegative.')
    order=np.argsort(data[:,0]);pos=data[order,0];val=data[order,1];result={}
    for name,width,step in CONFIGS:
        st,mid,means,count=windows(pos,val,width,step,minimum,int(pos.max())+1)
        result[name]=pd.DataFrame({'chrom':chrom,'start':st,'end':st+width,'mid':mid.astype(int),'value_mean':means,'count':count})
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True);ap.add_argument('--chrom',required=True)
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/empirical'))
    ap.add_argument('--min-snps',type=int,default=1)
    args=ap.parse_args()
    if args.min_snps<1:ap.error('--min-snps must be positive')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    for name,df in scan(args.input,args.chrom,args.min_snps).items():
        target=args.output_dir/f'{args.input.stem}_{name}.tsv';df.to_csv(target,sep='\t',index=False)
        print(f'{name}: {len(df):,} supported windows -> {target}')

if __name__=='__main__':main()
