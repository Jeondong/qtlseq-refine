"""Reaggregate saved F2 and empirical results without simulating new populations."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from sliding_window_analysis import scan
from assemble_raw_data import assemble_archives

ROOT=Path(__file__).resolve().parents[1]
DATASETS=[('abci8_Chr11','Chr11','OsABCI8',17329608,17332992,'min'),
          ('qHT1_Chr1','Chr1','qHT1',4665786,4667096,'min'),
          ('lesv_Chr11','Chr11','OsLESV',22175496,22178936,'max')]

def verify_f2(root):
    base=root/'data/supplementary/Supplementary_Data_S5';raw=pd.read_csv(base/'benchmark_replicates.tsv.gz',sep='\t',keep_default_na=False,na_values=[''])
    saved=pd.read_csv(base/'benchmark_summary.tsv',sep='\t').set_index(['scenario','method'])
    selected=pd.read_csv(base/'benchmark_verified_summary.tsv',sep='\t').set_index(['scenario','method'])
    assert len(raw)==364000
    assert len(raw[['scenario','stage','replicate']].drop_duplicates())==26000
    assert not raw.duplicated(['scenario','stage','replicate','method']).any()
    checks=0
    for key,g in raw.groupby(['scenario','method']):
        q=g[g.stage=='qtl'];null=g[g.stage=='null']
        assert len(q)==len(null)==1000
        vals=dict(coverage=q.covered.mean(),detection_rate=q.detected.mean(),false_positive_rate=null.detected.mean(),
                  coverage_given_detection=q.loc[q.detected,'covered'].mean(),median_peak_error_mb=q.peak_error_mb.median(),
                  median_called_width_mb=q.loc[q.detected,'total_width_mb'].median(),no_finite_peak_rate=q.peak_error_mb.isna().mean())
        for col,v in vals.items():
            assert np.isclose(v,saved.loc[key,col],rtol=1e-10,atol=1e-10,equal_nan=True),(key,col)
            if key in selected.index:assert np.isclose(v,selected.loc[key,col],rtol=1e-10,atol=1e-10,equal_nan=True),(key,col)
            checks+=1
    return dict(replicate_method_rows=len(raw),evaluation_populations=26000,calibration_populations=26000,summary_values_verified=checks,new_f2_populations_generated=0)

def summarize_empirical(root):
    rows=[];nchecked=0;supported=[]
    for stem,chrom,label,g1,g2,direction in DATASETS:
        inp=root/'data'/f'{stem}_p99.xlsx';data=pd.read_excel(inp).iloc[:,:2].dropna()
        regenerated=scan(inp,chrom);saved=pd.read_excel(root/'data'/f'{stem}_p99_sliding_window_multi.xlsx',sheet_name=None)
        for name,calc in regenerated.items():
            allrows=saved[name];valid=allrows[(allrows['count']>0)&allrows.value_mean.notna()].sort_values('start').reset_index(drop=True)
            assert np.array_equal(calc[['start','end','mid','count']],valid[['start','end','mid','count']]),(label,name,'support')
            assert np.allclose(calc.value_mean,valid.value_mean,rtol=0,atol=1e-10),(label,name,'means')
            export=calc.copy();export.insert(0,'window',name);export.insert(0,'dataset',label);supported.append(export)
            extreme=valid.value_mean.min() if direction=='min' else valid.value_mean.max()
            tied=valid[valid.value_mean==extreme];peak=tied.iloc[0]
            dist=lambda x:max(g1-x,0,x-g2)/1e6
            rows.append(dict(dataset=label,window=name,n_input=len(data),input_min_bp=int(data.iloc[:,0].min()),input_max_bp=int(data.iloc[:,0].max()),
                             supported_windows=len(valid),empty_rows=int((allrows['count']==0).sum()),median_count=float(valid['count'].median()),
                             singleton_percent=float((valid['count']==1).mean()*100),peak_direction=direction,extreme_mean=float(extreme),tied_windows=len(tied),
                             peak_start_bp=int(peak.start),peak_center_bp=int(peak.mid),peak_snp_count=int(peak['count']),gene_start_bp=g1,gene_end_bp=g2,
                             start_to_gene_mb=dist(peak.start),center_to_gene_mb=dist(peak.mid)))
            nchecked+=len(valid)
    return pd.DataFrame(rows),nchecked,pd.concat(supported,ignore_index=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/verification'))
    ap.add_argument('--verify-raw',action='store_true',help='Also verify the compressed and decompressed Data S3/S4 hashes and record counts.')
    args=ap.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    report=verify_f2(args.root);emp,n,supported=summarize_empirical(args.root)
    expected=pd.read_csv(args.root/'data/supplementary/Supplementary_Data_S6/empirical_summary.tsv',sep='\t')
    pd.testing.assert_frame_equal(emp,expected,check_exact=False,rtol=1e-10,atol=1e-10)
    expected_windows=pd.read_csv(args.root/'data/supplementary/Supplementary_Data_S6/empirical_supported_windows.tsv',sep='\t')
    pd.testing.assert_frame_equal(supported,expected_windows,check_exact=False,rtol=1e-10,atol=1e-10)
    assert supported['count'].gt(0).all() and supported.value_mean.notna().all()
    supported.to_csv(args.output_dir/'empirical_supported_windows.tsv',sep='\t',index=False)
    emp.to_csv(args.output_dir/'empirical_summary.tsv',sep='\t',index=False)
    report.update(empirical_supported_windows_verified=n,empirical_settings_verified=len(emp),empirical_summary_matches=True,empirical_supported_export_matches=True)
    if args.verify_raw:report.update(assemble_archives(args.root,args.output_dir/'raw_archives'))
    (args.output_dir/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
