from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from descriptive_core import windows,make_chrom
from peak_shift_sim import simulate
from population_benchmark import evaluate

class AnalysisTests(unittest.TestCase):
    def test_half_open_support_and_missing_windows(self):
        st,mid,mean,count=windows(np.array([0,10,20]),np.array([1.,2.,3.]),10,10,minimum=1,stop=40,include_empty=True)
        np.testing.assert_array_equal(count,[1,1,1,0]);np.testing.assert_allclose(mean[:3],[1,2,3]);self.assertTrue(np.isnan(mean[3]))
    def test_fast_means_match_direct_means(self):
        pos,val=make_chrom(np.random.RandomState(9));st,_,means,counts=windows(pos,val,100000,10000)
        for i in [0,len(st)//2,len(st)-1]:
            keep=(pos>=st[i])&(pos<st[i]+100000);self.assertEqual(counts[i],keep.sum());self.assertEqual(means[i],val[keep].mean())
    def test_first_example_and_coordinate_identity(self):
        frame=simulate(1);np.testing.assert_allclose(frame.error_center-frame.error_start,frame.window_mb/2,atol=1e-12)
        self.assertAlmostEqual(frame.loc[frame.window_mb==.1,'error_start'].iloc[0],1.38)
    def test_f2_union_and_unsupported_profile(self):
        result=evaluate((np.array([2.,2.,0.]),np.array([0.,5.,20.]),np.array([10.,15.,30.]),np.array([10,10,10])),1.,12.)
        self.assertTrue(result['covered']);self.assertEqual(result['total_width_mb'],15/1e6)
        empty=evaluate((np.array([np.nan]),np.array([0.]),np.array([10.]),np.array([0])),-np.inf,5.)
        self.assertFalse(empty['detected']);self.assertFalse(empty['covered']);self.assertTrue(np.isnan(empty['peak_error_mb']))

if __name__=='__main__':unittest.main()
