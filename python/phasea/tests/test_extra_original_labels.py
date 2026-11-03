import csv, io, sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from export_phasea import compare_family, best_offset

class TestOriginalExtraLabels(unittest.TestCase):
    def test_background_and_neighbour_extra_separately(self):
        rng=np.random.default_rng(221)
        image=rng.normal(size=(50,50))
        L=np.zeros((50,50),dtype=int)
        L[10:14,10:14]=8  # target
        L[14:17,12:15]=9  # neighbour
        row,col=8,8
        stamp=image[row:row+14,col:col+14]
        reg=np.zeros(stamp.shape,float)
        reg[10-row:14-row,10-col:14-col]=1
        reg[9-row,10-col]=1  # original background extra
        reg[14-row,12-col]=1  # original neighbour extra
        cap=io.StringIO();w=csv.writer(cap)
        fit={'fit_FS':None}
        arrays={('fit_FS','image'):stamp,('fit_FS','region'):reg}
        res=compare_family(L,[8],reg,best_offset(image,stamp),arrays,'fit_FS',None,None,w,'g')
        self.assertEqual(res['N_target'],16)
        self.assertEqual(res['N_missing'],0)
        self.assertEqual(res['N_extra'],2)
        self.assertEqual(res['N_extra_original_background'],1)
        self.assertEqual(res['N_extra_original_neighbour'],1)
        rows=list(csv.reader(io.StringIO(cap.getvalue())))
        self.assertEqual({r[5] for r in rows},{'extra_original_background','extra_original_neighbour'})
        self.assertTrue(any('DILATION_CONSISTENT_NOT_CAUSAL_PROOF' in r[6] for r in rows))
        self.assertTrue(any('WAVES_LABEL_9' in r[6] for r in rows))
    def test_zero_extra(self):
        rng=np.random.default_rng(222);image=rng.normal(size=(30,30));L=np.zeros((30,30),int);L[12:16,12:16]=8
        stamp=image[9:20,9:20];reg=(L[9:20,9:20]==8).astype(float)
        z=compare_family(L,[8],reg,best_offset(image,stamp),{('fit_FS','image'):stamp,('fit_FS','region'):reg},'fit_FS',None,None,csv.writer(io.StringIO()),'g')
        self.assertEqual(z['N_extra_original_background'],0)
        self.assertEqual(z['N_extra_original_neighbour'],0)
if __name__=='__main__':unittest.main()
