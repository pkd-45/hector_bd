import csv, io, json, sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from export_phasea import best_offset, compare_family, robust_verification, PILOTS

class TestPixelAudit(unittest.TestCase):
    def setUp(self):
        rg=np.random.default_rng(811)
        self.image=rg.normal(0,1,(64,64))
        self.L=np.zeros((64,64),dtype=int)
        self.L[21:29,25:34]=7
        self.L[45,47]=7 # disconnected island must be included
    def test_unique_affine_sky_registration(self):
        cut=self.image[19:39,22:41]*1.35-0.75
        reg=best_offset(self.image,cut)
        self.assertEqual(reg['status'],'VERIFIED_PIXEL_REGISTRATION')
        self.assertEqual(reg['best'][1:4],(19,22,'direct'))
    def test_transposed_registration(self):
        cut=(self.image[19:39,22:41]*1.35-0.75).T
        reg=best_offset(self.image,cut)
        self.assertEqual(reg['status'],'VERIFIED_PIXEL_REGISTRATION')
        self.assertEqual(reg['best'][1:4],(19,22,'transpose'))
    def test_unobserved_island_never_attributed_to_cut_and_extra_counted(self):
        im=self.image[19:39,22:41]
        region=np.isin(self.L[19:39,22:41],[7]).astype(float)
        region[0,0]=1 # extra background pixel
        m=best_offset(self.image,im)
        arr={('fit_FS','image'):im,('fit_FS','region'):region,
             ('fit_FS','sigma'):np.ones(region.shape),('fit_FS','mask'):np.zeros(region.shape),
             ('fit_FS','calcregion'):np.ones(region.shape)}
        capture=io.StringIO();writer=csv.writer(capture)
        d=compare_family(self.L,[7],region,m,arr,'fit_FS',None,None,writer,'1')
        self.assertEqual(d['N_target'],73)
        self.assertEqual(d['N_missing'],1)
        self.assertEqual(d['N_extra'],1)
        self.assertEqual(d['class_counts']['unexplained'],1)
        self.assertEqual(d['N_calc_only'],int(region.size-region.sum()))
        self.assertIn('45,47,unexplained',capture.getvalue())
        self.assertFalse(d['matches_target_label'])
    def test_nonunique_region_rejected(self):
        flat=np.zeros((64,64))
        self.assertIsNone(best_offset(flat,flat[1:10,1:10]))
    def test_pilot_n(self):
        self.assertEqual(len(PILOTS),10)

if __name__=='__main__':unittest.main()
