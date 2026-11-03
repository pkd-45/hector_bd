"""New synthetic byte-level safeguards; these are NOT Friday observations."""
import csv
import io
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from export_phasea import compare_family,best_offset,partition_comparison,unique_seg_offset

class TestProFoundInputs(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(194)
        self.im=rng.normal(size=(67,81))
        self.L=np.zeros(self.im.shape,dtype=int)
        self.L[12:36,18:44]=7
        self.L[53:57,60:63]=7
        self.L[42:47,7:11]=3

    def test_partition_relabel_invariant_and_background_not_relabelled(self):
        p=self.L.copy()
        p[p==7]=99;p[p==3]=21
        a=partition_comparison(self.L,p)
        self.assertEqual(a['status'],'PARTITION_IDENTICAL_RELABELED')
        self.assertEqual(a['mismatch_pixels'],0)
        self.assertEqual(a['map_new_to_orig'][99],7)
        p[54,61]=0
        a=partition_comparison(self.L,p)
        self.assertEqual(a['status'],'PARTITION_DIFFERS')
        self.assertEqual(a['mismatch_pixels'],1)

    def test_saved_profound_domain_prevents_island_misattribution(self):
        # ProFound parent unmodified; Data comes from a downstream island filter.
        fitted=self.L.copy();fitted[53:57,60:63]=0
        stamp=self.im[12:36,18:44]
        reg=np.isin(fitted[12:36,18:44],[7]).astype(float)
        arrays={('fit_FS','image'):np.asfortranarray(stamp.T),
                ('fit_FS','region'):np.asfortranarray(reg.T),
                ('fit_FS','segim'):np.asfortranarray(fitted[12:36,18:44].T)}
        match=best_offset(self.im,arrays[('fit_FS','image')])
        writer=csv.writer(io.StringIO())
        res=compare_family(self.L,[7],arrays[('fit_FS','region')],match,arrays,'fit_FS',None,None,writer,'T',
                original_image=self.im,passed_image=self.im,passed_seg=fitted,
                passed_reg=(0,0,*self.im.shape,'direct'))
        self.assertEqual(res['N_missing'],12)
        self.assertEqual(res['class_counts'],{'segim_differs_upstream':12})
        self.assertEqual(res['N_extra'],0)
        self.assertTrue(res['stamp_exact_copy'])
        self.assertEqual(res['stamp_max_abs_original_delta'],0.)
        self.assertEqual(res['profound_seg_partition'],'PARTITION_DIFFERS')

    def test_saved_profound_autoclip_outside_stamp(self):
        # ProFound image contains finite-source -> NaN; tightcrop then excludes it.
        pimg=self.im.copy()
        pimg[53:57,60:63]=np.nan
        stamp=self.im[12:36,18:44]
        arrays={('fit_FS','image'):np.asfortranarray(stamp.T),
                ('fit_FS','region'):np.asfortranarray(np.ones(stamp.T.shape))}
        m=best_offset(self.im,arrays[('fit_FS','image')]);w=csv.writer(io.StringIO())
        res=compare_family(self.L,[7],arrays[('fit_FS','region')],m,arrays,'fit_FS',None,None,w,'T',
                           original_image=self.im,passed_image=pimg,passed_seg=self.L,
                           passed_reg=(0,0,*self.im.shape,'direct'))
        self.assertEqual(res['class_counts'],{'autoclip_NA':12})

    def test_unique_exact_segmap_window_transposed(self):
        # Distinctive non-rectangular original segments: no orientation guess.
        cut=self.L[8:40,11:50]
        off=unique_seg_offset(self.L,np.asfortranarray(cut.T))
        self.assertEqual(off,(8,11,'transpose'))

    def test_unknown_parent_domain_not_cropped(self):
        cut=self.im[12:36,18:44]
        arr={('fit_FS','image'):cut,('fit_FS','region'):np.ones(cut.shape)}
        w=csv.writer(io.StringIO());m=best_offset(self.im,cut)
        res=compare_family(self.L,[7],arr[('fit_FS','region')],m,arr,'fit_FS',None,None,w,'T',original_image=self.im)
        self.assertEqual(res['N_missing'],12)
        self.assertEqual(res['class_counts'],{'unexplained':12})

if __name__=='__main__':unittest.main()
