#!/usr/bin/env python3
"""Synthetic in-memory guard tests; no historical file access, no fitting."""
import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import step1_segmap_audit as a
import wcs_registration_blockxcorr as c

class TranslationWCS:
    def __init__(self,dx=0.,dy=0.):self.dx=float(dx);self.dy=float(dy)
    def all_pix2world(self,pts,origin):return np.asarray(pts)
    def all_world2pix(self,pts,origin):return np.asarray(pts)+[self.dx,self.dy]

class SyntheticTests(unittest.TestCase):
    def test_orientation_never_uses_blank_sky_dominance(self):
        L=np.zeros((48,48),np.int32)
        L[5:13,27:36]=35;L[34:43,8:14]=2
        orig,way,score,amb=a.orient_saved(L.T.copy(),L)
        self.assertTrue(np.array_equal(orig,L))
        self.assertEqual(way,'transpose');self.assertEqual(score,1.)
        self.assertFalse(amb)

    def test_frozen_roster_832_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'roster.tsv'
            p.write_text('id\n'+''.join(f'{i:014d}\n' for i in range(1,833)))
            ids,col=a.roster_ids(p)
            self.assertEqual((len(ids),col),(832,'id'))
            p.write_text('id\n'+''.join(f'{i:014d}\n' for i in range(1,832)))
            with self.assertRaisesRegex(ValueError,'FROZEN_R61_832_GATE_FAIL'):
                a.roster_ids(p)

    def test_registration_sign_and_candidate_verdicts(self):
        rng=np.random.default_rng(2)
        yy,xx=np.ogrid[:512,:512]
        image=np.zeros((512,512),float);seg=np.zeros((512,512),np.int32)
        for i in range(150):
            cy,cx=rng.integers(25,487,2)
            rr=(yy-cy)**2+(xx-cx)**2
            image+=35*np.exp(-rr/15.0)
            seg[rr<19]=i+1
        image+=rng.normal(0,0.2,image.shape)
        moved=np.roll(np.roll(seg,-7,axis=0),10,axis=1)
        result=c.analyse(image,moved,TranslationWCS(),TranslationWCS(10,-7),nblock=4,search=15)
        self.assertEqual(result['verdict'],'MISREGISTERED_CANDIDATE')
        self.assertGreaterEqual(result['n_usable'],8)
        self.assertLess(result['median_residual_from_header_px'],.35)
        result=c.analyse(image,seg,TranslationWCS(),TranslationWCS(10,-7),nblock=4,search=15)
        self.assertEqual(result['verdict'],'REGISTERED_CANDIDATE')
        result=c.analyse(image,seg,TranslationWCS(),TranslationWCS(),nblock=4,search=15)
        self.assertEqual(result['verdict'],'HEADERS_AGREE_NO_DISCRIMINATION')

if __name__=='__main__':unittest.main(verbosity=2)
