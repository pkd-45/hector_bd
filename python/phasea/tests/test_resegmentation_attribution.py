"""Re-segmentation attribution tests for D1 (Claude review of V02, 31 Oct 2026).

Scenario: production passed segim=NULL, so ProFound re-segmented the image. Its own map
labels a tight core 31 (the fitted target) and splits part of the outskirts into a separate
label 32; the rest of the WAVES segment is background in ProFound's map.

Rule under test: a missing target pixel was excluded by segmentation if, in the segmentation
Found2Fit actually used, it does NOT carry the label(s) of the fitted region. Mapping saved
labels to their majority WAVES label is the wrong frame: label 32 lies mostly inside the WAVES
target, so a majority map sends it back to the target and hides the exclusion.

V02 fails both tests: 832 of 1680 pixels 'unexplained' with profound saved; all 1680 without it,
although Data$segim shows the excluding label for every missing pixel inside the stamp.
"""
import csv, io, sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from export_phasea import best_offset, compare_family

N = 220; yy, xx = np.mgrid[0:N, 0:N].astype(float)
def sersic(yc, xc, re, n, q, pa, amp):
    t = np.deg2rad(pa); dx, dy = xx - xc, yy - yc
    a = dx*np.cos(t) + dy*np.sin(t); b = (-dx*np.sin(t) + dy*np.cos(t))/q
    return amp*np.exp(-(2*n - 1/3)*((np.hypot(a, b)/re)**(1/n) - 1))
GAL = sersic(110, 105, 14, 1.5, 0.6, 30, 40); IMG = GAL + np.random.default_rng(7).normal(0, 1, (N, N))
L = np.zeros((N, N), int); L[GAL > 2.0] = 5                                   # WAVES target, label 5
P = np.zeros((N, N), int); P[GAL > 6.0] = 31; P[(GAL > 2.0) & (GAL <= 6.0) & (xx > 105)] = 32
REGION = (P == 31)
YS, XS = np.nonzero(REGION); SL = (slice(YS.min(), YS.max()+1), slice(XS.min(), XS.max()+1))

def run(with_profound):
    arr = {('fit_FS', 'image'):  np.asfortranarray(IMG[SL].T),
           ('fit_FS', 'region'): np.asfortranarray(REGION[SL].astype(float).T),
           ('fit_FS', 'segim'):  np.asfortranarray(P[SL].astype(float).T)}
    m = best_offset(IMG, arr[('fit_FS', 'image')])
    kw = dict(original_image=IMG)
    if with_profound: kw.update(passed_image=IMG, passed_seg=P, passed_reg=(0, 0, N, N, 'direct'))
    return compare_family(L, [5], arr[('fit_FS', 'region')], m, arr, 'fit_FS', None, None,
                          csv.writer(io.StringIO()), 'T', **kw)

class TestResegmentationAttribution(unittest.TestCase):
    def test_profound_saved_attributes_every_exclusion(self):
        r = run(True)
        self.assertEqual(r['N_missing'], int((L == 5).sum() - REGION.sum()))
        self.assertEqual(r['class_counts'], {'segim_differs_upstream': r['N_missing']},
            'with the full passed map known, every missing target pixel carries a non-target passed label')

    def test_profound_stripped_uses_data_segim_inside_stamp(self):
        r = run(False)
        missing = (L == 5) & ~REGION
        inside = np.zeros_like(missing); inside[SL] = True
        n_in, n_out = int((missing & inside).sum()), int((missing & ~inside).sum())
        self.assertEqual(r['class_counts'].get('segim_differs_upstream', 0), n_in,
            'inside the stamp, Data$segim is the passed map cropped and shows the excluding label')
        self.assertEqual(r['class_counts'].get('unexplained', 0), n_out,
            'outside the stamp with no parent domain, the cause must stay unresolved')

if __name__ == '__main__':
    unittest.main(verbosity=2)
