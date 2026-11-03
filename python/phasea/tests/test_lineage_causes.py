"""Adversarial lineage tests for the D1 exporter (Claude review, 31 Oct 2026).

Simulates ProFuse 2a88cb7 profuseFound2Fit: optional cutbox -> region = target label & mask==0
-> tightcrop to the region's bounding box. Arrays are delivered in R/Rfits [x,y] layout, as the
R extractor writes them. Each scenario applies ONE known historical cause.

The registration and counting assertions pass on V01. The cause assertions for the island-filter
and autoclip scenarios FAIL on V01 by design: V01 labels every pixel outside the saved stamp as
'cutbox_side_fit_crop', but under tightcrop the stamp is the bounding box of whatever survived,
so upstream removals also fall outside it.
"""
import csv, io, sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from export_phasea import best_offset, compare_family

N = 220
yy, xx = np.mgrid[0:N, 0:N].astype(float)

def sersic(yc, xc, re, n, q, pa, amp):
    t = np.deg2rad(pa); dx, dy = xx - xc, yy - yc
    a = dx*np.cos(t) + dy*np.sin(t); b = (-dx*np.sin(t) + dy*np.cos(t))/q
    return amp*np.exp(-(2*n - 1/3)*((np.hypot(a, b)/re)**(1/n) - 1))

GAL = sersic(110, 105, 14, 1.5, 0.6, 30, 40)
IMG = GAL + sersic(60, 170, 3, 1, 0.9, 0, 25) + np.random.default_rng(42).normal(0, 1, (N, N))
L = np.zeros((N, N), int); L[GAL > 2.0] = 5; L[sersic(60, 170, 3, 1, 0.9, 0, 25) > 2.0] = 9
ISLAND = (slice(150, 154), slice(150, 155)); L[ISLAND] = 5   # detached island of the target label
TARGET = [5]

def found2fit_like(segim, cutbox=None, autoclip_px=()):
    im = IMG.copy(); seg = segim.copy()
    for y, x in autoclip_px: im[y, x] = np.nan
    oy = ox = 0
    if cutbox is not None:
        (y0, y1), (x0, x1) = cutbox; im, seg = im[y0:y1, x0:x1], seg[y0:y1, x0:x1]; oy, ox = y0, x0
    mask = ~np.isfinite(im); region = np.isin(seg, TARGET) & ~mask
    ys, xs = np.nonzero(region); sl = (slice(ys.min(), ys.max()+1), slice(xs.min(), xs.max()+1))
    data = {'image': im[sl], 'region': region[sl].astype(float), 'mask': mask[sl].astype(float),
            'sigma': np.ones(im[sl].shape), 'segim': seg[sl].astype(float)}
    r = {k: np.asfortranarray(v.T) for k, v in data.items()}          # R/Rfits [x,y] layout
    passed = {'image': im, 'segim': seg, 'offset': (oy, ox)}          # what Found2Fit received
    return r, (oy + ys.min(), ox + xs.min()), passed

def audit(segim, **kw):
    arrs, true_off, passed = found2fit_like(segim, **kw)
    fam = {('fit_FS', k): v for k, v in arrs.items()}
    fam[('inputs', 'image')] = np.asfortranarray(passed['image'].T)   # image given to Found2Fit, R layout
    fam[('inputs', 'segim')] = np.asfortranarray(passed['segim'].T)
    oy, ox = passed['offset']; h, w = passed['segim'].shape
    m = best_offset(IMG, arrs['image'])
    res = compare_family(L, TARGET, arrs['region'], m, fam, 'fit_FS', passed['segim'], (oy, ox, h, w),
                         csv.writer(io.StringIO()), 'T')
    return m, res, true_off

class TestLineageCauses(unittest.TestCase):
    def test_static_nothing_removed(self):
        m, res, off = audit(L)
        self.assertEqual(m['status'], 'VERIFIED_PIXEL_REGISTRATION')
        self.assertEqual((res['offset_row0'], res['offset_col0']), tuple(map(int, off)))
        self.assertEqual(res['N_missing'], 0); self.assertEqual(res['N_extra'], 0)

    def test_island_filter_is_not_a_cut(self):
        L1 = L.copy(); L1[ISLAND] = 0
        m, res, off = audit(L1)
        self.assertEqual((res['offset_row0'], res['offset_col0']), tuple(map(int, off)))
        self.assertEqual(res['N_missing'], 20)
        self.assertNotIn('cutbox_side_fit_crop', res['class_counts'],
            'island removed upstream lies inside the image passed to Found2Fit; it is not a crop')

    def test_real_cutbox_is_a_cut(self):
        m, res, off = audit(L, cutbox=((100, 220), (0, 220)))
        self.assertEqual((res['offset_row0'], res['offset_col0']), tuple(map(int, off)))
        self.assertGreater(res['N_missing'], 0)
        self.assertEqual(set(res['class_counts']), {'cutbox_side_fit_crop'})

    def test_autoclip_is_not_a_cut(self):
        tys, txs = np.nonzero(np.isin(L, TARGET)); left = txs.min()
        px = [(y, x) for y, x in zip(tys, txs) if x == left]
        m, res, off = audit(L, autoclip_px=px)
        self.assertEqual(res['N_missing'], len(px))
        self.assertNotIn('cutbox_side_fit_crop', res['class_counts'],
            'autoclip NA pixels lie inside the image passed to Found2Fit; they are not a crop')

if __name__ == '__main__':
    unittest.main(verbosity=2)
