#!/usr/bin/env python3
"""Numerical checks for spatial diagnostics; execute through SLURM."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.ndimage import gaussian_filter

path = Path(__file__).resolve().parents[1]/'step2_spatial_calibration.py'
spec = importlib.util.spec_from_file_location('step2_spatial_calibration', path)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class SpatialCalibrationTests(unittest.TestCase):
    def test_fft_moments_equal_direct_masked_pairs_with_edges_and_missing_pixels(self):
        rng = np.random.default_rng(411)
        image = rng.normal(13, 2, (45, 53)) + np.arange(53)[None, :]/7
        mask = rng.random(image.shape) > .3
        mask[:7, :8] = False
        image[~mask] = np.nan
        got = m.pairwise_acf(image, mask, lag=8)
        for dy in range(-8, 9):
            for dx in range(-8, 9):
                want = m.direct_lag(image, mask, dx, dy)
                self.assertEqual(got['pairs'][dy+8, dx+8], want['pairs'])
                self.assertAlmostEqual(got['covariance'][dy+8, dx+8], want['covariance'], places=10)
                self.assertAlmostEqual(got['rho'][dy+8, dx+8], want['rho'], places=10)

    def test_correlation_invariant_to_offset_and_scale(self):
        rng = np.random.default_rng(812)
        image = rng.normal(size=(48, 64))
        mask = rng.random(image.shape) > .2
        a = m.pairwise_acf(image, mask, lag=7)
        b = m.pairwise_acf(17 + 3*image, mask, lag=7)
        np.testing.assert_allclose(a['rho'], b['rho'], atol=1e-12)
        np.testing.assert_allclose(9*a['covariance'], b['covariance'], atol=1e-12)

    def test_heldout_spatial_plane_recovers_injected_gradient(self):
        rng = np.random.default_rng(222)
        shape = (256, 288)
        truth = m.plane_image(shape, (0.17, 2., -.8))
        image = truth + rng.normal(0, .04, shape)
        labels = np.zeros(shape, dtype=np.int32)
        labels[90:125, 120:155] = 31
        image[labels != 0] += 50
        original = labels.copy()
        mask = m.sky_mask(image, labels, 10)
        cells = m.cell_statistics(image, mask, 32)
        got = m.plane_diagnostic(cells, rng, bootstrap=100)
        np.testing.assert_allclose(got['coefficients'], (.17, 2., -.8), atol=.005)
        self.assertLess(got['heldout_rmse_plane'], .04)
        self.assertGreater(got['heldout_rmse_constant'], .2)
        np.testing.assert_array_equal(labels, original)
        raw = m.pairwise_acf(image, mask, lag=16)
        detrended = m.pairwise_acf(image - m.plane_image(shape, got['coefficients']), mask, lag=16)
        self.assertGreater(m.acf_summary(raw)[2]['median_rho'], .9)
        self.assertLess(abs(m.acf_summary(detrended)[2]['median_rho']), .02)

    def test_aperture_variance_distinguishes_white_and_correlated_noise(self):
        rng = np.random.default_rng(926)
        shape = (256, 256)
        white = rng.normal(size=shape)
        correlated = gaussian_filter(rng.normal(size=shape), sigma=1)
        mask = np.ones(shape, dtype=bool)
        white_rows = m.aperture_statistics(white, mask, m.pairwise_acf(white, mask, lag=8),
                                           rng, sides=(1, 4, 8), bootstrap=100)
        corr_rows = m.aperture_statistics(correlated, mask, m.pairwise_acf(correlated, mask, lag=8),
                                          rng, sides=(1, 4, 8), bootstrap=100)
        self.assertAlmostEqual(white_rows[-1]['variance_ratio_to_independent'], 1., delta=.15)
        self.assertGreater(corr_rows[-1]['variance_ratio_to_independent'], 6.)
        for row in (white_rows[-1], corr_rows[-1]):
            self.assertAlmostEqual(row['variance']/row['covariance_sum_prediction_conditional_stationarity'],
                                   1., delta=.20)

    def test_apertures_require_all_pixels_blank(self):
        rng = np.random.default_rng(110)
        image = rng.normal(size=(128, 128))
        mask = np.ones(image.shape, dtype=bool)
        mask[::4, ::4] = False
        image[~mask] = 1e9
        got = m.aperture_statistics(image, mask, m.pairwise_acf(image, mask, 8),
                                    rng, sides=(4,), bootstrap=100)[0]
        self.assertEqual(got['status'], 'INSUFFICIENT_BLANK_APERTURES')
        self.assertEqual(got['n_apertures'], 0)

    def test_original_labels_and_radius_guard(self):
        image = np.ones((64, 64))
        with self.assertRaisesRegex(ValueError, 'integers'):
            m.sky_mask(image, np.zeros_like(image), 10)
        with self.assertRaisesRegex(ValueError, 'positive'):
            m.sky_mask(image, np.zeros_like(image, dtype=int), 0)
        with self.assertRaisesRegex(ValueError, 'insufficient'):
            m.sky_mask(image, np.ones_like(image, dtype=int), 10)

    def test_sparse_cells_do_not_prevent_raw_covariance_and_apertures(self):
        rng = np.random.default_rng(751)
        image = rng.normal(size=(128, 128))
        mask = np.zeros(image.shape, dtype=bool)
        mask[:48, :48] = True
        cells = m.cell_statistics(image, mask, 32)
        plane = m.plane_diagnostic(cells, rng, bootstrap=100)
        self.assertEqual(plane['status'], 'UNAVAILABLE_INSUFFICIENT_SPATIAL_COVERAGE')
        self.assertEqual(plane['ncells'], 3)
        self.assertEqual(plane['quadrant_coverage'][1]['ntest'], 0)
        acf = m.pairwise_acf(image, mask, lag=8)
        apertures = m.aperture_statistics(image, mask, acf, rng, sides=(4,), bootstrap=100)
        self.assertEqual(apertures[0]['n_apertures'], 144)
        self.assertGreater(apertures[0]['variance'], 0)
        self.assertNotIn('coefficients', plane)

    def test_empty_and_rank_deficient_plane_coverage_is_explicit(self):
        rng = np.random.default_rng(752)
        empty = m.plane_diagnostic([], rng)
        self.assertEqual(empty['design_rank'], 0)
        rows = [dict(accepted=True, xnorm=x, ynorm=0., median=x, fold=i % 4)
                for i, x in enumerate(np.linspace(-1, 1, 16))]
        result = m.plane_diagnostic(rows, rng)
        self.assertEqual(result['status'], 'UNAVAILABLE_INSUFFICIENT_SPATIAL_COVERAGE')
        self.assertEqual(result['design_rank'], 2)

    def test_missing_covariance_lags_do_not_become_zero_or_abort_json(self):
        rng = np.random.default_rng(753)
        image = rng.normal(size=(128, 128))
        mask = np.ones(image.shape, dtype=bool)
        acf = m.pairwise_acf(image, mask, lag=8)
        acf['covariance'][8, 9] = np.nan
        row = m.aperture_statistics(image, mask, acf, rng, sides=(4,), bootstrap=100)[0]
        self.assertIsNone(row['covariance_sum_prediction_conditional_stationarity'])
        self.assertEqual(row['covariance_prediction_missing_required_lags'], 1)
        self.assertEqual(row['covariance_prediction_status'],
                         'UNAVAILABLE_REQUIRED_LAGS_HAVE_INSUFFICIENT_PAIRS')
        self.assertGreater(row['variance'], 0)
        __import__('json').dumps(row, allow_nan=False)

    def test_covariance_prediction_requires_full_aperture_lag_support(self):
        rng = np.random.default_rng(754)
        image = rng.normal(size=(128, 128))
        mask = np.ones(image.shape, dtype=bool)
        acf = m.pairwise_acf(image, mask, lag=2)
        row = m.aperture_statistics(image, mask, acf, rng, sides=(4,), bootstrap=100)[0]
        self.assertEqual(row['covariance_prediction_status'], 'UNAVAILABLE_ACF_LAG_SUPPORT_TOO_SMALL')
        self.assertIsNone(row['covariance_sum_prediction_conditional_stationarity'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
