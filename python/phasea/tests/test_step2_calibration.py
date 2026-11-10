#!/usr/bin/env python3
"""Small synthetic tests; safe on a login node (no FITS I/O)."""
import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

modpath = Path(__file__).resolve().parents[1] / 'step2_calibration.py'
if not modpath.is_file():  # local authoring layout; deployed file lives in tests/
    modpath = Path(__file__).with_name('step2_calibration.py')
spec = importlib.util.spec_from_file_location('step2_calibration', modpath)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SkyTests(unittest.TestCase):
    def test_sky_recovers_pedestal_and_sigma_and_excludes_objects(self):
        rng = np.random.default_rng(441)
        image = rng.normal(0.17, 0.023, (176, 192))
        labels = np.zeros_like(image, dtype=np.int32)
        labels[65:91, 74:104] = 999
        image[labels > 0] += 7
        blank, got = module.compute_sky(image, labels)
        self.assertGreater(got['n_sky'], 20000)
        self.assertAlmostEqual(got['pedestal'], 0.17, delta=0.001)
        self.assertAlmostEqual(got['sigma_mad'], 0.023, delta=0.001)
        self.assertTrue(np.all(blank[labels != 0] == 0))
        self.assertAlmostEqual(got['z_std'], 1., delta=0.06)

    def test_acf_white_noise_near_zero(self):
        rng = np.random.default_rng(200)
        img = rng.normal(0, 1, (250, 240))
        blank = np.ones(img.shape, dtype=bool)
        ac = module.sky_acf(img, blank, 0, 2)
        self.assertTrue(all(abs(r['rho']) < 0.04 for r in ac
                            if r['rho'] is not None and (r['dx'] or r['dy'])))

    def test_blank_sky_insufficient_raises(self):
        with self.assertRaisesRegex(ValueError, 'insufficient'):
            module.compute_sky(np.ones((30,30)), np.ones((30,30), dtype=int))


if __name__ == '__main__':
    unittest.main(verbosity=2)
