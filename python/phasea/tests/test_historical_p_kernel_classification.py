"""Synthetic checks of the historical-P kernel value classification (Claude review P2-1).

Only the function definitions of export_historical_p_metadata.R are evaluated; no
historical RDS is read and no SLURM guard is bypassed for data access.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import unittest

HELPER = Path(__file__).resolve().parents[1] / 'src' / 'export_historical_p_metadata.R'

R_PROGRAM = r'''
args <- commandArgs(trailingOnly=TRUE)
for (e in parse(file=args[[1]])) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name('function'))) eval(e, envir=globalenv())
}
moffat <- function(d, fwhm=3.2, beta=2.7, dx=0, dy=0) {
  alpha <- (fwhm/2)/sqrt(2^(1/beta)-1)
  x <- seq_len(d)-.5-d/2-dx; y <- seq_len(d)-.5-d/2-dy
  z <- (1+outer(x^2,y^2,'+')/alpha^2)^(-beta); z/sum(z)
}
crop <- function(z, d) { o <- (nrow(z)-d)/2; z[o+seq_len(d),o+seq_len(d),drop=FALSE] }
a <- moffat(27L)
ulp <- a; ulp[14,14] <- ulp[14,14]*(1+2^-52)
cls <- function(x, y) classify_psf_values(x, y)$class
out <- list(
  exact=cls(a, a + 0),
  scaled=cls(a, 3.7*a),
  ulp=cls(a, ulp),
  different=cls(a, moffat(27L, fwhm=3.5)),
  resized=cls(a, crop(a, 21L)),
  resized_different=cls(a, crop(moffat(27L, fwhm=3.5), 21L)),
  parity=cls(a, matrix(1, 20L, 20L)),
  non_matrix=cls(a, NULL),
  non_finite=cls(a, replace(a, 1, NA)),
  centred_offset=kernel_geometry(a)$centroid_offset_from_grid_centre_xy,
  shifted_offset=kernel_geometry(moffat(61L, dx=.5, dy=.5))$centroid_offset_from_grid_centre_xy,
  ulp_metrics=classify_psf_values(a, ulp)[c('max_abs_normalised_difference','pixels_differing_common')]
)
cat(jsonlite::toJSON(out, auto_unbox=TRUE, digits=NA))
'''


def find_rscript():
    for candidate in (os.environ.get('BD_PHASEA_RSCRIPT'),
                      '/home/uqpdas/.conda/envs/waves_r/bin/Rscript', shutil.which('Rscript')):
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


@unittest.skipUnless(find_rscript(), 'Rscript unavailable')
class KernelClassificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        result = subprocess.run([find_rscript(), '--vanilla', '-e', R_PROGRAM, str(HELPER)],
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                timeout=300)
        if result.returncode:
            raise AssertionError(result.stderr)
        cls.out = json.loads(result.stdout.strip().splitlines()[-1])

    def test_classes(self):
        expected = {'exact': 'EXACT', 'scaled': 'RENORMALISED', 'ulp': 'RENORMALISED',
                    'different': 'DIFFERENT_KERNEL', 'resized': 'RESIZED',
                    'resized_different': 'DIFFERENT_KERNEL_AND_DIMENSIONS',
                    'parity': 'NOT_COMPARABLE_DIMENSION_PARITY',
                    'non_matrix': 'NOT_COMPARABLE_NON_MATRIX',
                    'non_finite': 'NOT_COMPARABLE_NON_FINITE'}
        self.assertEqual({key: self.out[key] for key in expected}, expected)

    def test_ulp_difference_is_recorded(self):
        self.assertEqual(self.out['ulp_metrics']['pixels_differing_common'], 1)
        self.assertLess(self.out['ulp_metrics']['max_abs_normalised_difference'], 1e-15)

    def test_centroid_geometry_uses_physical_grid_centre(self):
        for value in self.out['centred_offset']:
            self.assertLess(abs(value), 1e-12)
        for value in self.out['shifted_offset']:
            self.assertAlmostEqual(value, 0.5, delta=0.01)


if __name__ == '__main__':
    unittest.main(verbosity=2)
