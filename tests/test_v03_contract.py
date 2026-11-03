"""Structural guard for V03 PSF and source-audit decisions."""
import unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]
class V03Guards(unittest.TestCase):
 def code(self,part):return (P/part).read_text()
 def test_pa_sign_and_bound_clamp(self):
  self.assertIn('pa <- -atan2(',self.code('R/models.R'))
  s=self.code('R/fit_one.R')
  self.assertIn('pmin(pmax(p,spec$bounds$lower),spec$bounds$upper)',s)
  self.assertNotIn('spec$bounds$lower+1e-6',s)
 def test_psf_mode(self):
  s=self.code('R/fit_one.R')
  self.assertIn("spec$modellist$psf <- analytic",s)
  self.assertIn("family=='P' && fs!=1L",s)
  self.assertIn("psf_type=psf_mode$psf_type",s)
  self.assertIn("raw_ml$psf <- NULL",s)
 def test_f2_pinned_regression(self):
  s=self.code('R/f2_finesample_test.R')
  self.assertIn('F2_POINT_NUCLEAR_SNR 60',s)
  self.assertIn("ml$psf <- list(moffat",s)
  self.assertIn("m35$dchisq<1",s)
  self.assertIn('7L,9L',s)
  self.assertIn("identifiability>1",s)
 def test_audit_psf_fields(self):
  a=self.code('python/phasea/src/export_rds_arrays.R')
  b=self.code('python/phasea/src/export_phasea.py')
  for s in ('Data_finesample','Data_psftype'):
   self.assertIn(s,a);self.assertIn(s,b)
  self.assertIn('historical_finesample_by_fit_json',b)
  self.assertIn('historical_psftype_by_fit_json',b)
 def test_recovery_noise(self):
  s=self.code('tests/test_synthetic_recovery.R')
  self.assertIn('rnorm(prod(d)',s)
  self.assertNotIn('sin(X*.63',s)
  self.assertIn('optim_maxit=500L',s)
 def test_reproducer_quadrants_invoked(self):
  s=self.code('tests/run_checks.sh')
  self.assertIn('Rscript tests/test_moment_pa_convention.R',s)
  self.assertIn('Rscript tests/test_psf_branch_contract.R',s)
 def test_uncertified_recording(self):
  s=self.code('R/fit_one.R')
  self.assertIn('uncertified_starts=uncertified',s)
  self.assertIn("dpath <- file.path(diagnostic_dir",s)
  self.assertIn('bound_hits=which(',s)
if __name__=='__main__':unittest.main()
