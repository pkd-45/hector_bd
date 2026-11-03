"""Static F2 contract checks. They do NOT substitute for R/ProFit execution."""
import csv, hashlib, math, re, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class F2Contract(unittest.TestCase):
 def test_frozen_psf_sha_and_flux(self):
  path=ROOT/'tests/fixtures/f2_moffat_21px.csv'
  actual=hashlib.sha256(path.read_bytes()).hexdigest()
  self.assertEqual(actual,'d569969a59f40745c323250f888e3e122a457f9eb63e78029fd684982b2eb1d0')
  with path.open() as f:
   x=[[float(i) for i in row] for row in csv.reader(f)]
  self.assertEqual(len(x),21)
  self.assertTrue(all(len(row)==21 and all(math.isfinite(y) and y>0 for y in row) for row in x))
  self.assertAlmostEqual(sum(sum(row) for row in x),1,places=12)
 def test_real_profits_setup_not_raw_oversampling(self):
  text=(ROOT/'src/f2_finesample_test.R').read_text()
  self.assertIn('psffinesampled=FALSE',text)
  self.assertIn('ProFit::profitSetupData(',text)
  self.assertIn('ProFit::profitLikeModel(',text)
  self.assertIn("c(1L,3L,5L)",text)
  self.assertIn('m35$dchisq<1',text)
  self.assertIn('m35$max<.1',text)
  self.assertIn('F2_PARAMETER_PROXY',text)
  self.assertIn('F2_METRIC',text)
  self.assertNotIn('optim(',text)
 def test_runner_executes_f2(self):
  r=(ROOT/'src/renderer_tests.R').read_text()
  shell=(ROOT/'run_synthetic_tests.sh').read_text()
  self.assertIn("source('src/f2_finesample_test.R')",r)
  self.assertIn("F2=.bd_test('F2',function(){ bd_f2_run() })",r)
  self.assertIn('Rscript src/renderer_tests.R',shell)
if __name__=='__main__':unittest.main()
