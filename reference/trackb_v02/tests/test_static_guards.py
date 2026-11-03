"""Source-contract checks only: NOT R execution or scientific renderer validation."""
import hashlib,re,unittest
from pathlib import Path
SRC=Path(__file__).resolve().parents[1]/'src'
class TestStatic(unittest.TestCase):
 def test_v13_patch_stable_and_fixes_all_reported_mad_sites(self):
  s=(SRC/'profuse_psf_v14.R').read_text()
  for x in ('stats::mad(v[keep]) * 1.4826',
      'stats::mad(sky_pix2, na.rm = TRUE) * 1.4826',
      'stats::mad(as.numeric(data_stamp), na.rm = TRUE) * 1.4826'):
   self.assertNotIn(x,s)
  self.assertEqual(s.count('.v14_axial_pa_clip(angs)'),2)
  self.assertEqual(s.count('BLOCKED_MEASURE_PSF_EE_BEFORE_PILOT'),2)
  self.assertIn('min_dim = 21L',s)
 def test_four_modules_and_no_implicit_region(self):
  d2=(SRC/'segment_contract.R').read_text(); d3=(SRC/'build_inputs.R').read_text()
  self.assertIn('T_g & Q',d2);self.assertIn('original_region, contract$F_g',d2)
  self.assertIn('region=region',d3);self.assertIn('g2_full[y,x] <- input$region',d3)
  self.assertNotIn('profuseFound2Fit(',d3)
  self.assertNotIn('profitSetupData(',d3)
 def test_renderer_gates_declared_and_d7_nonexecution(self):
  r=(SRC/'renderer_tests.R').read_text()
  for n in range(1,10): self.assertRegex(r,r'\bF'+str(n)+r'\s*=')
  self.assertIn("F2=.bd_test('F2',function(){ bd_f2_run() })",r)
  f2=(SRC/'f2_finesample_test.R').read_text()
  self.assertIn('ProFit::profitSetupData(',f2)
  self.assertIn('psffinesampled=FALSE',f2)
  self.assertIn("ProFit::profitLikeModel(",f2)
  self.assertNotIn('optim(',f2)
  self.assertIn('m35$dchisq<1',f2)
  d7=(SRC/'fit_runner_D7_DRAFT.R').read_text()
  self.assertNotIn('profitSetupData(',d7)
  self.assertNotIn('Highlander::',d7)
 def test_all_sources_have_balanced_brackets_lexically(self):
  # Approximate lexical balance only (does NOT replace parse(file=...)).
  for path in SRC.glob('*.R'):
   s=path.read_text();stack=[]; quote=None; escape=False; comment=False
   pairs={')':'(',']':'[','}':'{'}
   for ch in s:
    if comment:
     if ch=='\n':comment=False
     continue
    if quote:
     if escape:escape=False
     elif ch=='\\':escape=True
     elif ch==quote:quote=None
     continue
    if ch=='#':comment=True;continue
    if ch in ('"',"'",'`'):quote=ch;continue
    if ch in '([{':stack.append(ch)
    elif ch in ')]}':
     self.assertTrue(stack,f'{path}: unexpected {ch}')
     self.assertEqual(stack.pop(),pairs[ch],f'{path}: unmatched {ch}')
   self.assertFalse(stack,f'{path}: unclosed {stack[-10:]}')
   self.assertIsNone(quote,f'{path}: unclosed string')
if __name__=='__main__':unittest.main()
