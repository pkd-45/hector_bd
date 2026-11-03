"""Static guard tests; not substitutes for the R/ProFit runtime suite."""
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1]
class ReviewGuards(unittest.TestCase):
    def read(self,name): return (ROOT/name).read_text()
    def test_original_label_matrix_dimensions(self):
        self.assertIn('which(matrix(L %in% A_g, nrow(L), ncol(L)), arr.ind=TRUE)',
                      self.read('R/segment_contract.R'))
    def test_nested_free_start_and_convergence(self):
        s=self.read('R/fit_one.R')
        self.assertIn('nested_from_D_FIXED4',s)
        self.assertIn("BD_FREE_N worse than nested D_FIXED4",s)
        self.assertIn('convergence_message=',s)
        self.assertIn('fit$stationary',s)
        self.assertNotIn('1e100',s)
    def test_profits_correct_api(self):
        combined='\n'.join(self.read('R/'+x+'.R') for x in
            ('fit_one','f2_finesample_test','renderer_tests'))
        self.assertNotIn('profitRemakeModelList(',combined)
        self.assertNotIn("algo.func=''",combined)
        self.assertIn('profitRemakeModellist',combined)
    def test_input_segmentation_not_dilated(self):
        self.assertNotIn('profoundMakeSegimExpand',self.read('R/build_inputs.R'))
        self.assertIn('bd_sky_exclusion',self.read('R/sky_noise.R'))
        self.assertIn('original_L_sha256=',self.read('R/sky_noise.R'))
    def test_f2_analytic_moffat_accuracy(self):
        s=self.read('R/f2_finesample_test.R')
        self.assertIn('.bd_f2_analytic_moffat',s)
        self.assertIn('9L,psf,d,region,sig',s)
        self.assertIn("'analytical'",s)
        self.assertIn("m35$dchisq<1",s)
    def test_all_requested_pins(self):
        import yaml
        pins=yaml.safe_load(self.read('config/package_pins.yaml'))['packages']
        self.assertEqual(pins['ProFit']['repo'],'ICRAR/ProFit')
        for x in ['ProFit','ProFound','Highlander','LaplacesDemon','Rfits',
                  'Rwcs','digest','magicaxis','celestial','ParmOff','cmaeshpc']:
            self.assertEqual(len(pins[x]['sha']),40,x)
    def test_recovery_is_invoked_in_r_suite(self):
        self.assertIn('Rscript tests/test_synthetic_recovery.R',
                      self.read('tests/run_checks.sh'))
        self.assertIn('bd_fit_four_models(',
                      self.read('tests/test_synthetic_recovery.R'))
if __name__=='__main__': unittest.main()
