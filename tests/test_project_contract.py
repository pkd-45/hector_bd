import json, pathlib, re, unittest, hashlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
class Layout(unittest.TestCase):
 def test_single_project_marker(self):
  self.assertTrue((ROOT/'.hector_bd_root').exists())
  for p in ['config/step0.yaml','config/paths_mac.yaml','config/paths_friday.yaml',
   'config/run_pilot.yaml','R/00_load_all.R','R/io_fits.R','R/segment_contract.R',
   'R/build_inputs.R','R/sky_noise.R','R/psf_v14.R','R/models.R',
   'R/fit_one.R','R/postprocess.R','R/cards.R','R/renderer_tests.R',
   'scripts/run_one_galaxy.R','slurm/submit_array.sh','env/install_pinned.R',
   '.vscode/settings.json','notebooks/01_explore_one_galaxy.ipynb']:
   self.assertTrue((ROOT/p).is_file(),p)
 def test_shared_code_path(self):
  nb=json.loads((ROOT/'notebooks/01_explore_one_galaxy.ipynb').read_text())
  code='\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type']=='code')
  engine=(ROOT/'R/fit_one.R').read_text()
  cli=(ROOT/'scripts/run_one_galaxy.R').read_text()
  self.assertIn('bd_run_one(id,cfg,run_id,mode,',cli)
  self.assertIn('for (stage in stages) bd_stage(ctx,stage)',engine)
  for stage in ['inputs','g2','f1_f9','fit','postprocess','card']:
   self.assertRegex(code,rf'bd_stage\(ctx,\s*"{stage}"\)')
  self.assertNotIn('ProFit::profitSetupData(',code)
  self.assertNotIn('stats::optim(',code)
 def test_amendment_a1_boundary_flag(self):
  src=(ROOT/'R/segment_contract.R').read_text()
  self.assertIn('boundary_flag <- !all(',src)
  self.assertIn("boundary_flag=boundary_flag",src)
  self.assertNotIn('catalogue centre falls on a 3x3 original-label boundary',src)
 def test_d1_source_byte_identical_to_approved(self):
  # The archive's sources copied unchanged to python/phasea/src.
  ref=ROOT/'python/phasea/src/export_phasea.py'
  self.assertIn('def compare_family(',ref.read_text())
  self.assertIn('original_background',ref.read_text())
 def test_no_approval_by_default(self):
  s=(ROOT/'config/run_pilot.yaml').read_text()
  for key in ['pilot_approved','full_832_approved','psf_g5_approved','models_g6_approved']:
   self.assertRegex(s,rf'{key}: false')
 def test_pilot_has_only_10_avoids_invented_ids(self):
  s=(ROOT/'config/run_pilot.yaml').read_text()
  self.assertIn('frozen_roster_required: R61_832_ONLY',s)
  self.assertEqual(len(re.findall(r'^  - W[0-9]+$',s,re.M)),20) # 10 reference + 10 working
 def test_json_notebooks_are_valid(self):
  for f in (ROOT/'notebooks').glob('*.ipynb'):
   data=json.loads(f.read_text());self.assertEqual(data['nbformat'],4)
   self.assertEqual(data['metadata']['kernelspec']['name'],'waves_r_bd')
if __name__=='__main__':unittest.main()
