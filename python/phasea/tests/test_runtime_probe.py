import json, os, shutil, stat, sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from export_phasea import create_runtime, RUNTIME_PROBE_SENTINEL

PACKAGES=('ProFit','ProFuse','ProFound','Highlander','LaplacesDemon','Rfits')


def find_rscript():
    for c in (os.environ.get('BD_PHASEA_RSCRIPT'),'/home/uqpdas/.conda/envs/waves_r/bin/Rscript',shutil.which('Rscript')):
        if c and os.path.isfile(c) and os.access(c,os.X_OK):return c
    return None


class TestRuntimeProbeGate(unittest.TestCase):
    def setUp(self):
        self.tmp=Path(tempfile.mkdtemp(prefix='bd_runtime_probe_'))
        self.addCleanup(shutil.rmtree,self.tmp,ignore_errors=True)

    def fake_rscript(self,body):
        p=self.tmp/'fake_Rscript';p.write_text('#!/bin/sh\n'+body+'\n')
        p.chmod(p.stat().st_mode|stat.S_IXUSR);return str(p)

    def report(self):
        with open(self.tmp/'runtime.json') as f:return json.load(f)

    def test_nonzero_exit_fails_closed_and_keeps_evidence(self):
        # Mirrors the Phase A V04 failure: first lines print, then a parse error.
        rs=self.fake_rscript("echo 'R_VERSION=R version 4.3.3'\necho 'Error: unexpected symbol' >&2\nexit 1")
        with self.assertRaisesRegex(RuntimeError,'RUNTIME_PROBE_NONZERO_EXIT:1'):
            create_runtime(rs,self.tmp,str(self.tmp))
        r=self.report()
        self.assertEqual(r['runtime_gate'],'FAIL');self.assertEqual(r['returncode'],1)
        self.assertIn('unexpected symbol',r['stderr'])

    def test_zero_exit_without_sentinel_fails_closed(self):
        rs=self.fake_rscript("echo 'R_VERSION=R version 4.3.3'\nexit 0")
        with self.assertRaisesRegex(RuntimeError,'RUNTIME_PROBE_INCOMPLETE_NO_SENTINEL'):
            create_runtime(rs,self.tmp,str(self.tmp))
        self.assertEqual(self.report()['runtime_gate'],'FAIL')

    def test_real_probe_runs_to_completion(self):
        rs=find_rscript()
        if rs is None:self.skipTest('Rscript unavailable; set BD_PHASEA_RSCRIPT')
        pinned=os.environ.get('BD_PHASEA_HISTORICAL_RLIB',str(self.tmp/'no_pinned_lib'))
        rid=create_runtime(rs,self.tmp,pinned)
        r=self.report()
        self.assertEqual(r['returncode'],0,r['stderr'])
        self.assertEqual(r['runtime_gate'],'PASS');self.assertEqual(r['runtime_id'],rid)
        lines=r['stdout'].splitlines()
        self.assertEqual(lines[-1],RUNTIME_PROBE_SENTINEL)
        for p in PACKAGES:
            for scope in ('CURRENT','PINNED'):
                self.assertTrue(any(l.startswith(f'PACKAGE={p} SCOPE={scope} ') for l in lines),(p,scope))


if __name__=='__main__':
    unittest.main()
