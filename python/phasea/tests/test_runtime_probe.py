import json
import os
import stat
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from export_phasea import create_runtime


class TestRuntimeProbe(unittest.TestCase):
    def run_stub(self, body):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            script = root / "Rscript"
            script.write_text("#!/bin/sh\n" + body + "\n")
            script.chmod(script.stat().st_mode | stat.S_IXUSR)
            output = root / "out"
            output.mkdir()
            with patch.dict(os.environ, {}, clear=False):
                return create_runtime(str(script), output, str(root / "rlib"))

    def test_nonzero_exit_fails_closed_after_writing_report(self):
        with self.assertRaisesRegex(RuntimeError, "NONZERO_EXIT"):
            self.run_stub("printf 'partial output\\n'; exit 1")

    def test_missing_sentinel_fails_closed(self):
        with self.assertRaisesRegex(RuntimeError, "SENTINEL_MISSING"):
            self.run_stub("printf 'probe output\\n'; exit 0")

    def test_success_requires_sentinel(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            script = root / "Rscript"
            script.write_text("#!/bin/sh\nprintf 'RUNTIME_PROBE=COMPLETE\\n'\n")
            script.chmod(script.stat().st_mode | stat.S_IXUSR)
            output = root / "out"
            output.mkdir()
            runtime_id = create_runtime(str(script), output, str(root / "rlib"))
            report = json.loads((output / "runtime.json").read_text())
            self.assertTrue(runtime_id)
            self.assertEqual(report["runtime_gate"], "PASS")
            self.assertEqual(report["runtime_gate_failures"], [])


if __name__ == "__main__":
    unittest.main()
