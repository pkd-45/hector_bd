"""Offline checks of the historical-P driver summary and source provenance (no RDS reads)."""
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest

DRIVER = Path(__file__).resolve().parents[1] / 'step3_historical_p_review.py'
spec = importlib.util.spec_from_file_location('step3_historical_p_review', DRIVER)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def record(uid, hash_comparison, value_class, cache, geometry):
    field = {'presence': 'PRESENT', 'value': 1}
    return {'uberID': uid, 'metadata_flags': ['source:Data_psffinesampled:ABSENT'],
            'accepted': {'Data_finesample': field, 'Data_psftype': {'presence': 'PRESENT', 'value': 'empirical'}},
            'accepted_source_psf_comparison': hash_comparison,
            'accepted_source_psf_value_comparison': {'class': value_class, 'same_dimensions': True,
                'sum_a': 1.0, 'sum_b': 1.0, 'max_abs_raw_difference_common': 1e-18,
                'max_abs_normalised_difference': 1e-18, 'pixels_differing_common': 3, 'pixels_common': 441},
            'accepted_cache_kernel_value_comparisons': cache, 'accepted_psf_geometry': geometry}


class SummaryTests(unittest.TestCase):
    def test_summary_counts_and_differences(self):
        records = [
            record('1', 'EXACT_NUMERIC_KERNEL_MATCH', 'EXACT',
                   {'psf': {'class': 'RENORMALISED'}, 'psf_default': {'class': 'DIFFERENT_KERNEL_AND_DIMENSIONS'}},
                   {'centroid_offset_from_grid_centre_xy': [0.49, 0.48]}),
            # R writes an empty named list as [] and a NULL geometry as null.
            record('2', 'DIFFERENT_NUMERIC_KERNEL_VALUES', 'RENORMALISED', [], None),
        ]
        summary = m.summarise({'job_id': 'test'}, records)
        self.assertEqual(summary['job_id'], 'test')
        self.assertEqual(summary['accepted_source_kernel_value_classification'],
                         {'EXACT': 1, 'RENORMALISED': 1})
        self.assertEqual(summary['accepted_cache_kernel_value_classification'],
                         {'psf:RENORMALISED': 1, 'psf_default:DIFFERENT_KERNEL_AND_DIMENSIONS': 1})
        self.assertEqual([d['uberID'] for d in summary['accepted_source_kernel_differences']], ['2'])
        self.assertEqual(summary['accepted_source_kernel_differences'][0]['class'], 'RENORMALISED')
        self.assertEqual(summary['accepted_kernel_centroid_offset_from_grid_centre_px'],
                         {'count': 1, 'x_min': 0.49, 'x_max': 0.49, 'y_min': 0.48, 'y_max': 0.48})
        self.assertFalse(summary['fit_authorized'])


class ProvenanceTests(unittest.TestCase):
    def test_uncommitted_helper_is_detected(self):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory).resolve()
            def git(*args):
                subprocess.run(['git', '-C', str(repo), *args], check=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            git('init', '-q')
            helper = repo / 'helper.R'
            helper.write_text('cat(1)\n')
            git('add', 'helper.R')
            git('-c', 'user.name=t', '-c', 'user.email=t@example.invalid', 'commit', '-q', '-m', 't')
            clean = m.git_provenance(repo, [helper])
            self.assertTrue(clean['executed_sources']['helper.R']['matches_HEAD'])
            self.assertEqual(clean['git_tracked_changes'], [])
            helper.write_text('cat(2)\n')
            dirty = m.git_provenance(repo, [helper])
            self.assertFalse(dirty['executed_sources']['helper.R']['matches_HEAD'])
            self.assertEqual(len(dirty['git_tracked_changes']), 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
