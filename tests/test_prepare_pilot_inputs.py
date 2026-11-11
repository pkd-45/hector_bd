#!/usr/bin/env python3
"""Synthetic gzip corruption, orientation and Rfits sample-comparison checks.

No FITS input is read. The optional evidence test reads only the two small JSON
records written by Friday job 1453635 when they are present.
"""
import copy
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('prepare_pilot_inputs',ROOT/'scripts/prepare_pilot_inputs.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class InputCopyTests(unittest.TestCase):
    def test_exact_decompression_and_original_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            original=Path(directory)/'science.fits.gz';destination=Path(directory)/'copy.fits'
            payload=bytes(range(256))*30
            original.write_bytes(gzip.compress(payload))
            before=m.source_identity(original)
            evidence=m.verified_decompress(original,destination)
            self.assertEqual(destination.read_bytes(),payload)
            self.assertEqual(evidence['sha256'],hashlib.sha256(payload).hexdigest())
            self.assertTrue(evidence['gzip_crc_verified'])
            self.assertEqual(m.source_identity(original),before)
            with self.assertRaises(FileExistsError):m.verified_decompress(original,destination)

    def test_corrupt_gzip_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            original=Path(directory)/'bad.gz';destination=Path(directory)/'copy'
            encoded=bytearray(gzip.compress(b'not astronomical science'*500))
            encoded[-8]^=1
            original.write_bytes(encoded)
            with self.assertRaises(gzip.BadGzipFile):m.verified_decompress(original,destination)

    def test_samples_use_fits_x_y_in_astropy_y_x_grid(self):
        data=np.arange(1333*1333,dtype=np.float32).reshape(1333,1333)
        evidence=m.array_evidence(data)
        a=next(p for p in evidence['samples'] if (p['fits_x_1based'],p['fits_y_1based'])==(2,5))
        b=next(p for p in evidence['samples'] if (p['fits_x_1based'],p['fits_y_1based'])==(5,2))
        self.assertEqual(a['value'],4*1333+1)
        self.assertEqual(b['value'],1333+4)
        self.assertNotEqual(a['value'],b['value'])

    def test_duplicate_target_table_entry_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'table.tsv'
            path.write_text('uberID\tstatus\n'+m.UID+'\tIDENTICAL\n'+m.WID+'\tIDENTICAL\n')
            with self.assertRaisesRegex(ValueError,'exactly one row'):
                m.selected(path,'uberID',m.UID,'\t')


def sample(x,y,value):
    return dict(fits_x_1based=x,fits_y_1based=y,value=value,finite=value is not None)


# Real Friday values (job 1453635, image HDU1): Astropy doubles and the 15-significant-digit
# values R's jsonlite wrote for the same pixels. (2,5)/(5,2) is the off-diagonal pair.
ASTROPY=[sample(2,5,0.012884007766842842),sample(5,2,-0.00876860786229372)]
RFITS_15_DIGIT=[sample(2,5,0.0128840077668428),sample(5,2,-0.00876860786229372)]
RFITS_TRANSPOSED=[sample(2,5,-0.00876860786229372),sample(5,2,0.0128840077668428)]
EVIDENCE=ROOT/'outputs'/'pilot_inputs_1453635'


class RfitsSampleComparisonTests(unittest.TestCase):
    """Regression for Claude review P1-2: exercises m.check_rfits_samples used by run()."""

    def test_run_uses_production_comparison(self):
        source=(ROOT/'scripts/prepare_pilot_inputs.py').read_text()
        self.assertIn("check_rfits_samples(name,copy['array']['samples'],observed['samples'])",source)
        self.assertNotIn("expected['value']!=actual['value']",source)

    def test_15_digit_rounding_accepted(self):
        self.assertNotEqual(ASTROPY[0]['value'],RFITS_15_DIGIT[0]['value'])
        m.check_rfits_samples('image',ASTROPY,RFITS_15_DIGIT)

    def test_transposed_read_rejected(self):
        with self.assertRaisesRegex(ValueError,'off-diagonal Rfits sample differs'):
            m.check_rfits_samples('image',ASTROPY,RFITS_TRANSPOSED)

    def test_coordinate_mismatch_rejected(self):
        moved=copy.deepcopy(RFITS_15_DIGIT);moved[1]['fits_y_1based']=3
        with self.assertRaisesRegex(ValueError,'off-diagonal Rfits sample differs'):
            m.check_rfits_samples('image',ASTROPY,moved)

    def test_missing_or_extra_samples_rejected(self):
        for actual in (RFITS_15_DIGIT[:1],RFITS_15_DIGIT+[sample(9,9,1.0)],[]):
            with self.assertRaisesRegex(ValueError,'sample count'):
                m.check_rfits_samples('image',ASTROPY,actual)
        with self.assertRaisesRegex(ValueError,'sample count'):
            m.check_rfits_samples('image',[],[])

    def test_nonfinite_samples(self):
        m.check_rfits_samples('image',[sample(1,1,None)],[sample(1,1,None)])
        for expected,actual in (([sample(1,1,None)],[sample(1,1,0.0)]),
                                ([sample(1,1,0.0)],[sample(1,1,None)])):
            with self.assertRaisesRegex(ValueError,'off-diagonal Rfits sample differs'):
                m.check_rfits_samples('image',expected,actual)
        flagged=[dict(sample(1,1,0.5),finite=False)]
        with self.assertRaisesRegex(ValueError,'off-diagonal Rfits sample differs'):
            m.check_rfits_samples('image',[sample(1,1,0.5)],flagged)

    def test_tolerance_is_far_below_physical_differences(self):
        value=ASTROPY[0]['value']
        self.assertTrue(m.sample_values_agree(value,float(f'{value:.15g}')))
        self.assertFalse(m.sample_values_agree(value,value*(1+1e-9)))
        self.assertTrue(m.sample_values_agree(0.0,0.0))
        self.assertFalse(m.sample_values_agree(0.0,1e-12))

    @unittest.skipUnless((EVIDENCE/'copy_manifest.json').is_file() and
                         (EVIDENCE/'Rfits_verification.json').is_file(),
                         'Friday job 1453635 evidence not present')
    def test_friday_1453635_evidence_passes_and_old_rule_failed(self):
        manifest=json.loads((EVIDENCE/'copy_manifest.json').read_text())
        verification=json.loads((EVIDENCE/'Rfits_verification.json').read_text())
        old_rule_failures=0
        for name in ('image','segmap'):
            expected=manifest['copies'][name]['array']['samples']
            actual=verification['products'][name]['samples']
            self.assertEqual(len(expected),10)
            m.check_rfits_samples(name,expected,actual)
            old_rule_failures+=sum(e['value']!=a['value'] for e,a in zip(expected,actual))
            # A transposed R read of the same products must still be rejected.
            swapped={(s['fits_x_1based'],s['fits_y_1based']):s['value'] for s in actual}
            transposed=[dict(s,value=swapped[(s['fits_y_1based'],s['fits_x_1based'])])
                        if (s['fits_y_1based'],s['fits_x_1based']) in swapped else s for s in actual]
            if any(not m.sample_values_agree(a['value'],t['value']) for a,t in zip(actual,transposed)):
                with self.assertRaises(ValueError):
                    m.check_rfits_samples(name,expected,transposed)
        self.assertGreater(old_rule_failures,0)


if __name__=='__main__':unittest.main(verbosity=2)
