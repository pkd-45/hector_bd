#!/usr/bin/env python3
"""Synthetic gzip corruption and orientation checks; no observed input reads."""
import csv
import gzip
import hashlib
import importlib.util
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


if __name__=='__main__':unittest.main(verbosity=2)
