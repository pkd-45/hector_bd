#!/usr/bin/env python3
"""Synthetic R61 table/eligibility tests; no science-image reads or fitting."""
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('freeze_pilot',ROOT/'scripts/freeze_pilot.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def write_table(path,rows,delimiter='\t'):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter=delimiter)
        writer.writeheader();writer.writerows(rows)


class PilotTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name)
        self.roster=self.folder/'roster.tsv'
        self.science=self.folder/'science.tsv'
        self.audit=self.folder/'audit.tsv'
        ids=list(m.ORIGINAL_ANCHORS)+[m.REPLACEMENT]
        ids += ['W'+str(90000000000000+i) for i in range(832-len(ids))]
        models=['S']*352+['P']*68+['D_FIXED4']*412
        self.r=[];self.s=[];self.a=[]
        for i,(uid,model) in enumerate(zip(ids,models)):
            self.r.append(dict(uberID=uid[1:],native_policy=m.POLICY,global_candidate=model))
            self.s.append(dict(uberID=uid,stratum=f'M{i%4+1}_BRIGHTNESS__R{(i//4)%4+1}_SIZE',
                input_resolution_quartile=('Q1_LOWEST','Q2','Q3','Q4_HIGHEST')[(i//16)%4],r61_global_candidate=model,
                r61_native_policy=m.POLICY,production_param__S__sersic_nser='-0.04'))
            self.a.append(dict(uberID=uid[1:],status='IDENTICAL',wcs_max_sep_arcsec='0.0',
                cat_boundary_flag='False',cat_agrees_with_historical='True',
                cat_label_is_sky='False',catalogue_row_found='True'))
        self.a[4]['wcs_max_sep_arcsec']='2.894172'
        self.a[20]['wcs_max_sep_arcsec']='0.010001'
        self.a[21]['cat_boundary_flag']='True'
        self.a[22]['cat_agrees_with_historical']='False'
        self.a[23]['wcs_max_sep_arcsec']='0.01'
        self.save()

    def save(self,delimiter='\t'):
        for path,rows in ((self.roster,self.r),(self.science,self.s),(self.audit,self.a)):
            write_table(path,rows,delimiter)

    def prepare(self):
        return m.prepare(self.roster,self.science,self.audit,expected_roster_sha=m.sha256(self.roster))

    def test_real_numeric_tsv_schema_and_recorded_covariates_join(self):
        before={p:m.sha256(p) for p in (self.roster,self.science,self.audit)}
        result=self.prepare()
        self.assertEqual(len(result['pilot_20']),20)
        self.assertEqual(len(set(result['pilot_20'])),20)
        self.assertEqual(result['known_10'],list(m.ANCHORS))
        self.assertIn(m.REPLACEMENT,result['known_10'])
        self.assertNotIn(m.REPLACED,result['pilot_20'])
        self.assertEqual(result['native_model_counts'],{'S':352,'P':68,'D_FIXED4':412})
        self.assertEqual(result['coverage_pilot']['native_model'],{'S':10,'P':5,'D_FIXED4':5})
        self.assertEqual(result['missing_original_covariates'],['R50','physical_Sersic_n','crowding'])
        self.assertFalse(result['config_modified'])
        self.assertFalse(result['approval_flags_modified'])
        self.assertEqual(before,{p:m.sha256(p) for p in before})
        self.assertTrue(all(r['stratum'].startswith('M') for r in result['selected_records']))
        self.assertEqual(set(result['coverage_eligible']['input_resolution_quartile']),
                         {'Q1_LOWEST','Q2','Q3','Q4_HIGHEST'})
        self.assertEqual(set(result['coverage_eligible']['resolution_quartile']),{'Q1','Q2','Q3','Q4'})

    def test_selection_stable_under_input_order_and_csv_tsv_encoding(self):
        a=self.prepare()
        self.r.reverse();self.s.reverse();self.a.reverse()
        self.save(delimiter=',')
        b=self.prepare()
        self.assertEqual(a['pilot_20'],b['pilot_20'])
        self.assertEqual(a['pilot_20_sha256'],b['pilot_20_sha256'])
        self.assertEqual(a['coverage_pilot'],b['coverage_pilot'])

    def test_strict_wcs_boundary_and_catalogue_eligibility(self):
        joined,excluded,_=m.join_inputs(self.roster,self.science,self.audit,
                                       expected_roster_sha=m.sha256(self.roster))
        self.assertFalse(joined[m.normalize_id(self.a[20]['uberID'])]['eligible'])
        self.assertFalse(joined[m.normalize_id(self.a[21]['uberID'])]['eligible'])
        self.assertFalse(joined[m.normalize_id(self.a[22]['uberID'])]['eligible'])
        self.assertTrue(joined[m.normalize_id(self.a[23]['uberID'])]['eligible'])
        self.assertEqual(len(excluded),4)
        result=self.prepare()
        excluded_ids={r['uberID'] for r in excluded}
        self.assertFalse(excluded_ids.intersection(result['pilot_20']))

    def test_sha_identity_and_native_policy_cannot_drift(self):
        with self.assertRaisesRegex(ValueError,'SHA256'):
            m.prepare(self.roster,self.science,self.audit)
        self.s[100]['r61_native_policy']='UNRELATED_POLICY'
        self.save()
        with self.assertRaisesRegex(ValueError,'policy'):
            self.prepare()

    def test_model_and_832_id_sets_must_agree(self):
        self.s[500]['r61_global_candidate']='P'
        self.save()
        with self.assertRaisesRegex(ValueError,'model disagrees'):
            self.prepare()
        self.s[500]['r61_global_candidate']='D_FIXED4'
        self.s[500]['uberID']='W999999999999999'
        self.save()
        with self.assertRaisesRegex(ValueError,'ID sets disagree'):
            self.prepare()

    def test_normalization_does_not_accept_rounded_or_duplicate_ids(self):
        self.assertEqual(m.normalize_id('11690092907692'),'W11690092907692')
        self.assertEqual(m.normalize_id('W11690092907692'),'W11690092907692')
        for bad in ('11690092907692.0','1.1690092907692e13','C11690092907692','W00123'):
            with self.assertRaises(ValueError):m.normalize_id(bad)
        self.r[12]['uberID']='W'+self.r[11]['uberID']
        self.save()
        with self.assertRaisesRegex(ValueError,'duplicate immutable'):
            self.prepare()

    def test_other_ineligible_anchor_cannot_be_silently_replaced(self):
        self.a[0]['cat_boundary_flag']='True'
        self.save()
        with self.assertRaisesRegex(ValueError,'Anchor missing/ineligible'):
            self.prepare()


if __name__=='__main__':
    unittest.main(verbosity=2)
