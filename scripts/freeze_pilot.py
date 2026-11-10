#!/usr/bin/env python3
"""Prepare a review-only pilot from exact R61 tables and existing segment audit.

Uses recorded original brightness/size strata, native model and resolution bin.
No FITS reads, observed fitting, configuration edits or approval changes.
"""
import argparse, csv, hashlib, json, math, os, re
from collections import Counter
from pathlib import Path

ROSTER_SHA256='874a350440b0361371db15b4f4e363bfc79eb924a31ff0ee5ae478ca41d275f3'
POLICY='R53_NATIVE_S_P_D_FIXED4_MODEL_SELECTION_V1'
MODEL_COUNTS={'S':352,'P':68,'D_FIXED4':412}
SEED=20261031
MAX_WCS_ARCSEC=.01
RESOLUTION_BINS={'Q1_LOWEST':'Q1','Q2':'Q2','Q3':'Q3','Q4_HIGHEST':'Q4'}
ORIGINAL_ANCHORS=('W42710765002459','W221990044410325','W351690654304280',
 'W176980350403709','W12690079611187','W41691060304604','W350680203312123',
 'W18680557409852','W49710254502112','W213991001205426')
REPLACED='W12690079611187'
REPLACEMENT='W11690092907692'
ANCHORS=tuple(REPLACEMENT if x==REPLACED else x for x in ORIGINAL_ANCHORS)
DEFAULT_RESULTS=Path('/scratch/uqpdas/WAVES_BD_PROFIT_PROFUSE_R61_R53_NATIVE_SPD_FIXED4_POLICY_EXECUTION_20260930/results')

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def normalize_id(value):
    match=re.fullmatch(r'W?([1-9][0-9]*)',str(value).strip())
    if not match:raise ValueError(f'Invalid immutable WAVES ID: {value!r}')
    return 'W'+match[1]

def truth(value):
    text=str(value).strip().lower()
    if text in ('true','1'):return True
    if text in ('false','0'):return False
    raise ValueError(f'Expected explicit Boolean: {value!r}')

def read_table(path,required):
    with Path(path).open(newline='',encoding='utf-8-sig') as stream:
        header=stream.readline();stream.seek(0)
        delimiter='\t' if header.count('\t')>header.count(',') else ','
        reader=csv.DictReader(stream,delimiter=delimiter)
        missing=set(required)-set(reader.fieldnames or [])
        if missing:raise ValueError(f'{path}: missing columns {sorted(missing)}')
        rows={}
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'{path}: malformed table row')
            uid=normalize_id(row['uberID'])
            if uid in rows:raise ValueError(f'{path}: duplicate immutable ID {uid}')
            rows[uid]=row
    return rows

def join_inputs(roster_path,science_path,audit_path,*,expected_roster_sha=ROSTER_SHA256):
    if sha256(roster_path)!=expected_roster_sha:raise ValueError('Frozen R61 roster SHA256 mismatch')
    roster=read_table(roster_path,('uberID','native_policy','global_candidate'))
    science=read_table(science_path,('uberID','stratum','input_resolution_quartile',
                                   'r61_global_candidate','r61_native_policy'))
    audit=read_table(audit_path,('uberID','status','wcs_max_sep_arcsec','cat_boundary_flag',
                    'cat_agrees_with_historical','cat_label_is_sky','catalogue_row_found'))
    if any(len(t)!=832 for t in (roster,science,audit)):
        raise ValueError('Expected exactly832 unique rows in roster/science/audit')
    if roster.keys()!=science.keys() or roster.keys()!=audit.keys():
        raise ValueError('832-object ID sets disagree across roster/science/audit')
    counts=Counter(r['global_candidate'] for r in roster.values())
    if dict(counts)!=MODEL_COUNTS:raise ValueError(f'Frozen native class counts differ: {dict(counts)}')
    joined={};exclusions=[]
    for uid in sorted(roster):
        r,s,a=roster[uid],science[uid],audit[uid]
        if r['native_policy']!=POLICY or s['r61_native_policy']!=POLICY:
            raise ValueError(f'{uid}: native R61 policy mismatch')
        if r['global_candidate']!=s['r61_global_candidate']:
            raise ValueError(f'{uid}: native model disagrees with science catalogue')
        match=re.fullmatch(r'M([1-4])(?:_[A-Z0-9]+)*__R([1-4])(?:_[A-Z0-9]+)*',s['stratum'])
        if not match or s['input_resolution_quartile'] not in RESOLUTION_BINS:
            raise ValueError(f'{uid}: unsupported recorded brightness/size or resolution stratum')
        reasons=[]
        if a['status']!='IDENTICAL':reasons.append('ORIGINAL_SEGMENTATION_NOT_IDENTICAL')
        try:wcs=float(a['wcs_max_sep_arcsec'])
        except (TypeError,ValueError):wcs=float('nan')
        if not math.isfinite(wcs) or wcs>MAX_WCS_ARCSEC or wcs<0:
            reasons.append('WCS_EXCEEDS_0P01_ARCSEC_OR_UNVERIFIED')
        for field,wanted,reason in (
            ('cat_boundary_flag',False,'CATALOGUE_BOUNDARY_FLAG'),
            ('cat_agrees_with_historical',True,'CATALOGUE_LABEL_DISAGREES_OR_UNVERIFIED'),
            ('cat_label_is_sky',False,'CATALOGUE_LABEL_IS_SKY_OR_UNVERIFIED'),
            ('catalogue_row_found',True,'CATALOGUE_ROW_MISSING_OR_UNVERIFIED')):
            try:observed=truth(a[field])
            except ValueError:observed=None
            if observed is not wanted:reasons.append(reason)
        item=dict(uberID=uid,numeric_uberID=uid[1:],native_model=r['global_candidate'],
            stratum=s['stratum'],original_brightness_bin='M'+match[1],original_size_bin='R'+match[2],
            input_resolution_quartile=s['input_resolution_quartile'],
            resolution_quartile=RESOLUTION_BINS[s['input_resolution_quartile']],
            wcs_max_sep_arcsec=wcs if math.isfinite(wcs) else None,
            audit_status=a['status'],eligible=not reasons,exclusion_reasons=reasons)
        joined[uid]=item
        if reasons:exclusions.append(item)
    sources={name:dict(path=str(Path(path).resolve()),sha256=sha256(path),rows=832)
        for name,path in (('r61_roster',roster_path),('r61_science_catalogue',science_path),
                          ('existing_original_segmap_audit',audit_path))}
    return joined,exclusions,sources

def coverage(rows):
    return {field:dict(sorted(Counter(r[field] for r in rows).items())) for field in
        ('native_model','original_brightness_bin','original_size_bin','stratum',
         'input_resolution_quartile','resolution_quartile')}

def choose_additional(joined,anchors=ANCHORS,count=10,seed=SEED):
    """Balance available native classes, original strata and bins after anchors.

    Deterministic SHA256(seed:ID) breaks ties independently of input row order.
    This purposeful pilot is not a probability sample of the survey population.
    """
    anchors=tuple(normalize_id(x) for x in anchors)
    if len(set(anchors))!=len(anchors):raise ValueError('Duplicate anchor ID')
    for uid in anchors:
        if uid not in joined or not joined[uid]['eligible']:raise ValueError(f'Anchor missing/ineligible: {uid}')
    selected=[joined[uid] for uid in anchors]
    pool={uid:r for uid,r in joined.items() if r['eligible'] and uid not in anchors}
    if len(pool)<count:raise ValueError('Insufficient eligible R61 objects')
    picks=[];records=[]
    for step in range(count):
        c=coverage(selected)
        def score(item):
            uid,r=item
            marginal=[c[f].get(r[f],0) for f in
                      ('original_brightness_bin','original_size_bin','resolution_quartile')]
            return (c['native_model'].get(r['native_model'],0),c['stratum'].get(r['stratum'],0),
                sum(marginal),max(marginal),c['resolution_quartile'].get(r['resolution_quartile'],0),
                hashlib.sha256(f'{seed}:{uid}'.encode('ascii')).hexdigest())
        uid,row=min(pool.items(),key=score)
        records.append(dict(selection_step=step+1,**row,selection_score=list(score((uid,row))),
            selection_reason='Least-represented available native model, then recorded original '
            'brightness/size stratum, marginal bins and resolution coverage; seeded SHA256 tie-break.'))
        picks.append(uid);selected.append(row);del pool[uid]
    return picks,records

def prepare(roster_path,science_path,audit_path,*,expected_roster_sha=ROSTER_SHA256):
    joined,exclusions,sources=join_inputs(roster_path,science_path,audit_path,
                                        expected_roster_sha=expected_roster_sha)
    if REPLACED not in joined or joined[REPLACED]['eligible']:
        raise ValueError('Specified anchor replacement requires original-anchor exclusion evidence')
    picks,records=choose_additional(joined)
    pilot=list(ANCHORS)+picks
    if len(pilot)!=20 or len(set(pilot))!=20:raise ValueError('Expected20 unique pilot objects')
    return dict(source_layer='R61_832_ONLY',sha256=sources['r61_roster']['sha256'],
        native_policy=POLICY,native_model_counts=MODEL_COUNTS,n_rows=832,sources=sources,random_seed=SEED,
        strata_columns=['global_candidate','stratum','input_resolution_quartile'],
        missing_original_covariates=['R50','physical_Sersic_n','crowding'],
        covariate_reconciliation='Use recorded original brightness/size bins, native R61 class '
        'and input resolution quartile. No nonexistent covariates are fabricated. '
        'Log-scaled production nser parameters are not interpreted as physical Sersic n.',
        eligibility_rule='Existing audit only: IDENTICAL labels; finite WCS<=0.01arcsec; '
        'independent catalogue row found; exact label agreement; not sky; no3x3 boundary flag.',
        original_known_10=list(ORIGINAL_ANCHORS),known_10=list(ANCHORS),
        anchor_replacement=dict(removed=REPLACED,added=REPLACEMENT,removed_audit=joined[REPLACED],
            replacement_audit=joined[REPLACEMENT],reason='Requested replacement of WCS2.894172arcsec '
            'anchor by catalogue-consistent W11690092907692 with WCS0arcsec.'),
        n_eligible=sum(r['eligible'] for r in joined.values()),exclusions=exclusions,
        stratified_10=picks,pilot_20=pilot,
        pilot_20_sha256=hashlib.sha256(''.join(x+'\n' for x in pilot).encode('ascii')).hexdigest(),
        selected_records=[dict(selection_kind='ANCHOR',**joined[uid]) for uid in ANCHORS]+records,
        coverage_eligible=coverage([r for r in joined.values() if r['eligible']]),
        coverage_anchors=coverage([joined[uid] for uid in ANCHORS]),
        coverage_additional=coverage([joined[uid] for uid in picks]),
        coverage_pilot=coverage([joined[uid] for uid in pilot]),
        design_limit='Purposeful balanced pilot conditional on anchors; not a random representative '
                     'survey sample or scientific approval.',
        status='CANDIDATE_20_FOR_REVIEW_NO_OBSERVED_FIT_AUTHORIZATION',
        config_modified=False,approval_flags_modified=False,script_sha256=sha256(__file__))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--r61-roster',type=Path,default=DEFAULT_RESULTS/'R61_NATIVE_MODEL_SELECTION_832.tsv')
    ap.add_argument('--r61-science',type=Path,default=DEFAULT_RESULTS/'R61_SCIENCE_CATALOGUE_WITH_NATIVE_SELECTION_832.tsv')
    ap.add_argument('--audit',type=Path,default=Path('outputs/step1_catalogue_1427444/step1_segmap_audit.tsv'))
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):ap.error('Run candidate selection in a SLURM job')
    root=Path(__file__).resolve().parents[1];output=args.output.expanduser().resolve()
    if not any(parent in output.parents for parent in (root/'outputs',root/'config')):
        ap.error('Output must be a NEW JSON below repository outputs/ or config/')
    if output.exists() or output.suffix.lower()!='.json':ap.error('Refusing overwrite or non-JSON output')
    result=prepare(args.r61_roster,args.r61_science,args.audit)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print('CANDIDATE_20_FOR_REVIEW_ONLY',output,'ROSTER_SHA256',result['sha256'],
          'PILOT_IDS_SHA256',result['pilot_20_sha256'])
if __name__=='__main__':main()
