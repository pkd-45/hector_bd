#!/usr/bin/env python3
"""Freeze ten stratified IDs from the 832-object R61 CSV. No fitting or data edits."""
import argparse, csv, hashlib, random, json
from collections import defaultdict
from pathlib import Path

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--r61-roster', required=True, type=Path)
    ap.add_argument('--id-column', default='uberID')
    ap.add_argument('--size-column', default='R50')
    ap.add_argument('--n-column',default='n')
    ap.add_argument('--crowding-column',default='crowding')
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    path=args.r61_roster
    raw=path.read_bytes()
    with path.open(newline='') as stream:
        reader=csv.DictReader(stream)
        cols=(args.id_column,args.size_column,args.n_column,args.crowding_column)
        if not all(x in (reader.fieldnames or []) for x in cols):
            ap.error('R61 roster lacks named strata columns: '+str(cols))
        rows=list(reader)
    if len(rows)!=832: ap.error(f'Expected frozen R61=832, got {len(rows)}')
    ids=[r[args.id_column] for r in rows]
    if len(set(ids))!=832 or not all(i.startswith('W') and i[1:].isdigit() for i in ids):
        ap.error('Nonunique/non-W IDs; do not infer another selection layer')
    import statistics
    vals={k:sorted(float(r[k]) for r in rows) for k in cols[1:]}
    known=['W42710765002459','W221990044410325','W351690654304280',
           'W176980350403709','W12690079611187','W41691060304604',
           'W350680203312123','W18680557409852','W49710254502112','W213991001205426']
    if not set(known).issubset(ids): ap.error('Known 10 are not all in the R61 roster; cannot freeze')
    def binidx(v,a):
        return sum(float(v)>a[int((len(a)-1)*q)] for q in (.25,.5,.75))
    groups=defaultdict(list)
    for row in rows:
        if row[args.id_column] in known: continue
        group=tuple(binidx(row[k],vals[k]) for k in cols[1:])
        groups[group].append(row[args.id_column])
    rng=random.Random(20261031)
    for group in groups:rng.shuffle(groups[group])
    keys=sorted(groups);rng.shuffle(keys)
    picks=[]
    while len(picks)<10:
        if not any(groups.values()):ap.error('Insufficient eligible IDs')
        for key in keys:
            if groups[key] and len(picks)<10:picks.append(groups[key].pop())
    result={'source_layer':'R61_832_ONLY','sha256':hashlib.sha256(raw).hexdigest(),
            'n_rows':832,'strata_columns':list(cols[1:]),'random_seed':20261031,
            'known_10':known,'stratified_10':picks,'pilot_20':known+picks,
            'status':'AWAITING_CLAUDE_AND_USER_APPROVAL_NO_REAL_FITS'}
    if args.output.exists():ap.error('Refusing to overwrite existing output')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('FROZEN_20_FOR_REVIEW_ONLY',args.output,'SHA256',result['sha256'])
if __name__=='__main__':main()
