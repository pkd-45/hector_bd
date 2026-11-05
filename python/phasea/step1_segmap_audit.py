#!/usr/bin/env python3
"""Read-only WAVES 832-roster delivered-vs-historical segmentation audit.

No fit, source modification, resegmentation or catalog-label invention.
Use --limit/--ids to preflight, then run the exact frozen 832 roster.
Original concept: Claude 2026-11-05; fail-closed revisions and checks by ChatGPT.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
try:
    from scipy import ndimage
except ImportError:
    ndimage = None

ROOT_DEFAULT = '/data/scratch/uqpdas/WAVES_bd_v32_production842/baseline/FIXED4_seed666'
ROSTER_DEFAULT = '/scratch/uqpdas/WAVES_BD_PROFIT_PROFUSE_R61_R53_NATIVE_SPD_FIXED4_POLICY_EXECUTION_20260930/results/R61_NATIVE_MODEL_SELECTION_832.tsv'
ID_FIELDS = ('uberID','uberid','uber_id','UBERID','id','ID','galaxy','galaxy_id','WAVES_ID','WAVESID')
FIELDS = ('uberID','status','reason','hist_segID','seed','image_path','segmap_path',
          'source_inputs_rds','image_hdu','seg_hdu','shape','orientation','orientation_score_active',
          'orientation_ambiguous','n_target_delivered','n_target_historical',
          'removed_from_target','added_to_target','added_from_sky','added_from_neighbours',
          'other_label_changes','islands_delivered','islands_historical',
          'wcs_max_sep_arcsec','wcs_centre_sep_arcsec','cat_pixel','cat_label',
          'cat_boundary_flag','cat_agrees_with_historical','cat_label_is_sky',
          'catalogue_row_found','seg_product_status')

R_DUMP = r'''a<-commandArgs(TRUE);x<-readRDS(a[1]);s<-x$segim
if(is.null(dim(s)) || length(dim(s)) != 2L || any(!is.finite(s)) ||
    any(s < 0) || any(s != floor(s)) || any(s > .Machine$integer.max))
  stop("BAD_SAVED_SEGMAP")
storage.mode(s)<-"integer"
writeBin(as.vector(s),a[2],size=4L,endian="little")
g<-function(v) if(is.null(v)||!length(v)) "" else as.character(v[[1]])
cat("BD_META\t",nrow(s),"\t",ncol(s),"\t",g(x$target$segID),"\t",
  g(x$target$seed_x),"\t",g(x$target$seed_y),"\n",sep="")
cat("BD_IMAGE\t",g(x$image_path),"\n",sep="")
cat("BD_SEGMAP\t",g(x$segmap_path),"\n",sep="")'''

def sha256(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda: f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()

def identifier(x):
    a = str(x or '').strip()
    if a.startswith('W'):
        a = a[1:]
    if not re.fullmatch(r'[0-9]+',a):
        raise ValueError('Invalid uberID: '+repr(x))
    return a

def roster_ids(path):
    with open(path,encoding='utf-8-sig',newline='') as f:
        rows = list(csv.DictReader(f,delimiter='\t'))
    if not rows:
        raise ValueError('Empty roster')
    available = [k for k in ID_FIELDS if k in rows[0]]
    if not available:
        raise ValueError('No recognised ID field in roster: '+','.join(rows[0]))
    key = available[0]
    ids = [identifier(r[key]) for r in rows]
    if len(ids)!=832 or len(set(ids))!=832:
        raise ValueError(f'FROZEN_R61_832_GATE_FAIL count={len(ids)} unique={len(set(ids))}')
    return ids,key

def output_dir(path):
    home=Path.home().resolve()
    out=Path(path).expanduser()
    if out.exists() or out.is_symlink():
        raise ValueError('--out must be a new, nonexistent output directory')
    out=out.resolve(strict=False)
    if home==out or home not in out.parents:
        raise ValueError('Output must be a new directory under your HOME, never /data or /scratch')
    for ancestor in (out, *out.parents):
        if ancestor.exists() and ancestor.is_symlink():
            raise ValueError('Symlink path component refused: '+str(ancestor))
    out.mkdir(parents=True,exist_ok=False)
    return out

def load_catalogue(path,id_col,ra_col,dec_col):
    if not path:
        return None
    if not all((id_col,ra_col,dec_col)):
        raise ValueError('Catalogue requires explicit --id-col, --ra-col, --dec-col')
    res={}
    with open(path,encoding='utf-8-sig',newline='') as f:
        for row in csv.DictReader(f):
            try:
                gid=identifier(row[id_col]); ra=float(row[ra_col]); dec=float(row[dec_col])
            except (ValueError,TypeError,KeyError) as e:
                raise ValueError('Invalid catalogue row: '+str(e))
            if gid in res:
                raise ValueError('Duplicate catalogue uberID '+gid)
            if not(np.isfinite(ra) and np.isfinite(dec) and -90<=dec<=90):
                raise ValueError('Nonfinite/invalid catalogue coordinates '+gid)
            res[gid]=(ra,dec)
    return res

def dump_inputs(rscript,rds,tmp,timeout):
    dst=Path(tmp)/'segim.bin'
    proc=subprocess.run([rscript,'--vanilla','-e',R_DUMP,str(rds),str(dst)],
                        capture_output=True,text=True,timeout=timeout)
    if proc.returncode:
        raise RuntimeError('RDS_READ_FAILED '+proc.stderr[-300:])
    dat={}
    for line in proc.stdout.splitlines():
        if line.startswith('BD_'):
            parts=line.split('\t')
            dat[parts[0]]=parts[1:]
    if not all(k in dat for k in ('BD_META','BD_IMAGE','BD_SEGMAP')):
        raise ValueError('RDS_DUMP_MISSING_METADATA '+proc.stdout[-300:])
    nr,nc,segid,sx,sy=dat['BD_META'][:5]
    nrow,ncol=int(nr),int(nc)
    if nrow<=0 or ncol<=0 or nrow*ncol>100000000:
        raise ValueError('BAD_RDS_DIMENSIONS')
    if dst.stat().st_size != nrow*ncol*4:
        raise ValueError('BAD_RDS_DUMP_SIZE')
    saved=np.fromfile(dst,dtype='<i4').reshape((nrow,ncol),order='F')
    dst.unlink()
    return saved,segid,sx,sy,dat['BD_IMAGE'][0],dat['BD_SEGMAP'][0]

def wcs_sep(imhead,seghead,shape):
    wi=WCS(imhead).celestial; ws=WCS(seghead).celestial
    ny,nx=shape
    xx=np.array([0.,nx-1.,0.,nx-1.,(nx-1.)/2])
    yy=np.array([0.,0.,ny-1.,ny-1.,(ny-1.)/2])
    s=wi.pixel_to_world(xx,yy).separation(ws.pixel_to_world(xx,yy)).arcsec
    if not np.all(np.isfinite(s)):
        raise ValueError('NONFINITE_WCS_SEPARATION')
    return round(float(np.max(s)),6),round(float(s[-1]),6)

def islands(mask):
    if ndimage is None:
        return None
    return int(ndimage.label(mask,structure=np.ones((3,3),dtype=int))[1])

def orient_saved(saved,original):
    """Score only pixels containing a label; never let blank sky dominate."""
    candidates={'direct':saved,'transpose':saved.T}
    candidates={k:v for k,v in candidates.items() if v.shape==original.shape}
    if not candidates:
        raise ValueError(f'SHAPE_MISMATCH delivered={original.shape} saved={saved.shape}')
    scores={}
    for key,v in candidates.items():
        active=(v!=0)|(original!=0)
        scores[key]=float(np.mean(v[active]==original[active])) if np.any(active) else 1.0
    best=max(scores,key=scores.get)
    ambiguous=(len(scores)>1 and abs(scores['direct']-scores['transpose'])<1e-8 and
               not np.array_equal(candidates['direct'],candidates['transpose']))
    return candidates[best],best,round(scores[best],8),ambiguous

def audit_one(gid,args,tmp,catalogue):
    out=dict(uberID=gid,seg_product_status='AS_DELIVERED_SUBTYPE_UNVERIFIED')
    rds=Path(args.root)/gid/args.subdir/'inputs.rds'
    out['source_inputs_rds']=str(rds)
    if not rds.is_file():
        return dict(out,status='NO_INPUTS_RDS')
    saved,segid,sx,sy,imgpath,segpath=dump_inputs(args.rscript,rds,tmp,args.r_timeout)
    out.update(hist_segID=segid,seed=f'{sx},{sy}',image_path=imgpath,segmap_path=segpath,
               image_hdu=args.image_hdu,seg_hdu=args.seg_hdu)
    if not Path(imgpath).is_file() or not Path(segpath).is_file():
        return dict(out,status='SOURCE_FITS_MISSING')
    with fits.open(segpath,memmap=False) as fs, fits.open(imgpath,memmap=False) as fi:
        shdu=fs[args.seg_hdu]; ihdu=fi[args.image_hdu]
        if shdu.data is None or ihdu.data is None or shdu.data.ndim!=2 or ihdu.data.ndim!=2:
            return dict(out,status='NOT_2D_DECLARED_HDU')
        L=np.asarray(shdu.data)
        ishape=ihdu.data.shape
        out['shape']=f'{L.shape[0]}x{L.shape[1]}'
        if L.shape!=ishape:
            return dict(out,status='IMAGE_SEG_GRID_SHAPE_MISMATCH')
        if not np.issubdtype(L.dtype,np.number) or not np.all(np.isfinite(L)) or np.any(L<0) or np.any(L!=np.floor(L)):
            return dict(out,status='INVALID_ORIGINAL_LABELS')
        hseg=shdu.header.copy(); himg=ihdu.header.copy()
    L=L.astype(np.int64,copy=False)
    out['wcs_max_sep_arcsec'],out['wcs_centre_sep_arcsec']=wcs_sep(himg,hseg,L.shape)
    S,orient,score,ambig=orient_saved(saved,L)
    out.update(orientation=orient,orientation_score_active=score,orientation_ambiguous=ambig)
    if ambig:
        return dict(out,status='ORIENTATION_AMBIGUOUS')
    try:
        t=int(segid)
        if t<=0 or str(float(segid))=='nan':
            raise ValueError()
    except (ValueError,TypeError):
        return dict(out,status='NO_HISTORICAL_SEGID')
    target=L==t; hist=S==t
    rem=int(np.count_nonzero(target&~hist)); add=int(np.count_nonzero(hist&~target))
    oth=int(np.count_nonzero((L!=S)&~target&~hist))
    out.update(n_target_delivered=int(target.sum()),n_target_historical=int(hist.sum()),
               removed_from_target=rem,added_to_target=add,
               added_from_sky=int(np.count_nonzero(hist&(L==0))),
               added_from_neighbours=int(np.count_nonzero(hist&(L>0)&~target)),
               other_label_changes=oth,islands_delivered=islands(target),islands_historical=islands(hist))
    if not np.any(target):
        out['status']='HIST_LABEL_MISSING_IN_DELIVERED'
    elif rem==0 and add==0 and oth==0:
        out['status']='IDENTICAL'
    elif add==0 and oth==0:
        out['status']='TARGET_REMOVED_ONLY'
    elif rem==0 and oth==0:
        out['status']='TARGET_ADDED_ONLY'
    else:
        out['status']='MIXED_OR_OTHER_LABEL_CHANGES'
    if catalogue is not None:
        out['catalogue_row_found']=gid in catalogue
        if gid in catalogue:
            ra,dec=catalogue[gid]
            x,y=WCS(hseg).celestial.world_to_pixel_values(ra,dec)
            if not (np.isfinite(x) and np.isfinite(y)):
                out['cat_label']='NONFINITE_WCS'
            else:
                xi=int(np.floor(x+0.5)); yi=int(np.floor(y+0.5))
                out['cat_pixel']=f'{xi},{yi}'
                if 1<=yi<L.shape[0]-1 and 1<=xi<L.shape[1]-1:
                    lab=int(L[yi,xi]); nb=np.unique(L[yi-1:yi+2,xi-1:xi+2])
                    out.update(cat_label=lab,cat_boundary_flag=bool(len(nb)>1),
                               cat_agrees_with_historical=bool(lab==t),cat_label_is_sky=bool(lab==0))
                else:
                    out['cat_label']='EDGE'
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--roster',default=ROSTER_DEFAULT)
    p.add_argument('--root',default=ROOT_DEFAULT)
    p.add_argument('--subdir',default='r_1333_v25')
    p.add_argument('--image-hdu',type=int,required=True)
    p.add_argument('--seg-hdu',type=int,required=True)
    p.add_argument('--rscript',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--limit',type=int,help='Smoke-test first N IDs from the exact 832 roster')
    p.add_argument('--ids',help='Comma-separated subset, must all belong to frozen 832')
    p.add_argument('--catalogue',help='Independent sample catalogue CSV; optional, unresolved if omitted')
    p.add_argument('--id-col');p.add_argument('--ra-col');p.add_argument('--dec-col')
    p.add_argument('--r-timeout',type=int,default=120)
    a=p.parse_args()
    if a.image_hdu<0 or a.seg_hdu<0 or (a.limit is not None and a.limit<1) or a.r_timeout<1:
        p.error('HDU indices >=0, limit >=1 and timeout >=1 required')
    if not Path(a.rscript).is_file():
        p.error('Rscript executable path not found')
    ids,id_column=roster_ids(a.roster)
    roster_sha=sha256(a.roster)
    if a.ids:
        want=[identifier(x) for x in a.ids.split(',')]
        if len(want)!=len(set(want)) or any(x not in ids for x in want):
            p.error('Specified --ids must be unique and in the frozen 832 roster')
        ids=want
    elif a.limit:
        ids=ids[:a.limit]
    cat=load_catalogue(a.catalogue,a.id_col,a.ra_col,a.dec_col)
    out=output_dir(a.out)
    script_sha=sha256(__file__)
    rows=[]
    with tempfile.TemporaryDirectory(dir=out,prefix='scratch_') as tmp, \
        open(out/'step1_segmap_audit.tsv','w',newline='') as dest:
        writer=csv.DictWriter(dest,fieldnames=FIELDS,delimiter='\t',extrasaction='ignore')
        writer.writeheader()
        for i,gid in enumerate(ids,1):
            try:
                r=audit_one(gid,a,tmp,cat)
            except Exception as e:
                r=dict(uberID=gid,status='ERROR',reason=f'{type(e).__name__}: {str(e)[:250]}')
            writer.writerow(r)
            dest.flush()
            rows.append(r)
            print(f'[{i}/{len(ids)}] {gid} {r.get("status")} removed={r.get("removed_from_target","")} added={r.get("added_to_target","")} WCS={r.get("wcs_max_sep_arcsec","")}',flush=True)
    counts={}
    for r in rows:
        counts[r['status']]=counts.get(r['status'],0)+1
    summary=dict(frozen_roster_count=832,roster_id_col=id_column,roster_sha256=roster_sha,
                 script_sha256=script_sha,run_count=len(rows),status_counts=counts,
                 WCS_gt_0p05_arcsec=sum(float(r.get('wcs_max_sep_arcsec') or 0)>0.05 for r in rows),
                 has_independent_catalogue=bool(cat is not None),
                 catalogue_agrees=sum(r.get('cat_agrees_with_historical') is True for r in rows),
                 catalogue_disagrees=sum(r.get('cat_agrees_with_historical') is False for r in rows),
                 eligible_for_fitting=False,
                 note='This audits saved vs delivered maps; it does NOT certify catalogue provenance, segmentation product subtype or observed G0/G2.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'SHA256SUMS.txt').write_text(f'{sha256(out/"step1_segmap_audit.tsv")}  step1_segmap_audit.tsv\n{sha256(out/"summary.json")}  summary.json\n')
    print('OUTPUT_DIR='+str(out));print(json.dumps(summary,indent=2))
    return 0

if __name__=='__main__':
    sys.exit(main())
