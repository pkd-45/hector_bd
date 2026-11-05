#!/usr/bin/env python3
"""Read-only block cross-correlation triage of WAVES image/segmap pixel registration.

Runs on one pair OR in batch from step1_segmap_audit.tsv (only WCS-mismatched pairs).
Positive shift: segmentation footprints lie +dx,+dy relative to image flux.
Verdicts are diagnostic ONLY, never automatic scientific eligibility.
Original concept: Claude 2026-11-05; batch/safety/uncertainty revisions by ChatGPT.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def block_shift(image, labels, radius):
    good=np.isfinite(image)
    if good.sum()<50:
        return None
    background=np.median(image[good])
    I=np.where(good,image,background).astype(float)
    I-=background
    # Compress bright contaminants; the segmentation is a *partition*, not a flux model.
    hi=np.percentile(I[good],99.9)
    if hi>0:
        I=np.minimum(I,hi)
    I=np.sqrt(np.clip(I,0,None))
    M=(labels>0).astype(float)
    if I.sum()<=0 or M.sum()<50:
        return None
    ny,nx=I.shape
    if radius>=min(nx,ny)//2:
        return None
    FI=np.fft.rfft2(I,s=(2*ny,2*nx));FM=np.fft.rfft2(M,s=(2*ny,2*nx))
    C=np.fft.irfft2(np.conj(FI)*FM,s=(2*ny,2*nx))
    C=np.roll(np.roll(C,radius,0),radius,1)[:2*radius+1,:2*radius+1]
    iy,ix=np.unravel_index(np.argmax(C),C.shape)
    if iy in (0,2*radius) or ix in (0,2*radius):
        return None
    median=np.median(C)
    sig=(C.max()-median)/(1.4826*np.median(np.abs(C-median))+1e-12)
    def refine(a,b,c):
        d=a-2*b+c
        return 0.0 if abs(d)<1e-12 else float(np.clip(.5*(a-c)/d,-.5,.5))
    dy=iy-radius+refine(C[iy-1,ix],C[iy,ix],C[iy+1,ix])
    dx=ix-radius+refine(C[iy,ix-1],C[iy,ix],C[iy,ix+1])
    return float(dx),float(dy),float(sig)


def analyse(image,seg,wimg,wseg,nblock=4,search=15,min_sig=3.0):
    if image.shape!=seg.shape or image.ndim!=2:
        raise ValueError('2D_GRID_SHAPE_MISMATCH')
    ny,nx=image.shape
    rows=[]
    for by in range(nblock):
        for bx in range(nblock):
            y0,y1=by*ny//nblock,(by+1)*ny//nblock
            x0,x1=bx*nx//nblock,(bx+1)*nx//nblock
            r=block_shift(image[y0:y1,x0:x1],seg[y0:y1,x0:x1],search)
            cx,cy=(x0+x1-1)/2,(y0+y1-1)/2
            ra,dec=wimg.all_pix2world([[cx,cy]],0)[0]
            px,py=wseg.all_world2pix([[ra,dec]],0)[0]
            rows.append(dict(block=f'{bx},{by}',pred_dx=float(px-cx),pred_dy=float(py-cy),
                             meas_dx=None if r is None else r[0],
                             meas_dy=None if r is None else r[1],
                             peak_sigma=None if r is None else r[2]))
    good=[r for r in rows if r['meas_dx'] is not None and r['peak_sigma']>=min_sig]
    out=dict(n_blocks=len(rows),n_usable=len(good),blocks=rows,
             verdict='AMBIGUOUS',note='Cross-correlation is diagnostic, not proof of registration.')
    if len(good)<max(4,len(rows)//2):
        out['verdict']='AMBIGUOUS_INSUFFICIENT'
        return out
    m=np.array([[r['meas_dx'],r['meas_dy']] for r in good]);p=np.array([[r['pred_dx'],r['pred_dy']] for r in good])
    pred_norm=np.hypot(p[:,0],p[:,1]);meas_norm=np.hypot(m[:,0],m[:,1])
    diff_norm=np.hypot(m[:,0]-p[:,0],m[:,1]-p[:,1])
    pred_size=float(np.median(pred_norm));to_zero=float(np.median(meas_norm));to_pred=float(np.median(diff_norm))
    out.update(median_residual_from_header_px=round(to_pred,3),
               median_measured_shift_px=round(to_zero,3),
               median_predicted_shift_px=round(pred_size,3),
               max_predicted_shift_px=round(float(np.max(pred_norm)),3),
               measured_shift_scatter_px=round(float(np.median(np.hypot(m[:,0]-np.median(m[:,0]),m[:,1]-np.median(m[:,1])))),3))
    if np.max(np.maximum(np.abs(p[:,0]),np.abs(p[:,1])))>=search-2:
        out['verdict']='AMBIGUOUS_SEARCH_RANGE'
    elif pred_size<0.3:
        out['verdict']='HEADERS_AGREE_NO_DISCRIMINATION'
    elif to_pred<0.3*pred_size and to_pred<0.35 and to_zero>max(.4,1.5*to_pred):
        out['verdict']='MISREGISTERED_CANDIDATE'
    elif to_zero<0.3*pred_size and to_zero<0.35 and to_pred>max(.4,1.5*to_zero):
        out['verdict']='REGISTERED_CANDIDATE'
    return out


def inspect_pair(imgpath,segpath,image_hdu,seg_hdu,nblock,search,min_sig):
    with fits.open(imgpath,memmap=False) as fi,fits.open(segpath,memmap=False) as fs:
        ih=fi[image_hdu];sh=fs[seg_hdu]
        if ih.data is None or sh.data is None or ih.data.ndim!=2 or sh.data.ndim!=2:
            raise ValueError('DECLARED_HDU_NOT_2D')
        image=np.asarray(ih.data,dtype=float)
        seg=np.asarray(sh.data)
        if not np.all(np.isfinite(seg)):
            raise ValueError('NONFINITE_SEGMAP')
        wimg=WCS(ih.header).celestial;wseg=WCS(sh.header).celestial
    return analyse(image,seg,wimg,wseg,nblock,search,min_sig)


def ensure_new_home_file(filename):
    path=Path(filename).expanduser()
    if path.exists() or path.is_symlink():
        raise ValueError('--out must be a new file')
    path=path.resolve(strict=False)
    home=Path.home().resolve()
    if home not in path.parents:
        raise ValueError('Output must be under HOME, outside historical /data or /scratch')
    path.parent.mkdir(parents=True,exist_ok=True)
    return path


def parse_num(value):
    try:return float(value)
    except (TypeError,ValueError):return None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('image',nargs='?');ap.add_argument('segmap',nargs='?')
    ap.add_argument('--audit-tsv',help='Batch scan audit TSV rows with WCS offset above threshold')
    ap.add_argument('--out',help='New HOME TSV output for batch mode; JSON for single-pair mode')
    ap.add_argument('--min-wcs-arcsec',type=float,default=0.05)
    ap.add_argument('--image-hdu',type=int,default=1)
    ap.add_argument('--seg-hdu',type=int,default=1)
    ap.add_argument('--search',type=int,default=15)
    ap.add_argument('--nblock',type=int,default=4)
    ap.add_argument('--min-sig',type=float,default=3)
    a=ap.parse_args()
    if a.image_hdu<0 or a.seg_hdu<0 or a.search<2 or not 2<=a.nblock<=8:
        ap.error('Invalid HDU, search radius or number of blocks')
    if a.audit_tsv:
        if a.image or a.segmap or not a.out:
            ap.error('Batch mode requires --audit-tsv and --out, but no positional files')
        output=ensure_new_home_file(a.out)
        with open(a.audit_tsv,encoding='utf-8-sig') as f:
            entries=list(csv.DictReader(f,delimiter='\t'))
        candidates=[r for r in entries if (parse_num(r.get('wcs_max_sep_arcsec')) or 0)>a.min_wcs_arcsec]
        results=[];detailpath=output.with_suffix('.blocks.jsonl')
        with open(detailpath,'x') as detail:
            for i,r in enumerate(candidates,1):
                item=dict(uberID=r['uberID'],image_path=r.get('image_path'),segmap_path=r.get('segmap_path'),
                          wcs_max_sep_arcsec=r.get('wcs_max_sep_arcsec'))
                try:
                    v=inspect_pair(item['image_path'],item['segmap_path'],a.image_hdu,a.seg_hdu,a.nblock,a.search,a.min_sig)
                    blocks=v.pop('blocks')
                    item.update(v)
                    detail.write(json.dumps(dict(uberID=r['uberID'],blocks=blocks))+'\n')
                except Exception as e:
                    item.update(verdict='ERROR',reason=f'{type(e).__name__}: {str(e)[:220]}')
                results.append(item)
                print(f'[{i}/{len(candidates)}] {r["uberID"]} {item["verdict"]} usable={item.get("n_usable","")}',flush=True)
        columns=('uberID','verdict','reason','wcs_max_sep_arcsec','n_blocks','n_usable',
                 'median_predicted_shift_px','median_measured_shift_px','median_residual_from_header_px',
                 'measured_shift_scatter_px','max_predicted_shift_px','image_path','segmap_path')
        with open(output,'x',newline='') as f:
            w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',extrasaction='ignore')
            w.writeheader();w.writerows(results)
        counts={}
        for r in results: counts[r['verdict']]=counts.get(r['verdict'],0)+1
        info=dict(rows_checked=len(results),threshold_arcsec=a.min_wcs_arcsec,counts=counts,
                  eligible_for_fitting=False,note='Candidate-only cross-correlation diagnostic')
        output.with_suffix('.summary.json').write_text(json.dumps(info,indent=2)+'\n')
        print('OUTPUT_TSV='+str(output));print(json.dumps(info,indent=2))
    else:
        if not a.image or not a.segmap:
            ap.error('Provide IMAGE SEGMAP or --audit-tsv')
        v=inspect_pair(a.image,a.segmap,a.image_hdu,a.seg_hdu,a.nblock,a.search,a.min_sig)
        print(json.dumps(v,indent=2))
        if a.out:
            output=ensure_new_home_file(a.out)
            output.write_text(json.dumps(v,indent=2)+'\n')
    return 0

if __name__=='__main__':
    sys.exit(main())
