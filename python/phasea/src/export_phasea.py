#!/usr/bin/env python3
"""D1 historical pixel-lineage audit. No fitting or mutation of inputs.

SOURCE REVIEW (verified): The contract is original WAVES L, target label set A_g,
T_g=(L in A_g); Data$region is only evidence of historical likelihood support.

ASSUMPTIONS/UNRESOLVED, deliberately not promoted into facts:
- Original image/segmentation HDUs MUST be explicitly supplied; first/primary is NOT assumed.
- Catalogue-to-label association MUST be explicitly supplied; proximity or centre labels
  are NEVER used to invent A_g. Multi-label strings are supported, not inferred.
- Original FITS pixels use numpy [row=y, col=x] coordinates, zero-based. R arrays
  are column-major; their orientation is recovered from the saved pixel data by
  direct/transposed image registration. No silent flips or rotations.
- FITS coordinate consistency across image and segmap is not established by shape
  alone: WCS identity is recorded and provenance remains an external Step 0 gate.
- Sky-subtraction and/or unit scaling are allowed in correlation registration;
  only unique, strongly validated image matches are accepted.
- An observed pixel omission is not proof of the historical branch responsible.
  Cause categories are evidence-ranked, and inferred causes are labelled INFERRED.
- It is not possible to define F_g until the original image, sigma and detector
  mask provenance are established; no false G2 pass is issued here.
"""
import argparse, csv, gzip, hashlib, json, math, os, re, shutil, subprocess, sys, traceback
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

PILOTS = ['42710765002459','221990044410325','351690654304280',
          '176980350403709','12690079611187','41691060304604',
          '350680203312123','18680557409852','49710254502112',
          '213991001205426']
KNOWN_ROOT = '/data/scratch/uqpdas/WAVES_bd_v32_production842/baseline/FIXED4_seed666'
KNOWN_R61 = '/scratch/uqpdas/WAVES_BD_PROFIT_PROFUSE_R61_R53_NATIVE_SPD_FIXED4_POLICY_EXECUTION_20260930/results/R61_NATIVE_MODEL_SELECTION_832.tsv'
FIELDS = ['id','inputs_path','status','errors','label_provenance','labels',
 'source_image_path','segmap_path','image_hdu','seg_hdu','original_grid','source_edge',
 'segmentation_product','source_image_sha256','segmap_sha256','image_wcs_hash','segmap_wcs_hash',
 'N_target','N_region','region_original_coords_hash','N_missing','N_extra','N_extra_original_background','N_extra_original_neighbour','N_missing_cutbox_side_fit',
 'N_missing_island_filter','N_missing_autoclip_NA','N_missing_detector_mask',
 'N_missing_invalid_sigma','N_missing_source_invalid','N_missing_segim_differs_upstream','N_missing_unexplained','N_calcregion','N_calc_only',
 'N_calc_outside_stamp','N_masked_target_in_stamp','N_invalid_sigma_target_in_stamp',
 'offset_row0','offset_col0','orientation','offset_method','offset_validation',
 'offset_uniqueness','input_folder','input_folders_found','found2fit_domain_provenance','profound_fields_json','data_seg_partition','profound_seg_partition','data_seg_partition_mismatch_pixels','profound_seg_partition_mismatch_pixels','stamp_exact_copy','stamp_max_abs_original_delta','stamp_scale','stamp_shift','selected_fit_source','fit_binding','fits_available','found2fit_branch',
 'matches_target_label','matches_any_individual_label','matches_union_original_labels','overlay_status','psf_kind','psf_d_use',
 'psf_fwhm','psf_beta','historical_finesample_by_fit_json','historical_psftype_by_fit_json','psf_sum','psf_ee_fraction','psf_ee_basis','psf_hash','fit_data_psf_hash','psf_matches_fit_data',
 'sigma_source','sky_source','gain_source','runtime_id','family_results_json','provenance']


def sha256_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def normalize_id(s):
    s=str(s).strip()
    if s.upper().startswith('W'):s=s[1:]
    return s


def read_tsv(path):
    with open(path, newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def id_field(row):
    for name in ('uberID','uberid','uber_id','UBERID','id','ID','galaxy','galaxy_id','WAVES_ID','WAVESID'):
        if name in row and row[name] is not None and str(row[name]).strip():return normalize_id(row[name])
    raise ValueError('No recognised ID field; require uberID/id/galaxy in frozen 832-row table')


def sample_ids(path):
    rows=read_tsv(path)
    ids=[id_field(r) for r in rows]
    if len(rows)!=832 or len(set(ids))!=832 or any(not re.fullmatch('[0-9]+',i) for i in ids):
        raise ValueError(f'FROZEN_SAMPLE_GATE_FAIL rows={len(rows)} unique={len(set(ids))}: never mix 842/832')
    return ids


def mappings(path):
    if not path:return {}
    result={}
    for r in read_tsv(path):
        g=id_field(r)
        label_s=next((r[x] for x in ('labels','label_ids','seg_labels','segid','segID','segment_labels') if x in r),None)
        if label_s is None:raise ValueError('label-map requires labels column (comma/semicolon-separated integer labels)')
        labels=[int(x.strip()) for x in re.split('[,;| ]+',label_s.strip()) if x.strip()]
        if not labels or any(x<=0 for x in labels):raise ValueError(f'Invalid nonpositive/empty label set for {g}')
        if g in result:raise ValueError(f'Duplicate galaxy in label-map {g}')
        result[g]=sorted(set(labels))
    return result


def path_map(path):
    if not path:return {}
    out={}
    for row in read_tsv(path):
        a=row.get('recorded_path');b=row.get('verified_actual_path')
        if not a or not b or not os.path.isabs(a) or not os.path.isabs(b):
            raise ValueError('Path map requires absolute recorded_path and verified_actual_path')
        if a in out:raise ValueError('Duplicate recorded path '+a)
        out[a]=b
    return out


def resolve_recorded(src,explicit):
    if not src or not os.path.isabs(src):return '', 'INVALID_RECORDED_PATH'
    if src in explicit:
        actual=explicit[src]
        return actual, 'EXPLICIT_USER_VERIFIED_PATH_MAP' if os.path.isfile(actual) else 'EXPLICIT_MAP_FILE_MISSING'
    return src, 'EXACT_RECORDED_PATH' if os.path.isfile(src) else 'ORIGINAL_NOT_FOUND_NO_GUESS'


def parse_metadata(path):
    out={}
    if not path.exists():return out
    for r in read_tsv(path):out.setdefault(r['source'],{})[r['key']]=r['value']
    return out


def parse_arrays(path,folder):
    out={}
    if not path.exists():return out
    for r in read_tsv(path):
        n,m=int(r['nrow']),int(r['ncol'])
        p=folder/r['file']
        if p.stat().st_size!=8*n*m:raise ValueError('array byte length differs from declared dimension: '+str(p))
        arr=np.fromfile(p,dtype='<f8').reshape((n,m),order='F')
        out[(r['source'],r['key'])]=arr
    return out


def find_field(meta,source,keys):
    obj=meta.get(source,{})
    for key in keys:
        if key in obj and obj[key].strip():return obj[key]
    return ''


def wcs_fingerprint(header):
    keys=[k for k in header if k.startswith(('CRPIX','CRVAL','CTYPE','CDELT','CD','PC','CUNIT','WCSAXES'))]
    payload='\n'.join(f'{k}={header[k]}' for k in sorted(set(keys)))
    return hashlib.sha256(payload.encode()).hexdigest() if payload else ''

@lru_cache(maxsize=5)
def fits_array(path,hdu):
    from astropy.io import fits
    with fits.open(path,mode='readonly',memmap=False) as hdus:
        obj=hdus[int(hdu)] if re.fullmatch(r'\d+',str(hdu)) else hdus[hdu]
        if obj.data is None:raise ValueError('selected HDU contains no image array')
        arr=np.asarray(obj.data)
        if arr.ndim!=2:raise ValueError(f'expected 2D FITS image, got {arr.shape}')
        return np.asarray(arr), wcs_fingerprint(obj.header)


def robust_verification(src,stamp,row,col):
    h,w=stamp.shape
    if row<0 or col<0 or row+h>src.shape[0] or col+w>src.shape[1]:return None
    a=np.asarray(src[row:row+h,col:col+w],dtype=float)
    b=np.asarray(stamp,dtype=float)
    use=np.isfinite(a)&np.isfinite(b)
    if np.count_nonzero(use)<max(12,int(.40*h*w)):return None
    x=a[use];y=b[use]
    xm=x.mean();ym=y.mean();xc=x-xm;yc=y-ym
    normx=np.linalg.norm(xc);normy=np.linalg.norm(yc)
    if normx<1e-12 or normy<1e-12:return None
    corr=float(np.dot(xc,yc)/normx/normy)
    scale=float(np.dot(xc,yc)/np.dot(xc,xc)); shift=float(ym-scale*xm)
    err=y-(scale*x+shift)
    residual=float(np.linalg.norm(err)/normy)
    return dict(corr=corr,residual=residual,scale=scale,sky_offset=shift,valid=int(len(x)))


def best_offset(src,stamp):
    """Uniquely validated matching with free multiplicative/additive calibration.
    cv2 fast path; scipy FFT fallback; never accept a 'best guess' below tolerance.
    """
    if stamp is None or src is None:return None
    candidates=[]
    for orientation in ('direct','transpose'):
        t=stamp if orientation=='direct' else stamp.T
        if t.shape[0]>src.shape[0] or t.shape[1]>src.shape[1]:continue
        if t.shape==src.shape:
            v=robust_verification(src,t,0,0)
            if v:candidates.append((v['corr'],0,0,orientation,v));continue
        finite=np.isfinite(t)
        if finite.sum()<max(12,int(.40*t.size)):continue
        med=float(np.nanmedian(t));tfill=np.where(finite,t,med).astype('float32')
        srcfill=np.where(np.isfinite(src),src,np.nanmedian(src)).astype('float32')
        try:
            import cv2
            s=cv2.matchTemplate(srcfill,tfill,cv2.TM_CCOEFF_NORMED)
        except ImportError:
            # Numpy-only correlation fallback for the existing Friday Conda env.
            # Avoid introducing a SciPy dependency into an otherwise read-only audit.
            def fftconv_valid(x,k):
                H,W=x.shape;h,w=k.shape
                dims=(H+h-1,W+w-1)
                z=np.fft.irfftn(np.fft.rfftn(x,s=dims)*np.fft.rfftn(k,s=dims),s=dims)
                return z[h-1:H,w-1:W]
            tcenter=tfill-tfill.mean()
            numerator=fftconv_valid(srcfill.astype(float),tcenter[::-1,::-1])
            ones=np.ones(t.shape,dtype=float)
            sumsq=fftconv_valid(srcfill.astype(float)**2,ones)
            sumx=fftconv_valid(srcfill.astype(float),ones)
            ss=np.maximum(0,sumsq-sumx*sumx/t.size)
            s=numerator/np.maximum(1e-12,np.sqrt(ss*np.sum(tcenter*tcenter)))
        s=np.where(np.isfinite(s),s,-np.inf)
        if not np.any(np.isfinite(s)):continue
        top=min(s.size,12)
        flat=np.argpartition(s.ravel(),-top)[-top:]
        for ix in flat:
            r,c=np.unravel_index(ix,s.shape)
            v=robust_verification(src,t,int(r),int(c))
            if v:candidates.append((v['corr'],int(r),int(c),orientation,v))
    candidates.sort(key=lambda x:(x[0],-x[4]['residual']),reverse=True)
    if not candidates:return None
    valid=[c for c in candidates if c[4]['corr']>=.995 and c[4]['residual']<=.10 and c[4]['scale']>0]
    if not valid:return dict(status='UNVERIFIED_IMAGE_REGISTRATION',best=candidates[0])
    unique={ (c[1],c[2],c[3]):c for c in valid }
    best=valid[0]
    competing=[c for c in unique.values() if (c[1],c[2],c[3])!=(best[1],best[2],best[3]) and c[0]>=best[0]-.0001]
    if competing:return dict(status='AMBIGUOUS_IMAGE_REGISTRATION',best=best,n_alternatives=len(competing))
    return dict(status='VERIFIED_PIXEL_REGISTRATION',best=best,n_alternatives=0)


def orient(arr,orientation):return arr if orientation=='direct' else arr.T


def psf_details(meta,arrays,selected_source):
    pm=meta.get('psf_rds',{})
    kind=pm.get('output_kind','')
    if not kind:kind='UNRESOLVED_KIND'
    a=arrays.get(('psf_rds','psf'))
    if a is None:a=arrays.get((selected_source,'psf'))
    combined={**pm,**meta.get('psf_rds',{})}
    # Prefer PSF attributes read directly from the saved, selected kernel.
    def scalar_pattern(pattern):
        for k,v in combined.items():
            if re.search(pattern,k,re.IGNORECASE):
                try:
                    f=float(v)
                    if np.isfinite(f):return f,k
                except (ValueError,TypeError):pass
        return None,''
    def attr_numeric(key):
        try:
            v=float(pm.get('attr.psf.'+key,''))
            return v if np.isfinite(v) else None
        except (ValueError,TypeError):return None
    fwhm=attr_numeric('fwhm')
    if fwhm is None:fwhm,fk=scalar_pattern(r'(^|\.)fwhm(_|\.|$)')
    beta=attr_numeric('beta')
    if beta is None:beta=attr_numeric('con')
    if beta is None:beta,bk=scalar_pattern(r'(^|\.)beta($|\.)')
    if beta is None:beta,bk=scalar_pattern(r'(^|\.)con($|\.)')
    d_use=attr_numeric('size_px')
    if d_use is None:d_use,dkey=scalar_pattern(r'(^|\.)d_use($|\.)')
    if d_use is None:d_use,dkey=scalar_pattern(r'(^|\.)size_px($|\.)')
    if a is not None:
        d_use=d_use or min(a.shape)
        psfhash=hashlib.sha256(np.asarray(a,dtype='<f8').tobytes()).hexdigest()
        psfsum=float(np.nansum(a))
    else:psfhash='';psfsum=None
    # A cropped, already renormalised PSF DOES NOT reveal lost pre-crop energy.
    # The enclosed energy is unknown unless an upstream full analytic profile is provided.
    ee=None;ee_basis='UNRESOLVED_NO_UNCROPPED_PROFILE'
    stored_ee,ee_key=scalar_pattern(r'(^|\.)(encircled_energy|ee_fraction)($|\.)')
    if stored_ee is not None and 0<=stored_ee<=1:
        ee=stored_ee;ee_basis='REPORTED_SAVED_FIELD:'+ee_key
    return dict(psf_kind=kind,psf_d_use=d_use,psf_fwhm=fwhm,psf_beta=beta,
                psf_sum=psfsum,psf_ee_fraction=ee,psf_ee_basis=ee_basis,psf_hash=psfhash)


def orient_inputs_image(arrays, input_seg, input_registration):
    """Orient a saved R input image using known aligned segmentation when available."""
    raw=arrays.get(('inputs','image'))
    if raw is None or input_registration is None:return None
    if len(input_registration)>4:
        z=orient(raw,input_registration[4])
        return z if z.shape==tuple(input_registration[2:4]) else None
    s=arrays.get(('inputs','segim'))
    if s is not None and input_seg is not None:
        for ori in ('direct','transpose'):
            if np.array_equal(orient(s,ori),input_seg,equal_nan=True):
                z=orient(raw,ori)
                if z.shape==input_seg.shape:return z
    if raw.shape==tuple(input_registration[2:4]):return raw
    return None


def partition_comparison(original, saved, ignored=None):
    """Compare partitions invariant under one-to-one relabelling, including background.

    A conflict is any pair of pixels with one shared label in either partition
    but differing labels in the other. This flags genuine splits/merges, while
    allowing arbitrary bijective renumbering. Pixel mismatch count is a
    mode-map diagnostic, *not* a globally minimal edit distance.
    """
    orig=np.asarray(original);new=np.asarray(saved)
    if orig.shape!=new.shape:return dict(status='GRID_MISMATCH',n_valid=0,mismatch_pixels=None,map_new_to_orig={})
    valid=np.isfinite(orig)&np.isfinite(new)
    if ignored is not None:
        if np.shape(ignored)!=orig.shape:raise ValueError('partition ignore grid mismatch')
        valid &= ~np.asarray(ignored,dtype=bool)
    if not np.any(valid):return dict(status='NO_COMPARABLE_PIXELS',n_valid=0,mismatch_pixels=None,map_new_to_orig={})
    a=orig[valid].astype(np.int64);b=new[valid].astype(np.int64)
    # Require a bijection for complete partition identity, with 0 (background)
    # fixed to itself. A relabelled source cannot silently relabel background.
    ab={};ba={};nconflict=0
    for k in np.unique(b):
        vals,cnt=np.unique(a[b==k],return_counts=True)
        chosen=int(vals[np.argmax(cnt)])
        ab[int(k)]=chosen
        nconflict+=int(np.sum(cnt)-np.max(cnt))
    for k in np.unique(a):
        vals,cnt=np.unique(b[a==k],return_counts=True)
        ba[int(k)]=int(vals[np.argmax(cnt)])
    coherent=all(np.unique(a[b==k]).size==1 for k in np.unique(b)) and all(np.unique(b[a==k]).size==1 for k in np.unique(a))
    coherent=bool(coherent and ab.get(0,0)==0 and ba.get(0,0)==0)
    identical=bool(np.array_equal(a,b))
    status='PARTITION_IDENTICAL' if identical else ('PARTITION_IDENTICAL_RELABELED' if coherent else 'PARTITION_DIFFERS')
    # Count modal-map discordance in both directions (split and merge sensitive);
    # not an exact per-pixel reconstruction of a changed historical segmentation.
    keys=np.array(sorted(ab),dtype=np.int64)
    expected=np.array([ab[int(k)] for k in keys],dtype=np.int64)
    mismatch_ab=int(np.count_nonzero(expected[np.searchsorted(keys,b)]!=a))
    reverse_keys=np.array(sorted(ba),dtype=np.int64)
    reverse_expected=np.array([ba[int(k)] for k in reverse_keys],dtype=np.int64)
    mismatch_ba=int(np.count_nonzero(reverse_expected[np.searchsorted(reverse_keys,a)]!=b))
    mismatch=max(mismatch_ab,mismatch_ba)
    return dict(status=status,n_valid=int(valid.sum()),mismatch_pixels=mismatch,
                map_new_to_orig=ab, n_discordant_under_modal_map=mismatch)


def classify_pixel(y,x,box,inputs_seg,inputs_registration,data_arrays,fit_source,orient_name,labels,
                   original_image=None,passed_image=None,passed_seg=None,passed_registration=None,
                   passed_mask=None,passed_sigma=None,pass_label_mapping=None, fitted_seg_labels=None):
    """Use Found2Fit input domain FIRST; NEVER equate a tightcrop stamp with cutbox.

    When profound is absent, registered inputs.rds is a proxy, explicitly
    INFERRED rather than verified to be exactly what Found2Fit received.
    """
    y0,x0,h,w=box
    dom=passed_registration if passed_registration is not None else inputs_registration
    seg=passed_seg if passed_seg is not None else inputs_seg
    im=passed_image if passed_image is not None else orient_inputs_image(data_arrays,inputs_seg,inputs_registration)
    msk=passed_mask; sig=passed_sigma
    confidence='VERIFIED_SAVED_PROFOUND_ARRAYS' if passed_registration is not None else 'INFERRED_REGISTERED_INPUTS_RDS_PROXY'
    if dom is not None:
        a,b,H,W=dom[:4]
        if not (a<=y<a+H and b<=x<b+W):
            return 'cutbox_side_fit_crop',confidence+'_OUTSIDE_PRE_TIGHTCROP_INPUT_DOMAIN'
        py,px=y-a,x-b
        if seg is not None and seg.shape==(H,W) and np.isfinite(seg[py,px]):
            sv=int(seg[py,px])
            # R9: classify in the segmentation frame USED BY Found2Fit.
            # Never majority-map a split or reassigned label back to WAVES:
            # this conceals real missing pixels. The fitted label set is
            # observed directly under Data$region (not inferred from L).
            if fitted_seg_labels is not None and sv not in fitted_seg_labels:
                return 'segim_differs_upstream',confidence+'_NOT_FITTED_REGION_LABEL_R9'
            # Only if no saved fitted-label evidence exists, preserve the
            # conservative legacy diagnostic, explicitly labelled inferred.
            if fitted_seg_labels is None:
                orig_lab=pass_label_mapping.get(sv) if pass_label_mapping is not None else sv
                if orig_lab is None or orig_lab not in labels:
                    return 'segim_differs_upstream',confidence+'_FALLBACK_LABEL_MAPPING_INFERRED'
        if im is not None and im.shape==(H,W) and not np.isfinite(im[py,px]):
            if original_image is not None and not np.isfinite(original_image[y,x]):
                return 'source_invalid','VERIFIED_ORIGINAL_SOURCE_NONFINITE'
            return 'autoclip_NA',confidence+'_NONFINITE_PRE_TIGHTCROP_IMAGE_AUTCLIP_UNCONFIRMED'
        if original_image is not None and not np.isfinite(original_image[y,x]):
            return 'source_invalid','VERIFIED_ORIGINAL_SOURCE_NONFINITE'
        if msk is not None and msk.shape==(H,W) and np.isfinite(msk[py,px]) and msk[py,px]!=0:
            return 'detector_mask',confidence+'_PASSED_MASK_NONZERO_DETECTOR_ORIGIN_UNPROVEN'
        if sig is not None and sig.shape==(H,W) and (not np.isfinite(sig[py,px]) or sig[py,px]<=0):
            return 'invalid_sigma',confidence+'_INVALID_PASSED_SIGMA'
    else:
        if original_image is not None and not np.isfinite(original_image[y,x]):
            return 'source_invalid','VERIFIED_ORIGINAL_SOURCE_NONFINITE_NO_INPUT_DOMAIN'
    inside=y0<=y<y0+h and x0<=x<x0+w
    if inside:
        i,j=y-y0,x-x0
        data_seg=data_arrays.get((fit_source,'segim'))
        if data_seg is not None and fitted_seg_labels is not None:
            ds=orient(data_seg,orient_name)
            if ds.shape==(h,w) and np.isfinite(ds[i,j]) and int(ds[i,j]) not in fitted_seg_labels:
                return 'segim_differs_upstream','OBSERVED_DATA_SEGIM_NOT_FITTED_LABEL_R9'
        dimage=data_arrays.get((fit_source,'image'))
        if dimage is not None:
            z=orient(dimage,orient_name)
            if z.shape==(h,w) and not np.isfinite(z[i,j]):
                return 'autoclip_NA','OBSERVED_SAVED_DATA_IMAGE_NONFINITE_CAUSE_UNVERIFIED'
        mask=data_arrays.get((fit_source,'mask'))
        if mask is not None:
            z=orient(mask,orient_name)
            if z.shape==(h,w) and np.isfinite(z[i,j]) and z[i,j]!=0:
                return 'detector_mask','OBSERVED_SAVED_DATA_MASK_NONZERO_DETECTOR_ORIGIN_UNPROVEN'
        sigma=data_arrays.get((fit_source,'sigma'))
        if sigma is not None:
            z=orient(sigma,orient_name)
            if z.shape==(h,w) and (not np.isfinite(z[i,j]) or z[i,j]<=0):
                return 'invalid_sigma','OBSERVED_SAVED_DATA_SIGMA_INVALID'
    return 'unexplained', ('OUTSIDE_TIGHTCROP_STAMP_NOT_EVIDENCE_OF_CUTBOX' if not inside
                           else 'INSIDE_SAVED_STAMP_NOT_IN_REGION_UNRESOLVED')


def unique_seg_offset(full_seg, sub_seg):
    """Find a unique EXACT integer-valued segmentation window in a parent map.

    Report ambiguity instead of choosing the highest-scoring location.
    """
    if full_seg is None or sub_seg is None:return None
    if not np.isfinite(full_seg).all() or not np.isfinite(sub_seg).all():return None
    H,W=full_seg.shape
    if not np.all(np.equal(full_seg,np.rint(full_seg))) or not np.all(np.equal(sub_seg,np.rint(sub_seg))):return None
    matches=[]
    for ori in ('direct','transpose'):
        t=orient(sub_seg,ori)
        h,w=t.shape
        if h>H or w>W:continue
        if h==H and w==W:
            if np.array_equal(full_seg,t):matches.append((0,0,ori))
            continue
        # TM_SQDIFF provides possible candidates; test equality on the original
        # double-valued labels before accepting (float32 can round IDs).
        try:
            import cv2
            scores=cv2.matchTemplate(np.asarray(full_seg,dtype=np.float32),np.asarray(t,dtype=np.float32),cv2.TM_SQDIFF)
            score_min=float(scores.min())
            eps=max(1e-3,abs(score_min)*1e-5)
            possible=np.argwhere(scores<=score_min+eps)
            if len(possible)>3000:return None
            for r,c in possible:
                if np.array_equal(full_seg[r:r+h,c:c+w],t):matches.append((int(r),int(c),ori))
        except ImportError:
            # Pure NumPy FFT fallback using exact-check of the candidate minima.
            def fftvalid(x,k):
                dim=(x.shape[0]+k.shape[0]-1,x.shape[1]+k.shape[1]-1)
                z=np.fft.irfftn(np.fft.rfftn(x,s=dim)*np.fft.rfftn(k,s=dim),s=dim)
                return z[k.shape[0]-1:x.shape[0],k.shape[1]-1:x.shape[1]]
            f=np.asarray(full_seg,dtype=float);q=np.asarray(t,dtype=float)
            sq=fftvalid(f*f,np.ones(q.shape,float))-2*fftvalid(f,q[::-1,::-1])+np.sum(q*q)
            idx=np.argwhere(sq<=max(1e-6,1e-9*np.sum(q*q)))
            if len(idx)>1000:return None
            for r,c in idx:
                if np.array_equal(full_seg[r:r+h,c:c+w],t):matches.append((int(r),int(c),ori))
    # A mathematically symmetric stamp can match in two orientations at the
    # same coordinates. In such cases image registration must disambiguate.
    return matches[0] if len(set(matches))==1 else None


def compare_family(L,labels,reg,method,arrays,fit_source,inputs_seg,inputs_reg,ledger_writer,gal_id,
                   original_image=None, passed_image=None, passed_seg=None, passed_reg=None,
                   passed_mask=None, passed_sigma=None):
    if method is None or method.get('status')!='VERIFIED_PIXEL_REGISTRATION':return dict(status='OFFSET_UNRESOLVED')
    _,y0,x0,ori,validation=method['best']; rr=orient(reg,ori)
    h,w=rr.shape
    if y0<0 or x0<0 or y0+h>L.shape[0] or x0+w>L.shape[1]:return dict(status='STAMP_OUTSIDE_ORIGINAL_GRID')
    target=np.isin(L,labels)
    region=np.isfinite(rr)&(rr!=0)
    iy,ix=np.nonzero(target)
    in_box=(iy>=y0)&(iy<y0+h)&(ix>=x0)&(ix<x0+w)
    missing_y=[];missing_x=[]
    for y,x,inside in zip(iy,ix,in_box):
        if not inside or not region[y-y0,x-x0]:
            missing_y.append(int(y));missing_x.append(int(x))
    reg_y,reg_x=np.nonzero(region)
    extra_y=reg_y+y0;extra_x=reg_x+x0
    extra_keep=~target[extra_y,extra_x]
    extra_original_labels=L[extra_y[extra_keep],extra_x[extra_keep]]
    n_extra_background=int(np.count_nonzero(extra_original_labels==0))
    n_extra_neighbour=int(np.count_nonzero(extra_original_labels!=0))
    support_coords=np.stack((reg_y+y0,reg_x+x0),axis=1).astype('<i4',copy=False)
    support_hash=hashlib.sha256(support_coords.tobytes()).hexdigest()
    local_L=L[y0:y0+h,x0:x0+w]
    label_match=bool(np.array_equal(region,np.isin(local_L,labels)))
    individual_match=any(np.array_equal(region,local_L==label) for label in np.unique(local_L) if label>0)
    covered_labels=np.unique(local_L[region]);covered_labels=covered_labels[covered_labels>0]
    union_match=bool(covered_labels.size and np.all(local_L[region]>0) and np.array_equal(region,np.isin(local_L,covered_labels)))
    dseg=arrays.get((fit_source,'segim'))
    data_partition=dict(status='NOT_SAVED',mismatch_pixels=None)
    if dseg is not None:
        ds=orient(dseg,ori)
        data_partition=partition_comparison(local_L,ds,ignored=None)
        mask=arrays.get((fit_source,'mask'))
        img=arrays.get((fit_source,'image'))
        if data_partition['status']=='PARTITION_DIFFERS' and mask is not None and img is not None:
            m=orient(mask,ori);q=orient(img,ori)
            if m.shape==ds.shape and q.shape==ds.shape:
                restriction=partition_comparison(local_L,ds,ignored=(~np.isfinite(q))|(np.isfinite(m)&(m!=0)))
                if restriction['status'] in ('PARTITION_IDENTICAL','PARTITION_IDENTICAL_RELABELED'):
                    data_partition['status']='PARTITION_DIFFERS_ONLY_NA_MASK'
        # Saved Data$segim is evidence, not itself proof of the Found2Fit branch.
    pass_partition=dict(status='NOT_SAVED',mismatch_pixels=None)
    pass_mapping=None
    if passed_seg is not None and passed_reg is not None:
        a,b,H,W=passed_reg[:4]
        if (a>=0 and b>=0 and a+H<=L.shape[0] and b+W<=L.shape[1] and passed_seg.shape==(H,W)):
            pass_partition=partition_comparison(L[a:a+H,b:b+W],passed_seg)
            pass_mapping=pass_partition.get('map_new_to_orig')
    elif inputs_seg is not None and inputs_reg is not None:
        a,b,H,W=inputs_reg[:4]
        if a>=0 and b>=0 and a+H<=L.shape[0] and b+W<=L.shape[1]:
            pass_partition=partition_comparison(L[a:a+H,b:b+W],inputs_seg)
            pass_mapping=pass_partition.get('map_new_to_orig')
    # Data$segim is the closest saved evidence of the labels actually fitted.
    # Recover labels only beneath Data$region, ignoring zero/background.
    # If Data$segim is absent, use passed ProFound segim inside its registered
    # domain. Never guess the fitted label from a WAVES majority projection.
    fitted_seg_labels=None
    if dseg is not None:
        ds=orient(dseg,ori)
        if ds.shape==region.shape:
            vals=ds[region & np.isfinite(ds)]
            if vals.size and np.all(vals==np.rint(vals)):
                fitted_seg_labels=set(map(int,vals[vals>0]))
                if not fitted_seg_labels:fitted_seg_labels=None
    if fitted_seg_labels is None and passed_seg is not None and passed_reg is not None:
        a,b,H,W=passed_reg[:4]
        if passed_seg.shape==(H,W) and a<=y0 and b<=x0 and y0+h<=a+H and x0+w<=b+W:
            px=passed_seg[y0-a:y0-a+h,x0-b:x0-b+w]
            vals=px[region & np.isfinite(px)]
            if vals.size and np.all(vals==np.rint(vals)):
                fitted_seg_labels=set(map(int,vals[vals>0])) or None
    class_counts=Counter();reason_counts=Counter()
    for y,x in zip(missing_y,missing_x):
        cls,reason=classify_pixel(y,x,(y0,x0,h,w),inputs_seg,inputs_reg,arrays,fit_source,ori,labels,
                original_image=original_image,passed_image=passed_image,passed_seg=passed_seg,
                passed_registration=passed_reg,passed_mask=passed_mask,passed_sigma=passed_sigma,
                pass_label_mapping=pass_mapping,fitted_seg_labels=fitted_seg_labels)
        class_counts[cls]+=1;reason_counts[reason]+=1
        ledger_writer.writerow([gal_id,fit_source,'missing',y,x,cls,reason])
    for y,x,original_label in zip(extra_y[extra_keep],extra_x[extra_keep],extra_original_labels):
        # This measures original-label provenance, NOT whether an expansion
        # function was executed. A background extra is consistent with dilation
        # but can also arise from other fitted-region choices.
        if original_label==0:
            extra_class='extra_original_background'
            evidence='VERIFIED_ORIGINAL_WAVES_L_ZERO;DILATION_CONSISTENT_NOT_CAUSAL_PROOF'
        else:
            extra_class='extra_original_neighbour'
            evidence='VERIFIED_ORIGINAL_WAVES_LABEL_'+str(int(original_label))+';ENCROACHMENT_OBSERVED_CAUSE_UNVERIFIED'
        ledger_writer.writerow([gal_id,fit_source,'extra',int(y),int(x),extra_class,evidence])
    branch='CONSISTENT_WITH_STATIC_PARTITION_NOT_PROOF' if pass_partition['status'] in ('PARTITION_IDENTICAL','PARTITION_IDENTICAL_RELABELED') else 'UNDETERMINED_SAVED_PARTITION_NOT_EXECUTED_BRANCH_PROOF'
    calc=arrays.get((fit_source,'calcregion'))
    calc_n=calc_only=calc_outside=None
    if calc is not None:
        calc=orient(calc,ori);calc_region=np.isfinite(calc)&(calc!=0)
        calc_n=int(calc_region.sum())
        if calc.shape==region.shape:calc_only=int(np.count_nonzero(calc_region&~region));calc_outside=0
        else:calc_outside='DIFFERENT_DIMENSIONS_OR_COORDINATE_DOMAIN'
    mask=arrays.get((fit_source,'mask'));sig=arrays.get((fit_source,'sigma'))
    tlocal=np.isin(local_L,labels)
    mask_n=invalid_n=None
    if mask is not None and orient(mask,ori).shape==rr.shape:
        m=orient(mask,ori);mask_n=int(np.count_nonzero(tlocal & np.isfinite(m) & (m!=0)))
    if sig is not None and orient(sig,ori).shape==rr.shape:
        s=orient(sig,ori);invalid_n=int(np.count_nonzero(tlocal & (~np.isfinite(s)|(s<=0))))
    source_dmax=None;is_exact=False
    if original_image is not None:
        a=np.asarray(original_image[y0:y0+h,x0:x0+w],float)
        z=arrays.get((fit_source,'image'))
        if z is not None:
            b=orient(z,ori)
            valid=np.isfinite(a)&np.isfinite(b)
            if valid.any():
                source_dmax=float(np.max(np.abs(a[valid]-b[valid])))
                is_exact=bool(np.array_equal(np.isfinite(a),np.isfinite(b)) and source_dmax==0.)
    return dict(status='COMPARE_COMPLETE_WITH_INFERRED_CAUSES',N_target=int(target.sum()),N_region=int(region.sum()),
      region_original_coords_hash=support_hash,N_missing=len(missing_x),N_extra=int(extra_keep.sum()),N_extra_original_background=n_extra_background,N_extra_original_neighbour=n_extra_neighbour,class_counts=dict(class_counts),reason_counts=dict(reason_counts),
      offset_row0=y0,offset_col0=x0,orientation=ori,offset_method=('EXACT_COPY_STAMP' if is_exact else 'UNIQUE_IMAGE_PIXEL_CORRELATION'),
      offset_validation=validation,offset_uniqueness=method['n_alternatives'],
      stamp_exact_copy=is_exact,stamp_max_abs_original_delta=source_dmax,
      stamp_scale=validation.get('scale'),stamp_shift=validation.get('sky_offset'),
      data_seg_partition=data_partition['status'],profound_seg_partition=pass_partition['status'],
      data_seg_partition_mismatch_pixels=data_partition.get('mismatch_pixels'),
      profound_seg_partition_mismatch_pixels=pass_partition.get('mismatch_pixels'),
      matches_target_label=label_match,matches_any_individual_label=individual_match,matches_union_original_labels=union_match,found2fit_branch=branch,
      N_calcregion=calc_n,N_calc_only=calc_only,N_calc_outside_stamp=calc_outside,
      N_masked_target_in_stamp=mask_n,N_invalid_sigma_target_in_stamp=invalid_n,
      missing_xy=(missing_y,missing_x),extra_xy=(extra_y[extra_keep].tolist(),extra_x[extra_keep].tolist()),box=(y0,x0,h,w))


def draw_overlay(out,path,gid,image,L,labels,comparison,source):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    target=np.isin(L,labels)
    if not target.any():return
    ys,xs=np.nonzero(target)
    bounds=[int(ys.min()),int(ys.max()),int(xs.min()),int(xs.max())]
    y0,x0,h,w=comparison['box'];bounds=[min(bounds[0],y0),max(bounds[1],y0+h-1),min(bounds[2],x0),max(bounds[3],x0+w-1)]
    pad=35;u=max(0,bounds[0]-pad);v=min(L.shape[0],bounds[1]+pad+1);s=max(0,bounds[2]-pad);t=min(L.shape[1],bounds[3]+pad+1)
    z=np.asarray(image[u:v,s:t],dtype=float)
    finite=z[np.isfinite(z)];lo,hi=(np.percentile(finite,[2,99]) if finite.size else (0,1))
    fig,ax=plt.subplots(figsize=(8,7));ax.imshow(z,cmap='gray',origin='lower',vmin=lo,vmax=hi,extent=(s,t,u,v))
    yy,xx=np.mgrid[u:v,s:t]
    ax.contour(xx+.5,yy+.5,target[u:v,s:t].astype(float),levels=[.5],colors=['cyan'],linewidths=[1.2])
    my,mx=comparison['missing_xy'];ey,ex=comparison['extra_xy']
    if my:ax.scatter(np.array(mx)+.5,np.array(my)+.5,s=3,c='red',alpha=.75,label='missing')
    if ey:ax.scatter(np.array(ex)+.5,np.array(ey)+.5,s=3,c='yellow',alpha=.75,label='extra')
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((x0,y0),w,h,fill=False,ec='lime',lw=1.6,ls='--',label='saved stamp'))
    ax.set_xlim(s,t);ax.set_ylim(u,v);ax.set_aspect('equal');ax.set_xlabel('Original FITS x (0-based)');ax.set_ylabel('Original FITS y (0-based)')
    ax.set_title(f'W{gid} {source}: missing={comparison["N_missing"]} extra={comparison["N_extra"]}')
    ax.legend(loc='upper right',fontsize=8);fig.tight_layout();fig.savefig(path,dpi=140);plt.close(fig)


def create_runtime(rscript,out,historical_rlib):
    code=("cat('R_VERSION=',R.version.string,'\\n',sep='');"
      "cat('LIBS=',paste(.libPaths(),collapse=';'),'\\n',sep='');"
      "pin=Sys.getenv('BD_PINNED_RLIB');"
      "cat('PINNED_LIB=',pin,' EXISTS=',dir.exists(pin),'\\n',sep='');"
      "for(p in c('ProFit','ProFuse','ProFound','Highlander','LaplacesDemon','Rfits')){"
      "for(scope in c('CURRENT','PINNED')){"
      "loc=if(scope=='PINNED') pin else .libPaths();"
      "ds=if(scope=='PINNED' && !dir.exists(pin)) character() else "
      "tryCatch(find.package(p,lib.loc=loc,quiet=TRUE),error=function(e)character());"
      "cat('PACKAGE=',p,' SCOPE=',scope,' PATH=',paste(ds,collapse=';'),'\\n',sep='');"
      "if(length(ds)){d=file.path(ds[1],'DESCRIPTION');if(file.exists(d)){"
      "a=tryCatch(read.dcf(d),error=function(e)NULL);if(!is.null(a)){"
      "for(k in intersect(c('Version','RemoteSha','GithubSHA1','Built'),colnames(a)))"
      "cat('FIELD=',p,'/',scope,'/',k,'/',a[1,k],'\\n',sep='');"
      "cat('DESCRIPTION_MD5=',p,'/',scope,'/',unname(tools::md5sum(d)),'\\n',sep='')}};"
      "dll=list.files(file.path(ds[1],'libs'),pattern='[.](so|dylib)$',recursive=TRUE,full.names=TRUE);"
      "if(length(dll))for(f in dll)cat('DLL_MD5=',p,'/',scope,'/',basename(f),'/',unname(tools::md5sum(f)),'\\n',sep='')"
      "};cat('RUNTIME_PROBE=COMPLETE\\n',sep='');}}}")
    env={**os.environ,'BD_PINNED_RLIB':historical_rlib}
    proc=subprocess.run([rscript,'--vanilla','-e',code],capture_output=True,text=True,env=env)
    report={'rscript':rscript,'returncode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr,
            'pinned_historical_rlib':historical_rlib,'python':sys.version,'numpy':np.__version__}
    failures=[]
    if proc.returncode != 0:
        failures.append('RUNTIME_PROBE_NONZERO_EXIT')
    if 'RUNTIME_PROBE=COMPLETE' not in proc.stdout:
        failures.append('RUNTIME_PROBE_SENTINEL_MISSING')
    report['runtime_gate']='PASS' if not failures else 'FAIL'
    report['runtime_gate_failures']=failures
    report['runtime_id']=hashlib.sha256(json.dumps(report,sort_keys=True).encode()).hexdigest()
    with open(out/'runtime.json','w') as f:json.dump(report,f,indent=2)
    if failures:
        raise RuntimeError('Runtime provenance probe failed: '+', '.join(failures))
    return report['runtime_id']


def main():
    ap=argparse.ArgumentParser(description='Read-only Hector/WAVES Phase A exporter; no fitting')
    ap.add_argument('--out',required=True,help='New, nonexistent output directory (never historical /data or /scratch)')
    ap.add_argument('--historical-rlib',default='/data/uqpdas/WAVES_v38_profit272_pilot/Rlib',help='Read-only pinned historical library inspected explicitly via lib.loc')
    ap.add_argument('--root',default=KNOWN_ROOT,help='Historical 842 parent root, for RDS evidence only')
    ap.add_argument('--sample-ids',default=KNOWN_R61,help='Frozen 832-row ID TSV (roster only; no selection mixing)')
    ap.add_argument('--labels',help='Verified catalogue-to-WAVES-label TSV; columns uberID and labels')
    ap.add_argument('--path-map',help='Optional exact, human-verified recorded_path to verified_actual_path TSV')
    ap.add_argument('--paths-tsv',help='Stage0 collector paths.tsv, recorded evidence only, NOT remapping authority')
    ap.add_argument('--stage0-log',help='Stage0 evidence log for provenance (never executed as source)')
    ap.add_argument('--fit-binding',help='Optional exact per-galaxy fit binding: uberID, selected_fit_path, selection_layer')
    ap.add_argument('--seg-hdu',help='Explicit upstream confirmed segmap HDU name/index (no default)')
    ap.add_argument('--image-hdu',help='Explicit upstream confirmed image HDU name/index (no default)')
    ap.add_argument('--seg-product',choices=['detection','dilated_total_photometry','regrouped'],help='Upstream confirmed WAVES segmentation product identity')
    ap.add_argument('--rscript',required=True)
    subset=ap.add_mutually_exclusive_group()
    subset.add_argument('--dry-run',action='store_true',help='Process first 3 frozen IDs only; explicitly marked not 832 complete')
    subset.add_argument('--pilot-13',action='store_true',help='Process first 3 plus available 10 frozen pilot IDs; never expand beyond 832 roster')
    ap.add_argument('--keep-intermediates',action='store_true')
    args=ap.parse_args()
    if not os.path.isdir(args.root):
        print('WARNING: historical root not present; rows will state missing files',file=sys.stderr)
    out=Path(args.out).expanduser().absolute()
    # Fail closed: user-owned new output path, no /data or /scratch or existing path.
    if out.exists():ap.error('output directory already exists; pick a NEW name')
    parent_resolved=out.parent.resolve()
    if parent_resolved==Path('/data') or parent_resolved==Path('/scratch') or str(parent_resolved).startswith(('/data/','/scratch/')):
        ap.error('output parent resolves to /data or /scratch: forbidden')
    if str(out).startswith(('/data/','/scratch/')) or str(out) in ('/data','/scratch'):
        ap.error('output under /data or /scratch forbidden')
    for forbidden in (Path(args.root).absolute(),Path('/data'),Path('/scratch')):
        if out==forbidden or forbidden in out.parents:ap.error('output inside a historical data tree forbidden')
    ids=sample_ids(args.sample_ids)
    label_map=mappings(args.labels);remap=path_map(args.path_map)
    bindings={}
    if args.fit_binding:
        for r in read_tsv(args.fit_binding):
            k=id_field(r)
            if k in bindings:raise ValueError('duplicate fit-binding '+k)
            v=r.get('selected_fit_path');layer=r.get('selection_layer','')
            if not v or not os.path.isabs(v):raise ValueError('selected_fit_path must be absolute for '+k)
            if layer not in ('R53','R54','R61'):raise ValueError('binding selection_layer must be R53/R54/R61 (not R74)')
            bindings[k]=(v,layer)
    if args.pilot_13:
        frozen=set(ids)
        pilot_missing=[x for x in PILOTS if x not in frozen]
        if pilot_missing:raise ValueError('FROZEN_PILOT_IDS_NOT_IN_832_ROSTER: '+','.join(pilot_missing))
        ids=list(dict.fromkeys(ids[:3]+PILOTS))
    elif args.dry_run:ids=ids[:3]
    out.mkdir(parents=True,exist_ok=False)
    (out/'overlays').mkdir();(out/'scratch').mkdir();(out/'tmp').mkdir();(out/'mplconfig').mkdir()
    os.environ['TMPDIR']=str(out/'tmp')
    os.environ['MPLCONFIGDIR']=str(out/'mplconfig')
    os.environ['XDG_CACHE_HOME']=str(out/'tmp')
    os.environ['PYTHONDONTWRITEBYTECODE']='1'
    parse_test=subprocess.run([args.rscript,'--vanilla','-e',
       'invisible(parse(file=commandArgs(trailingOnly=TRUE)[1L]));cat("R_SYNTAX_GATE=PASS\\n")',
       str(Path(__file__).with_name('export_rds_arrays.R'))],capture_output=True,text=True)
    with open(out/'r_parse_gate.txt','w') as f:f.write(parse_test.stdout+'\nSTDERR:\n'+parse_test.stderr)
    if parse_test.returncode!=0 or 'R_SYNTAX_GATE=PASS' not in parse_test.stdout:
        raise RuntimeError('R syntax parse gate failed before any RDS read: '+parse_test.stderr[-400:])
    print('R_SYNTAX_GATE=PASS',flush=True)
    stage0={}
    if args.paths_tsv:
        stage0['paths_tsv_sha256']=sha256_file(args.paths_tsv)
        stage0['paths_tsv_rows']=len(read_tsv(args.paths_tsv))
    if args.stage0_log:stage0['stage0_log_sha256']=sha256_file(args.stage0_log)
    runtime_id=create_runtime(args.rscript,out,args.historical_rlib)
    config={'schema':'BD_D1_PHASE_A_V04','sample_n':len(ids),'sample_ids_sha256':sha256_file(args.sample_ids),
       'sample_ids_selection_use':'ROSTER_ONLY_NO_ENDPOINT_SELECTION','root':args.root,
       'dry_run':args.dry_run,'pilot_13':args.pilot_13,'historical_rlib':args.historical_rlib,'runtime_scope':'INSTALLED_FILES_TODAY_NOT_PROOF_OF_HISTORICAL_EXECUTION','labels_file_sha256':sha256_file(args.labels) if args.labels else '',
       'path_map_sha256':sha256_file(args.path_map) if args.path_map else '',
       'fit_binding_sha256':sha256_file(args.fit_binding) if args.fit_binding else '',
       'seg_hdu':args.seg_hdu,'seg_product':args.seg_product or 'UNRESOLVED','image_hdu':args.image_hdu,'stage0':stage0,
       'evidence_levels':['VERIFIED_BYTES','REPORTED','INFERRED','UNRESOLVED'],
       'gate_G2':'NOT_TESTED_NEW_INPUTS_NOT_BUILT','runtime_id':runtime_id}
    with open(out/'run_config.json','w') as f:json.dump(config,f,indent=2)
    counts=Counter();extra_counts=Counter();all_family_pixel_counts=Counter();status=Counter();bad=Counter();rows=[];family_counts=Counter()
    inventory=[]
    with gzip.open(out/'pixel_ledger.tsv.gz','wt',newline='') as ledger:
        lw=csv.writer(ledger,delimiter='\t');lw.writerow(['uberID','fit_source','kind','row_y0','col_x0','class','evidence'])
        for i,g in enumerate(ids,1):
            candidates=sorted((Path(args.root)/g).glob('r_1333_v*/inputs.rds'))
            # Keep ambiguity visible. Never silently select one historical version.
            ip=candidates[0] if len(candidates)==1 else Path(args.root)/g/'UNRESOLVED_FAMILY_VERSION'/'inputs.rds'
            row={x:'' for x in FIELDS};row.update(id=g,inputs_path=str(ip),runtime_id=runtime_id,
              image_hdu=args.image_hdu or '',seg_hdu=args.seg_hdu or '',segmentation_product=args.seg_product or 'UNRESOLVED',fits_available=False,
              label_provenance='USER_VERIFIED_LABEL_MAP' if g in label_map else 'UNRESOLVED_NO_DOCUMENTED_LABEL_SET',
              labels=','.join(map(str,label_map.get(g,[]))),provenance='D1_HISTORICAL_AUDIT_NOT_G2_ACCEPTANCE')
            issues=[];fam={}
            row['input_folders_found']=';'.join(str(x.parent) for x in candidates)
            row['input_folder']=str(ip.parent) if len(candidates)==1 else ''
            if not candidates:issues.append('NO_R_1333_VERSION_FOUND')
            if len(candidates)>1:issues.append('MULTIPLE_R_1333_VERSIONS_NO_GUESS')
            try:
                if not ip.is_file():issues.append('INPUT_RDS_NOT_FOUND');continue_row=True
                else:continue_row=False
                if not continue_row:
                    extraction=out/'scratch'/g;extraction.mkdir(exist_ok=False)
                    ff=[ip.parent/x for x in ('fit_FS.rds','fit_PD.rds','fit_BD.rds')]
                    binding=bindings.get(g)
                    if binding and Path(binding[0]) not in ff:ff.append(Path(binding[0]))
                    p=subprocess.run([args.rscript,'--vanilla',str(Path(__file__).with_name('export_rds_arrays.R')),str(ip),str(extraction),'|'.join(map(str,ff))],capture_output=True,text=True)
                    if p.returncode!=0:issues.append('RDS_EXTRACTION_ERROR:'+p.stderr[-200:])
                    else:
                        meta=parse_metadata(extraction/'metadata.tsv');arrays=parse_arrays(extraction/'arrays.tsv',extraction)
                        # Always report saved fit/profound field inventory, even
                        # when labels/FITS/Step0 HDUs are not yet available.
                        row['historical_finesample_by_fit_json']=json.dumps({k:meta.get(k,{}).get('Data_finesample','NOT_EXTRACTED') for k in ('fit_FS','fit_PD','fit_BD')},sort_keys=True)
                        row['historical_psftype_by_fit_json']=json.dumps({k:meta.get(k,{}).get('Data_psftype','NOT_EXTRACTED') for k in ('fit_FS','fit_PD','fit_BD')},sort_keys=True)
                        prof_inventory={}
                        for fs in ('fit_FS','fit_PD','fit_BD'):
                            fields={k:meta.get(fs,{}).get('profound_presence_'+k,'NOT_EXTRACTED')
                                    for k in ('image','segim','segim_orig','mask')}
                            prof_inventory[fs]=fields
                            for key in ('image','segim','segim_orig','mask'):
                                arr=arrays.get((fs,'profound_'+key))
                                inventory.append(dict(id=g,fit_source=fs,key=key,presence=fields[key],
                                                      shape=('x'.join(map(str,arr.shape)) if arr is not None else ''),
                                                      rds_names=meta.get(fs,{}).get('__names__',''),
                                                      profound_names=meta.get(fs,{}).get('__profound_names__','')))
                        row['profound_fields_json']=json.dumps(prof_inventory,sort_keys=True)
                        row['source_image_path']=find_field(meta,'inputs',['image_path'])
                        row['segmap_path']=find_field(meta,'inputs',['segmap_path'])
                        row['sigma_source']=find_field(meta,'inputs',['sigma_path','sigma_source','weight_path']) or 'UNRESOLVED'
                        row['sky_source']=find_field(meta,'inputs',['sky','skyRMS','sky_rms']) or 'UNRESOLVED'
                        row['gain_source']=find_field(meta,'inputs',['gain','effective_gain']) or 'UNRESOLVED'
                        source,provenance=resolve_recorded(row['source_image_path'],remap)
                        seg,segprov=resolve_recorded(row['segmap_path'],remap)
                        if provenance!='EXACT_RECORDED_PATH' or segprov!='EXACT_RECORDED_PATH':issues.extend([provenance,segprov])
                        fits_ok=bool(source and seg and os.path.isfile(source) and os.path.isfile(seg))
                        row['fits_available']=fits_ok
                        row.update(psf_details(meta,arrays,'fit_FS'))
                        if fits_ok and args.image_hdu is not None and args.seg_hdu is not None:
                            im,iwcs=fits_array(source,args.image_hdu);L,swcs=fits_array(seg,args.seg_hdu)
                            row['original_grid']='x'.join(map(str,L.shape));row['image_wcs_hash']=iwcs;row['segmap_wcs_hash']=swcs
                            row['source_image_sha256']=sha256_file(source);row['segmap_sha256']=sha256_file(seg)
                            if L.shape!=im.shape:issues.append('FITS_IMAGE_SEG_DIFFERENT_GRIDS_NO_COMPARISON')
                            elif g not in label_map:issues.append('DOCUMENTED_LABEL_MAPPING_MISSING')
                            else:
                                labels=label_map[g]
                                if L.shape!=(1333,1333):issues.append('WARN_NONSTANDARD_ORIGINAL_GRID_DO_NOT_DROP_EDGE_CASE')
                                if np.any(~np.isfinite(L)) or not np.issubdtype(L.dtype,np.number):issues.append('INVALID_LABEL_MAP_VALUES')
                                else:
                                    target=np.isin(L,labels)
                                    row['N_target']=int(target.sum())
                                    if not target.any():issues.append('TARGET_LABELS_NOT_PRESENT')
                                    else:
                                        yy,xx=np.nonzero(target);row['source_edge']=bool(np.any((yy==0)|(yy==L.shape[0]-1)|(xx==0)|(xx==L.shape[1]-1)))
                                        # The input RDS itself may be an island-filtered or cut image.
                                        inputim=arrays.get(('inputs','image'))
                                        inreg=best_offset(im,inputim) if inputim is not None else None
                                        input_seg=arrays.get(('inputs','segim'))
                                        input_seg_reg=None
                                        if inreg and inreg.get('status')=='VERIFIED_PIXEL_REGISTRATION' and input_seg is not None:
                                            _,iy0,ix0,iori,_=inreg['best'];input_seg=orient(input_seg,iori)
                                            if input_seg.shape==(inputim.shape if iori=='direct' else inputim.T.shape):input_seg_reg=(iy0,ix0,*input_seg.shape,iori)
                                            else:input_seg=None
                                        else:input_seg=None
                                        available=[x for x in ('fit_FS','fit_PD','fit_BD') if (x,'region') in arrays]
                                        if binding:
                                            bound='fit_'+Path(binding[0]).stem if False else Path(binding[0]).stem
                                            if bound not in available and (bound,'region') in arrays:available.append(bound)
                                        row['fit_binding']=binding[1] if binding else 'NONE_SUPPORT_ONLY_NOT_SELECTED_ENDPOINT'
                                        results={}
                                        for fit in available:
                                            reg=arrays[(fit,'region')];stamp=arrays.get((fit,'image'))
                                            if stamp is None:
                                                results[fit]={'status':'MISSING_DATA_IMAGE_OFFSET_UNRESOLVED'};continue
                                            match=best_offset(im,stamp)
                                            # If original image was sky corrected, try original input image
                                            # and compose its VERIFIED offset. Never derive from dimensions alone.
                                            if (match is None or match.get('status')!='VERIFIED_PIXEL_REGISTRATION') and inreg and inreg.get('status')=='VERIFIED_PIXEL_REGISTRATION':
                                                mid=best_offset(inputim,stamp)
                                                if mid and mid.get('status')=='VERIFIED_PIXEL_REGISTRATION':
                                                    _,ar,ac,ao,av=mid['best'];_,br,bc,bo,bv=inreg['best']
                                                    if ao==bo=='direct':match={'status':'VERIFIED_PIXEL_REGISTRATION','best':(min(av['corr'],bv['corr']),br+ar,bc+ac,'direct',av),'n_alternatives':0,'composed_input_match':True}
                                                    # Transpose compositions are not silently guessed.
                                            # Prefer a saved ProFound pre-tightcrop domain.
                                            # Registration is independent of the cropped Data arrays.
                                            profraw=arrays.get((fit,'profound_image'))
                                            profsegraw=arrays.get((fit,'profound_segim'))
                                            profmaskraw=arrays.get((fit,'profound_mask'))
                                            pmatch=best_offset(im,profraw) if profraw is not None else None
                                            pdom=None;pim=None;pseg=None;pmask=None
                                            if pmatch and pmatch.get('status')=='VERIFIED_PIXEL_REGISTRATION':
                                                _,py,px,pori,_=pmatch['best']
                                                pim=orient(profraw,pori)
                                                pdom=(py,px,*pim.shape,pori)
                                                for raw,key in ((profsegraw,'seg'),(profmaskraw,'mask')):
                                                    if raw is not None:
                                                        oriented=orient(raw,pori)
                                                        if oriented.shape==pim.shape:
                                                            if key=='seg':pseg=oriented
                                                            else:pmask=oriented
                                            # An exact segmentation match may determine tightcrop
                                            # offset from the saved ProFound map, provided its
                                            # original-grid placement was independently established.
                                            if pdom is not None and pseg is not None and (fit,'segim') in arrays:
                                                local=unique_seg_offset(pseg,arrays[(fit,'segim')])
                                                if local is not None:
                                                    lr,lc,lori=local
                                                    suggested=(pdom[0]+lr,pdom[1]+lc)
                                                    if match is not None and match.get('status')=='VERIFIED_PIXEL_REGISTRATION':
                                                        if match['best'][1:3]!=suggested or match['best'][3]!=lori:
                                                            results[fit]={'status':'PROFOUND_EXACT_SEGIM_VERSUS_IMAGE_OFFSET_CONFLICT'}
                                                            issues.append('PROFOUND_EXACT_SEGIM_VERSUS_IMAGE_OFFSET_CONFLICT:'+fit)
                                                            continue
                                                    if match is None or match.get('status')!='VERIFIED_PIXEL_REGISTRATION':
                                                        val=robust_verification(im,orient(stamp,lori),*suggested)
                                                        if val is not None and val['corr']>=.995 and val['residual']<=.10 and val['scale']>0:
                                                            match={'status':'VERIFIED_PIXEL_REGISTRATION','best':(val['corr'],*suggested,lori,val),'n_alternatives':0,'exact_segim_verified':True}
                                            if profraw is not None:
                                                row['found2fit_domain_provenance']='SAVED_PROFOUND_PRE_TIGHTCROP' if pdom is not None else 'PROFOUND_IMAGE_PRESENT_UNREGISTERED'
                                            else:
                                                row['found2fit_domain_provenance']='INPUTS_RDS_PROXY_IF_REGISTERED_NOT_PROVEN_FOUND2FIT_INPUT'
                                            res=compare_family(L,labels,reg,match,arrays,fit,input_seg,input_seg_reg,lw,g,
                                                   original_image=im,passed_image=pim,passed_seg=pseg,passed_reg=pdom,passed_mask=pmask)
                                            res['found2fit_domain_provenance']=row['found2fit_domain_provenance']
                                            res['profound_fields']=prof_inventory.get(fit,{})
                                            results[fit]=res
                                            if res['status'].startswith('COMPARE_COMPLETE'):
                                                family_counts[fit]+=1
                                                for cl,n in res['class_counts'].items():all_family_pixel_counts[cl]+=n
                                                if (g in PILOTS or res['N_missing'] or res['N_extra']):
                                                    label=f'{g}_{fit}.png'
                                                    try:
                                                        draw_overlay(out,out/'overlays'/label,g,im,L,labels,res,fit)
                                                        row['overlay_status']='WRITTEN_'+fit
                                                    except Exception as overlay_error:
                                                        issues.append('OVERLAY_ERROR:'+str(overlay_error)[:180])
                                        # Never infer an R53 endpoint from the R61 family label.
                                        selected=Path(binding[0]).stem if binding else (available[0] if len(available)==1 else None)
                                        if selected not in results:
                                            if results and all(r.get('status')=='COMPARE_COMPLETE_WITH_INFERRED_CAUSES' for r in results.values()):
                                                signatures={r['region_original_coords_hash'] for r in results.values()}
                                                if len(signatures)==1:selected=next(iter(results));row['fit_binding']='UNBOUND_ALL_FAMILY_COUNTS_AGREE'
                                            if selected not in results:issues.append('SELECTED_ENDPOINT_UNBOUND_OR_DISAGREES')
                                        if selected in results:
                                            res=results[selected];row['selected_fit_source']=selected
                                            if res['status'].startswith('COMPARE_COMPLETE'):
                                                fitted_psf=arrays.get((selected,'psf'))
                                                saved_psf=arrays.get(('psf_rds','psf'))
                                                if fitted_psf is not None:
                                                    row['fit_data_psf_hash']=hashlib.sha256(np.asarray(fitted_psf,dtype='<f8').tobytes()).hexdigest()
                                                    row['psf_matches_fit_data']=(bool(np.array_equal(saved_psf,fitted_psf,equal_nan=True)) if saved_psf is not None else 'NO_SAVED_PSF_ARRAY')
                                                row['found2fit_domain_provenance']=res.get('found2fit_domain_provenance','')
                                                for k in ('N_region','region_original_coords_hash','N_missing','N_extra','N_extra_original_background','N_extra_original_neighbour','N_calcregion','N_calc_only','N_calc_outside_stamp','N_masked_target_in_stamp','N_invalid_sigma_target_in_stamp','offset_row0','offset_col0','orientation','offset_method','offset_uniqueness','matches_target_label','matches_any_individual_label','matches_union_original_labels','found2fit_branch','data_seg_partition','profound_seg_partition','data_seg_partition_mismatch_pixels','profound_seg_partition_mismatch_pixels','stamp_exact_copy','stamp_max_abs_original_delta','stamp_scale','stamp_shift'):row[k]=res.get(k,'')
                                                row['offset_validation']=json.dumps(res.get('offset_validation'),sort_keys=True)
                                                for cl,n in res['class_counts'].items():counts[cl]+=n
                                                extra_counts['extra_original_background']+=res['N_extra_original_background']
                                                extra_counts['extra_original_neighbour']+=res['N_extra_original_neighbour']
                                                for cls in ('cutbox_side_fit_crop','island_filter','segim_differs_upstream','autoclip_NA','source_invalid','detector_mask','invalid_sigma','unexplained'):
                                                    column='N_missing_cutbox_side_fit' if cls=='cutbox_side_fit_crop' else 'N_missing_'+cls
                                                    row[column]=res['class_counts'].get(cls,0)
                                            else:issues.append(res['status'])
                                        # Ensure row is JSON serializable; pixel coordinates already recorded in ledger.
                                        for fit,res in results.items():
                                            res.pop('missing_xy',None);res.pop('extra_xy',None);res.pop('box',None)
                                        fam=results
                        else:issues.append('SOURCE_FITS_NOT_AVAILABLE_OR_HDU_NOT_DECLARED')
                        if not args.keep_intermediates:shutil.rmtree(extraction)
            except Exception as e:
                issues.append('UNHANDLED_PER_GALAXY:'+type(e).__name__+':'+str(e)[:300]);bad[g]=traceback.format_exc()
            if g in PILOTS and not row['overlay_status']:
                row['overlay_status']='BLOCKED_NO_SAFE_PIXEL_REGISTRATION_OR_SEGMENT'
            row['errors']=';'.join(dict.fromkeys(issues));row['family_results_json']=json.dumps(fam,default=lambda x:x.item() if isinstance(x,np.generic) else str(x),sort_keys=True)
            row['status']='COMPARED' if row['N_region']!='' and not issues else ('COMPARED_WITH_WARNINGS' if row['N_region']!='' else 'BLOCKED_PARTIAL')
            status[row['status']]+=1;rows.append(row)
            print(f'[{i}/{len(ids)}] {g} {row["status"]} {row["N_missing"]} missing / {row["N_extra"]} extra {row["errors"][:100]}',flush=True)
    with open(out/'profound_presence_inventory.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['id','fit_source','key','presence','shape','rds_names','profound_names'],delimiter='\t');w.writeheader();w.writerows(inventory)
    with open(out/'phase_a_832_rows.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
    summ={'schema':'BD_D1_PHASE_A_V04','observed_rows':len(rows),'expected_rows':len(ids),'pilot_13':args.pilot_13,
          'status_counts':dict(status),'per_missing_pixel_class_selected_only':dict(counts),'per_extra_pixel_class_selected_only':dict(extra_counts),'per_missing_pixel_class_all_historical_families':dict(all_family_pixel_counts),'family_compared_counts':dict(family_counts),
          'unhandled_errors':bad,'stage0':stage0,'runtime_id':runtime_id,
          'warning':'Historical support audit only; G2 requires a separate T_g ∧ Q input builder and original-grid comparison.',
          'gate_G2':'NOT_EVALUATED_FOR_NEW_INPUTS','dry_run':args.dry_run}
    with open(out/'summary.json','w') as f:json.dump(summ,f,indent=2,sort_keys=True)
    # Hash all deliverables while avoiding self-reference.
    with open(out/'SHA256SUMS.txt','w') as f:
        for path in sorted(out.rglob('*')):
            if path.is_file() and path.name!='SHA256SUMS.txt':
                f.write(f'{sha256_file(path)}  {path.relative_to(out)}\n')
    print('OUTPUT_DIRECTORY='+str(out))
    print('PHASE_A_ROWS='+str(len(rows)))
    print('G2=NOT_TESTED_NOT_NEW_INPUTS')
    if len(rows)!=len(ids):raise RuntimeError('cardinality mismatch')

if __name__=='__main__':main()
