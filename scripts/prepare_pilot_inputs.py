#!/usr/bin/env python3
"""Prepare byte-verified uncompressed copies for one explicitly requested target.

Requires SLURM. Writes a NEW directory below outputs/ only. The resulting path
index is a blocked candidate: PSF, detector-mask, units/gain/sky/noise provenance
and scientific approvals are unresolved. No observed fitting is performed.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np

UID='11690092907692'
WID='W'+UID
RA=11.091131080899515
DEC=-31.1410252394602
SEGID=4211


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def source_identity(path):
    stat=Path(path).stat()
    return dict(path=str(Path(path).resolve()),sha256=sha256(path),
                size_bytes=stat.st_size,mtime_ns=stat.st_mtime_ns)


def selected(path,column,uid,delimiter):
    with Path(path).open(newline='',encoding='utf-8-sig') as stream:
        rows=[row for row in csv.DictReader(stream,delimiter=delimiter)
              if row[column].removeprefix('W')==uid]
    if len(rows)!=1:raise ValueError(f'{path}: target must have exactly one row')
    return rows[0]


def verified_decompress(source,destination):
    """Check gzip CRC, hash output and independently stream original again."""
    first=hashlib.sha256()
    with gzip.open(source,'rb') as original,Path(destination).open('xb') as copy:
        for block in iter(lambda:original.read(1024*1024),b''):
            first.update(block);copy.write(block)
    second=hashlib.sha256()
    with gzip.open(source,'rb') as original:
        for block in iter(lambda:original.read(1024*1024),b''):second.update(block)
    output=sha256(destination)
    if first.hexdigest()!=second.hexdigest() or output!=first.hexdigest():
        raise ValueError('Uncompressed copy differs from independently decompressed original bytes')
    return dict(path=str(Path(destination).resolve()),sha256=output,
                independent_gzip_uncompressed_sha256=second.hexdigest(),
                size_bytes=Path(destination).stat().st_size,gzip_crc_verified=True)


def array_evidence(data):
    ny,nx=data.shape
    points=[(1,1),(2,5),(5,2),(nx,1),(1,ny),(nx,ny),(667,667),
            (nx//3,ny//4),(nx//4,ny//3),(nx-3,ny-7)]
    samples=[]
    for x,y in points:
        value=float(data[y-1,x-1])
        samples.append(dict(fits_x_1based=x,fits_y_1based=y,
                            value=value if np.isfinite(value) else None,
                            finite=bool(np.isfinite(value))))
    return dict(astropy_shape_yx=[ny,nx],expected_Rfits_shape_xy=[nx,ny],
                numeric_sha256_f64le_R_columnmajor=hashlib.sha256(
                    np.asarray(data,dtype='<f8',order='C').tobytes(order='C')).hexdigest(),
                finite_count=int(np.isfinite(data).sum()),samples=samples,
                orientation='Rfits[x,y] equals Astropy[y-1,x-1]; R column-major bytes '
                            'equal Astropy row-major bytes after float64 little-endian conversion')


# R's jsonlite (digits=NA) writes 15 significant digits, so the R-side sample values are
# rounded copies of the Astropy doubles. The tolerance is far below any physical difference
# but above serialisation rounding (Claude review of 489e737: all Friday samples agree to
# <4e-15 relative), and far below the difference between a value and its transpose.
SAMPLE_REL_TOL=1e-12
SAMPLE_ABS_TOL=1e-15


def sample_values_agree(expected,actual):
    """True when two serialised sample values are the same double up to JSON rounding."""
    if expected is None or actual is None:return expected is None and actual is None
    return math.isclose(expected,actual,rel_tol=SAMPLE_REL_TOL,abs_tol=SAMPLE_ABS_TOL)


def check_rfits_samples(name,expected_samples,actual_samples):
    """Production comparison of Astropy (expected) and Rfits (actual) point samples.

    Raises ValueError unless both lists have the same, nonzero length and every sample
    has identical FITS coordinates, identical finiteness and numerically equal values.
    """
    if not expected_samples or len(expected_samples)!=len(actual_samples):
        raise ValueError(f'{name}: Rfits sample count differs from Astropy sample count')
    for expected,actual in zip(expected_samples,actual_samples):
        if (expected['fits_x_1based']!=actual['fits_x_1based'] or
            expected['fits_y_1based']!=actual['fits_y_1based'] or
            expected.get('finite')!=actual.get('finite') or
            not sample_values_agree(expected['value'],actual['value'])):
            raise ValueError(f'{name}: explicit off-diagonal Rfits sample differs')


def write_json(path,value):
    with Path(path).open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def run(args):
    from astropy.io import fits
    from astropy.wcs import WCS
    import astropy
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Actual FITS preparation requires SLURM')
    root,out=Path(args.repo).resolve(),Path(args.out).resolve()
    if not (root/'.hector_bd_root').is_file() or root/'outputs' not in out.parents or out.exists():
        raise ValueError('Require a NEW target directory below repository outputs/')
    audit=selected(args.audit,'uberID',UID,'\t')
    cat=selected(args.catalogue,'UBERID',UID,',')
    if audit['status']!='IDENTICAL' or not 0<=float(audit['wcs_max_sep_arcsec'])<=.01:
        raise ValueError('Target requires existing IDENTICAL/WCS<=0.01arcsec audit')
    if (audit['cat_boundary_flag']!='False' or audit['cat_agrees_with_historical']!='True'
        or audit['cat_label_is_sky']!='False' or audit['catalogue_row_found']!='True'):
        raise ValueError('Target catalogue-label/boundary evidence is not clean')
    if int(audit['hist_segID'])!=SEGID or float(cat['RA'])!=RA or float(cat['DEC'])!=DEC:
        raise ValueError('Requested target coordinate/segmentation identity changed')
    if int(audit['image_hdu'])!=1 or int(audit['seg_hdu'])!=1:
        raise ValueError('Expected explicit source HDU1 for science and original labels')
    originals={name:Path(audit[key]).resolve() for name,key in
               (('image','image_path'),('segmap','segmap_path'))}
    expected_dir=Path('/data/uqpdas/WAVES')/UID
    for name,path in originals.items():
        if path.parent!=expected_dir or not path.name.endswith('.fits.gz'):
            raise ValueError(f'{name}: unexpected original source path')
    before={name:source_identity(path) for name,path in originals.items()}
    out.mkdir(parents=True)
    copies={name:verified_decompress(path,out/path.name[:-3]) for name,path in originals.items()}
    arrays={};headers={}
    for name,copy in copies.items():
        with fits.open(copy['path'],mode='readonly',memmap=False) as hdus:
            if len(hdus)<=1 or int(hdus[1].header.get('NAXIS',0))!=2:
                raise ValueError(f'{name}: declared HDU1 is not explicitly2D')
            data=np.asarray(hdus[1].data)
            if data.shape!=(1333,1333):raise ValueError('Unexpected original grid size')
            arrays[name]=data
            headers[name]=hdus[1].header.copy()
            copy['array']=array_evidence(data)
            (out/f'{name}_headers.txt').write_text('\n\n'.join(
                f'HDU{i}\n'+h.header.tostring(sep='\n',padding=False) for i,h in enumerate(hdus))+'\n')
    labels=arrays['segmap']
    if not np.all(np.isfinite(labels)) or np.any(labels<0) or np.any(labels!=np.floor(labels)):
        raise ValueError('Original segmentation is not finite nonnegative integer-valued')
    x,y=WCS(headers['segmap']).world_to_pixel_values(RA,DEC)
    col,row=int(np.floor(x+.5)),int(np.floor(y+.5))
    if int(labels[row,col])!=SEGID:raise ValueError('Exact original WCS catalogue pixel has wrong label')
    nearby=labels[row-1:row+2,col-1:col+2]
    if nearby.shape!=(3,3) or np.any(nearby!=SEGID):raise ValueError('Actual target lies on original-label boundary')
    manifest=dict(status='BLOCKED_CANDIDATE_INPUT_COPIES_ONLY',uberID=WID,numeric_uberID=UID,
        ra_deg=RA,dec_deg=DEC,segID=SEGID,science_hdu=1,segmap_hdu=1,
        target_astropy_pixel_xy_0based=[float(x),float(y)],
        target_exact_label=SEGID,originals_before=before,copies=copies,
        evidence_sources=dict(audit=dict(path=str(Path(args.audit).resolve()),sha256=sha256(args.audit)),
            catalogue=dict(path=str(Path(args.catalogue).resolve()),sha256=sha256(args.catalogue))),
        runtime=dict(python=os.sys.version,executable=os.sys.executable,numpy=np.__version__,astropy=astropy.__version__,
                     job=os.environ['SLURM_JOB_ID']),
        source_sha256=dict(python=sha256(__file__),Rhelper=sha256(root/'scripts/verify_pilot_fits.R'),
                           pipeline_io=sha256(root/'R/io_fits.R')),
        calibration_fields={name:{key:str(header[key]) if key in header else None
            for key in ('BUNIT','MAGZERO','MAGZP','PHOTZP','GAIN','EXPTIME','DMAG','FLXSCALE')}
            for name,header in headers.items()})
    manifest_path=out/'copy_manifest.json'
    write_json(manifest_path,manifest)
    result=subprocess.run([args.rscript,'--vanilla',str(root/'scripts/verify_pilot_fits.R'),
                           str(manifest_path),str(out/'Rfits_verification.json')],
                          cwd=root,capture_output=True,text=True,timeout=1200)
    (out/'Rfits_verification.out').write_text(result.stdout)
    (out/'Rfits_verification.err').write_text(result.stderr)
    if result.returncode:raise RuntimeError('Pinned Rfits orientation/read verification failed: '+result.stderr[-1500:])
    verification=json.loads((out/'Rfits_verification.json').read_text())
    for name,copy in copies.items():
        observed=verification['products'][name]
        if (observed['numeric_sha256_f64le_R_columnmajor']!=copy['array']['numeric_sha256_f64le_R_columnmajor']
            or observed['matrix_dims']!=copy['array']['expected_Rfits_shape_xy']
            or observed['finite_count']!=copy['array']['finite_count']):
            raise ValueError(f'{name}: full Astropy/Rfits numeric orientation identity failed')
        check_rfits_samples(name,copy['array']['samples'],observed['samples'])
        if sha256(copy['path'])!=copy['sha256']:raise ValueError('R read modified copied FITS bytes')
    after={name:source_identity(path) for name,path in originals.items()}
    if before!=after:raise RuntimeError('Original compressed SHA/mtime/size changed during preparation')
    write_json(out/'originals_after.json',dict(originals_unchanged=True,originals_after=after))
    candidate=dict(uberID=WID,image_path=copies['image']['path'],segmap_path=copies['segmap']['path'],
        image_hdu=1,segmap_hdu=1,psf_path='',ra_deg=RA,dec_deg=DEC,segID=SEGID,
        psf_fwhm_px='',detector_mask_policy='UNRESOLVED',
        candidate_status='BLOCKED_PSF_MASK_CALIBRATION_AND_SCIENTIFIC_GATES',
        image_sha256=copies['image']['sha256'],segmap_sha256=copies['segmap']['sha256'])
    with (out/'candidate_path_index.tsv').open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(candidate),delimiter='\t')
        writer.writeheader();writer.writerow(candidate)
    write_json(out/'candidate_step0.json',dict(uberID=WID,status='BLOCKED',
        units_zeropoint='UNRESOLVED',effective_gain='UNRESOLVED',upstream_sky='UNRESOLVED',
        noise_protocol='UNRESOLVED',detector_mask_policy='UNRESOLVED',
        psf_path_provenance='UNRESOLVED',psf_original_E1_E7_G5='NOT_PROVIDED_NOT_RUN',
        input_copy_verification='EXACT_GZIP_BYTES_AND_ASTROPY_RFITS_ORIENTATION',
        observed_fitting_authorized=False,config_modified=False,approval_flags_modified=False))
    files=sorted(p for p in out.iterdir() if p.is_file())
    (out/'SHA256SUMS.txt').write_text(''.join(f'{sha256(p)}  {p.name}\n' for p in files))
    print('ONE_TARGET_INPUT_COPIES_PREPARED_GATES_REMAIN_BLOCKED',out,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',default='.')
    parser.add_argument('--audit',default='outputs/step1_catalogue_1427444/step1_segmap_audit.tsv')
    parser.add_argument('--catalogue',default='outputs/step1_catalogue_input_832.csv')
    parser.add_argument('--rscript',default='/home/uqpdas/.conda/envs/waves_r/bin/Rscript')
    parser.add_argument('--out',required=True)
    run(parser.parse_args())


if __name__=='__main__':main()
