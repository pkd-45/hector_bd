#!/usr/bin/env python3
"""Read-only WAVES three-object calibration *diagnostics*, not an approval gate.

Source maps are opened only for reading. Measurements and their SHA256s go under
hector_bd/outputs. A ProFound gain estimate is COARSE/PROVISIONAL (not a detector
gain certification). No Step0 configuration is edited.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.ndimage import binary_dilation

IDS = ('11690092907692', '11690105411641', '41691060304604')
R_GAIN = r'''
args <- commandArgs(trailingOnly=TRUE)
if (!requireNamespace('ProFound',quietly=TRUE)) stop('ProFound missing')
ny <- as.integer(args[4]); nx <- as.integer(args[5]); n <- ny*nx
read_d <- function(p) { con<-file(p,'rb'); on.exit(close(con));
  z<-readBin(con,'double',n=n,size=8,endian='little');
  if(length(z)!=n) stop('image data incomplete');matrix(z,ny,nx) }
read_u <- function(p) { con<-file(p,'rb'); on.exit(close(con));
  z<-readBin(con,'integer',n=n,size=1,signed=FALSE,endian='little');
  if(length(z)!=n) stop('mask data incomplete'); matrix(as.logical(z),ny,nx) }
img<-read_d(args[1]); obj<-read_u(args[2]); reject<-read_u(args[3])
g<-ProFound::profoundGainEst(image=img,objects=obj,mask=reject,
                              sky=as.numeric(args[6]),skyRMS=as.numeric(args[7]))
if(length(g)!=1 || !is.finite(g) || g<=0) stop('invalid gain estimate')
cat(sprintf('COARSE_GAIN\t%.15g\n',g))
'''


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def robust(x):
    med = float(np.median(x))
    sig = float(1.4826 * np.median(np.abs(x - med)))
    return med, sig


def compute_sky(image, labels, radius=10):
    if image.ndim != 2 or labels.shape != image.shape:
        raise ValueError('image/segmap 2D grid mismatch')
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError('segmap labels must be integers')
    if np.any(labels < 0):
        raise ValueError('negative segmentation labels')
    contaminated = binary_dilation(labels != 0, iterations=radius)
    blank = (labels == 0) & ~contaminated & np.isfinite(image)
    v = np.asarray(image[blank], float)
    if len(v) < 1000:
        raise ValueError('insufficient uncontaminated sky pixels')
    med, sig = robust(v)
    if not np.isfinite(sig) or sig <= 0:
        raise ValueError('nonpositive/invalid robust sky RMS')
    z = (v - med) / sig
    std = float(np.std(z))
    skew = float(np.mean(((z - float(np.mean(z))) / std) ** 3)) if std > 0 else None
    props = {'n_sky': int(len(v)), 'pedestal': med, 'sigma_mad': sig,
             'sigma_p50_p159': float(np.percentile(v, 50) - np.percentile(v, 15.8655)),
             'z_mean': float(np.mean(z)), 'z_std': std,
             'z_skew': skew, 'z_abs_gt3_fraction': float(np.mean(np.abs(z) > 3))}
    block_medians = []
    ny, nx = image.shape
    for by in range(2):
        for bx in range(4):
            y0, y1 = by * ny // 2, (by + 1) * ny // 2
            x0, x1 = bx * nx // 4, (bx + 1) * nx // 4
            vv = image[y0:y1, x0:x1][blank[y0:y1, x0:x1]]
            block_medians.append(None if vv.size < 100 else float(np.median(vv)))
    props['sky_block_medians_2x4'] = block_medians
    return blank, props


def sky_acf(image, blank, pedestal, lag=5):
    # Pairs must be jointly in the blank set; no FFT mask-bleed artefact.
    x = np.asarray(image, dtype=float) - pedestal
    a = x[blank]
    v = float(np.mean((a - np.mean(a)) ** 2))
    if not (v > 0 and np.isfinite(v)):
        raise ValueError('zero sky variance')
    results = []
    ny, nx = x.shape
    for dy in range(-lag, lag + 1):
        for dx in range(-lag, lag + 1):
            y1 = slice(max(0, dy), min(ny, ny + dy))
            y2 = slice(max(0, -dy), min(ny, ny - dy))
            x1 = slice(max(0, dx), min(nx, nx + dx))
            x2 = slice(max(0, -dx), min(nx, nx - dx))
            valid = blank[y1, x1] & blank[y2, x2]
            n = int(valid.sum())
            rho = None
            if n >= 100:
                av = x[y1, x1][valid]; bv = x[y2, x2][valid]
                rho = float(np.mean((av - np.mean(av)) * (bv - np.mean(bv))) / v)
            results.append({'dx': dx, 'dy': dy, 'pairs': n, 'rho': rho})
    return results


def estimated_gain(image, labels, blank, pedestal, sigma, rscript, folder):
    folder = Path(folder)
    with tempfile.TemporaryDirectory(prefix='gain_scratch_', dir=folder) as tmp:
        tmp = Path(tmp)
        image_path, objects_path, reject_path = tmp/'image.bin', tmp/'objects.bin', tmp/'reject.bin'
        np.asarray(image, dtype='<f8', order='F').ravel(order='F').tofile(image_path)
        np.asarray(labels != 0, dtype='u1', order='F').ravel(order='F').tofile(objects_path)
        np.asarray(~blank, dtype='u1', order='F').ravel(order='F').tofile(reject_path)
        command = [rscript, '--vanilla', '-e', R_GAIN,
                   str(image_path), str(objects_path), str(reject_path),
                   str(image.shape[0]), str(image.shape[1]), str(pedestal), str(sigma)]
        try:
            r = subprocess.run(command, capture_output=True, text=True, timeout=900)
        except (OSError, subprocess.TimeoutExpired) as exc:
            return None, f'{type(exc).__name__}: {exc}'
        if r.returncode:
            return None, (r.stderr or r.stdout)[-1200:].replace('\n', ' | ')
        lines = [s.split('\t', 1)[1].strip() for s in r.stdout.splitlines()
                 if s.startswith('COARSE_GAIN\t')]
        if len(lines) != 1:
            return None, 'ProFound did not print exactly one gain'
        try:
            val = float(lines[0])
        except ValueError:
            return None, 'ProFound returned nonnumeric gain'
        return (val, '') if math.isfinite(val) and val > 0 else (None, 'nonfinite gain')


def selected_row(path, gid):
    with open(path, newline='', encoding='utf-8-sig') as f:
        rows = [r for r in csv.DictReader(f, delimiter='\t') if r.get('uberID') == gid]
    if len(rows) != 1:
        raise ValueError(f'{gid}: audit TSV must have exactly one row')
    return rows[0]


def run(args):
    from astropy.io import fits
    out = Path(args.out).expanduser().resolve()
    root = Path(args.repo).expanduser().resolve()
    if not (root / '.hector_bd_root').is_file():
        raise ValueError('Not the HectorBD repository')
    if root / 'outputs' not in out.parents or out.exists():
        raise ValueError('Output must be a NEW directory below repository outputs/')
    out.mkdir(parents=True)
    data = []
    details = []
    for gid in args.ids.split(','):
        gid = gid.strip().lstrip('W')
        row = selected_row(args.audit, gid)
        if row['status'] not in ('IDENTICAL', 'TARGET_REMOVED_ONLY'):
            raise ValueError(f'{gid}: unsupported historical segmentation status')
        wcs = float(row['wcs_max_sep_arcsec'])
        if wcs > 0.05:
            raise ValueError(f'{gid}: WCS discrepancy {wcs} arcsec; not safe for this 3-object control')
        ip, sp = Path(row['image_path']), Path(row['segmap_path'])
        with fits.open(ip, memmap=False) as f:
            image = np.asarray(f[args.image_hdu].data, float)
            header = f[args.image_hdu].header
        with fits.open(sp, memmap=False) as f:
            labels = np.asarray(f[args.seg_hdu].data)
        blank, sky = compute_sky(image, labels, args.dilation)
        acf = sky_acf(image, blank, sky['pedestal'], args.acf_lag)
        label = int(row['hist_segID'])
        target = (labels == label) & np.isfinite(image)
        if not np.any(target):
            raise ValueError(f'{gid}: historical label absent from delivered segmentation')
        flux = float(np.sum(image[target] - sky['pedestal']))
        m = None if flux <= 0 else 23.9 - 2.5 * math.log10(flux)
        if args.skip_gain:
            gain, err = None, 'SKIPPED by --skip-gain'
        else:
            gain, err = estimated_gain(image, labels, blank, sky['pedestal'], sky['sigma_mad'], args.rscript, out)
        bunit = str(header.get('BUNIT', '')).strip()
        item = {'uberID': gid, 'target_label': label, 'wcs_header_diff_arcsec': wcs,
                'image_sha256': sha(ip), 'segmap_sha256': sha(sp), 'image_bunit': bunit,
                'n_target': int(target.sum()), 'target_flux_native_units': flux,
                'm_at_assumed_zp23p9': m, 'gain_coarse_lower_estimate': gain,
                'gain_heuristic_error': err, **{k: v for k, v in sky.items() if not isinstance(v, list)}}
        data.append(item)
        details.append({'uberID': gid, 'sky': sky, 'acf': acf,
                        'interpretation': 'BUNIT and image zero point NOT independently certified; '
                          'profoundGainEst is a rough gain estimator and must NOT set the Step0 gain status. '
                          'ACF lags are descriptive, not an automatically valid BIC Neff.'})
        print(gid, 'sky=', sky['sigma_mad'], 'coarse_gain=', gain, 'bunit=', repr(bunit), flush=True)
    with (out / 'calibration.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(data)
    with (out / 'calibration.json').open('w') as f:
        json.dump({'status': 'DIAGNOSTIC_ONLY_NOT_GAIN_OR_UNITS_VERIFICATION',
                   'n_galaxies': len(data), 'audit_sha256': sha(args.audit),
                   'per_galaxy': details, 'script_sha256': sha(__file__)}, f, indent=2)
        f.write('\n')
    with (out / 'SHA256SUMS.txt').open('w') as f:
        for name in ('calibration.tsv', 'calibration.json'):
            f.write(f'{sha(out / name)}  {name}\n')
    print('OUTPUT=', out, 'GAIN_STATUS=UNVERIFIED', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', default='.')
    p.add_argument('--audit', default='outputs/step1_catalogue_1427444/step1_segmap_audit.tsv')
    p.add_argument('--ids', default=','.join(IDS))
    p.add_argument('--image-hdu', type=int, default=1)
    p.add_argument('--seg-hdu', type=int, default=1)
    p.add_argument('--dilation', type=int, default=10)
    p.add_argument('--acf-lag', type=int, default=5)
    p.add_argument('--rscript', default='/home/uqpdas/.conda/envs/waves_r/bin/Rscript')
    p.add_argument('--skip-gain', action='store_true', help='Synthetic smoke tests only')
    p.add_argument('--out', required=True)
    args = p.parse_args()
    if args.dilation < 1 or args.acf_lag < 1 or args.image_hdu < 0 or args.seg_hdu < 0:
        p.error('Invalid dilation, lag or HDU index')
    run(args)


if __name__ == '__main__':
    main()
