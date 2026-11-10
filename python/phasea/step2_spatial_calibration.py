#!/usr/bin/env python3
"""Three-object spatial background diagnostics; never an observed-fit approval.

Run measurement and synthetic tests through SLURM. Science images and original
segmentations are opened read-only. A fitted plane is diagnostic only: it never
replaces an input image or changes a scientific gate. No gain is estimated here.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
from scipy import fft
from scipy.ndimage import binary_dilation


IDS = ('11690092907692', '11690105411641', '41691060304604')
PROVENANCE_KEYS = ('BUNIT', 'MAGZERO', 'MAGZP', 'PHOTZP', 'GAIN', 'EXPTIME',
                   'DMAG', 'FLXSCALE', 'SOFTNAME', 'SOFTVERS', 'ORIGIN',
                   'HISTORY', 'COMMENT')
STATUS = 'DIAGNOSTIC_ONLY_UNITS_GAIN_NOISE_PROVENANCE_UNVERIFIED'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def array_sha256(a):
    a = np.ascontiguousarray(a)
    return hashlib.sha256(a.tobytes()).hexdigest()


def robust(a):
    a = np.asarray(a, dtype=float)
    median = float(np.median(a))
    return median, float(1.4826 * np.median(np.abs(a - median)))


def sky_mask(image, labels, radius):
    if image.ndim != 2 or image.shape != labels.shape:
        raise ValueError('image/segmentation must share a 2D grid')
    if not np.issubdtype(labels.dtype, np.integer) or np.any(labels < 0):
        raise ValueError('original labels must be nonnegative integers')
    if radius < 1:
        raise ValueError('positive sky-exclusion dilation required')
    # scipy's default connectivity is the same Manhattan dilation as stage2.
    mask = ((labels == 0) & np.isfinite(image) &
            ~binary_dilation(labels != 0, iterations=radius))
    if np.count_nonzero(mask) < 1000:
        raise ValueError('insufficient uncontaminated sky pixels')
    return mask


def cell_statistics(image, mask, side, minimum_fraction=0.4):
    """Equal-weight, robust local summaries, including masked coverage."""
    ny, nx = image.shape
    rows = []
    for y0 in range(0, ny, side):
        for x0 in range(0, nx, side):
            y1, x1 = min(y0 + side, ny), min(x0 + side, nx)
            values = image[y0:y1, x0:x1][mask[y0:y1, x0:x1]]
            fraction = values.size / ((y1 - y0) * (x1 - x0))
            accepted = values.size >= 200 and fraction >= minimum_fraction
            median, sigma = robust(values) if values.size else (None, None)
            rows.append(dict(cell_side=side, grid_x=x0 // side, grid_y=y0 // side,
                             x=(x0 + x1 - 1) / 2, y=(y0 + y1 - 1) / 2,
                             xnorm=2 * ((x0 + x1 - 1) / 2) / (nx - 1) - 1,
                             ynorm=2 * ((y0 + y1 - 1) / 2) / (ny - 1) - 1,
                             n=int(values.size), fraction=fraction,
                             median=median, sigma_mad=sigma,
                             accepted=accepted,
                             fold=int((x0 + x1 - 1) / 2 >= nx / 2) +
                                  2 * int((y0 + y1 - 1) / 2 >= ny / 2)))
    return rows


def fit_plane(design, values):
    """Huber IRLS on cell medians; large contaminated cells get less weight."""
    design, values = np.asarray(design), np.asarray(values)
    if len(values) < 6 or np.linalg.matrix_rank(design) < 3:
        raise ValueError('too few well-distributed cells to fit a plane')
    coefficient = np.linalg.lstsq(design, values, rcond=None)[0]
    for _ in range(20):
        residual = values - design @ coefficient
        centre, scale = robust(residual)
        if scale <= np.finfo(float).eps:
            break
        weight = np.minimum(1, 1.345 * scale /
                            np.maximum(np.abs(residual - centre), scale * 1e-12))
        weight = np.sqrt(weight)
        new = np.linalg.lstsq(design * weight[:, None], values * weight,
                             rcond=None)[0]
        if np.allclose(new, coefficient, rtol=1e-9, atol=1e-14):
            coefficient = new
            break
        coefficient = new
    return coefficient


def plane_diagnostic(rows, rng, bootstrap=300):
    good = [r for r in rows if r['accepted']]
    design = np.asarray([[1., r['xnorm'], r['ynorm']] for r in good]).reshape(-1, 3)
    values = np.asarray([r['median'] for r in good])
    folds = np.asarray([r['fold'] for r in good])
    rank = lambda a: int(np.linalg.matrix_rank(a)) if len(a) else 0
    coverage = dict(ncells=len(good), design_rank=rank(design),
                    requirements=dict(minimum_cells=6, design_rank=3,
                                      minimum_training_cells=6, minimum_heldout_cells=3),
                    quadrant_coverage=[])
    reasons = []
    if len(good) < 6 or coverage['design_rank'] < 3:
        reasons.append('full design has fewer than six cells or rank below three')
    for fold in range(4):
        train, test = folds != fold, folds == fold
        counts = dict(fold=fold, ntrain=int(train.sum()), ntest=int(test.sum()),
                      train_rank=rank(design[train]), test_rank=rank(design[test]))
        coverage['quadrant_coverage'].append(counts)
        if counts['ntest'] < 3 or counts['ntrain'] < 6 or counts['train_rank'] < 3:
            reasons.append(f'quadrant {fold} lacks held-out or training spatial coverage')
    if reasons:
        return dict(status='UNAVAILABLE_INSUFFICIENT_SPATIAL_COVERAGE',
                    reason='; '.join(reasons), **coverage)
    coefficient = fit_plane(design, values)
    constant_prediction = np.empty(len(good))
    plane_prediction = np.empty(len(good))
    fold_rows = []
    for fold in range(4):
        train, test = folds != fold, folds == fold
        c = float(np.median(values[train]))
        p = fit_plane(design[train], values[train])
        constant_prediction[test] = c
        plane_prediction[test] = design[test] @ p
        fold_rows.append(dict(fold=fold, ntrain=int(train.sum()),
                              ntest=int(test.sum()), constant=c,
                              plane_coefficients=p.tolist(),
                              rmse_constant=float(np.sqrt(np.mean((values[test] - c)**2))),
                              rmse_plane=float(np.sqrt(np.mean((values[test] - plane_prediction[test])**2)))))
    e0, e1 = values - constant_prediction, values - plane_prediction
    # Resample 2x2-cell spatial clusters, not individual pixels. Repeat at three
    # cell scales; no claim is made that these blocks remove all dependence.
    group = np.asarray([(r['grid_y'] // 2) * 10000 + r['grid_x'] // 2 for r in good])
    keys = np.unique(group)
    indices = [np.flatnonzero(group == k) for k in keys]
    differences = []
    for _ in range(bootstrap):
        sample = np.concatenate([indices[i] for i in rng.integers(len(keys), size=len(keys))])
        differences.append(np.sqrt(np.mean(e0[sample]**2)) -
                           np.sqrt(np.mean(e1[sample]**2)))
    return dict(status='MEASURED', **coverage, coefficients=coefficient.tolist(),
                model='b0 + bx*xnorm + by*ynorm; normalized coordinates [-1,1]',
                constant=float(np.median(values)),
                cell_median_range=[float(np.min(values)), float(np.max(values))],
                cell_median_sigma_mad=robust(values)[1],
                heldout_quadrants=fold_rows,
                heldout_rmse_constant=float(np.sqrt(np.mean(e0**2))),
                heldout_rmse_plane=float(np.sqrt(np.mean(e1**2))),
                heldout_mae_constant=float(np.mean(np.abs(e0))),
                heldout_mae_plane=float(np.mean(np.abs(e1))),
                rmse_improvement_block_bootstrap_percentiles_2p5_50_97p5=
                np.percentile(differences, [2.5, 50, 97.5]).tolist(),
                bootstrap_clusters=len(keys), bootstrap_samples=bootstrap,
                interpretation='Descriptive held-out prediction. Clusters may remain correlated; '
                               'a plane is diagnostic only and is not an approved sky subtraction.')


def plane_image(shape, coefficient):
    y, x = np.ogrid[-1:1:complex(0, shape[0]), -1:1:complex(0, shape[1])]
    return coefficient[0] + coefficient[1] * x + coefficient[2] * y


def pairwise_acf(image, mask, lag=64):
    """Masked linear correlation with per-lag means AND variances.

    rho(h) = Cov(X[p+h],X[p] | both blank) / sqrt(Var(A) Var(B)).
    ACF covariance is in squared native image units. Zero padding prevents
    opposite-edge wraparound. FFT moment sums are tested against direct pairs.
    """
    if image.shape != mask.shape or image.ndim != 2 or not mask.any():
        raise ValueError('invalid image/mask')
    if lag < 0 or lag >= min(image.shape):
        raise ValueError('lag must be smaller than both image dimensions')
    centre = float(np.mean(image[mask]))
    x = np.where(mask, image - centre, 0.)
    m = mask.astype(float)
    shape = tuple(fft.next_fast_len(2 * n - 1) for n in image.shape)
    fm, fx, fx2 = (fft.rfftn(a, s=shape, workers=1) for a in (m, x, x*x))
    index = np.ix_(np.arange(-lag, lag + 1) % shape[0],
                   np.arange(-lag, lag + 1) % shape[1])

    def corr(a, b):
        return fft.irfftn(a * b.conjugate(), s=shape, workers=1)[index].copy()

    count = np.rint(corr(fm, fm)).astype(np.int64)
    denominator = np.maximum(count, 1)
    sa, sa2, sab = corr(fx, fm), corr(fx2, fm), corr(fx, fx)
    sb, sb2 = sa[::-1, ::-1], sa2[::-1, ::-1]
    covariance = sab / denominator - sa * sb / denominator**2
    va = np.maximum(sa2 / denominator - (sa / denominator)**2, 0.)
    vb = np.maximum(sb2 / denominator - (sb / denominator)**2, 0.)
    with np.errstate(divide='ignore', invalid='ignore'):
        rho = covariance / np.sqrt(va * vb)
    enough = count >= 100
    covariance[~enough], rho[~enough] = np.nan, np.nan
    return dict(dx=np.arange(-lag, lag + 1), dy=np.arange(-lag, lag + 1),
                pairs=count, covariance=covariance, rho=rho,
                pair_mean_a=sa / denominator + centre,
                pair_mean_b=sb / denominator + centre,
                pair_variance_a=va, pair_variance_b=vb)


def direct_lag(image, mask, dx, dy):
    ny, nx = image.shape
    ay, by = slice(max(dy, 0), min(ny, ny + dy)), slice(max(-dy, 0), min(ny, ny - dy))
    ax, bx = slice(max(dx, 0), min(nx, nx + dx)), slice(max(-dx, 0), min(nx, nx - dx))
    valid = mask[ay, ax] & mask[by, bx]
    a, b = image[ay, ax][valid], image[by, bx][valid]
    if len(a) < 100:
        return dict(dx=dx, dy=dy, pairs=int(len(a)), rho=None, covariance=None)
    covariance = float(np.mean((a - a.mean()) * (b - b.mean())))
    denominator = float(np.std(a) * np.std(b))
    return dict(dx=dx, dy=dy, pairs=int(len(a)), covariance=covariance,
                rho=covariance / denominator if denominator > 0 else None)


def acf_summary(acf):
    dx, dy = np.meshgrid(acf['dx'], acf['dy'])
    radius = np.hypot(dx, dy)
    rings = []
    for low, high in ((1, 2), (2, 4), (4, 8), (8, 16), (16, 32), (32, 65)):
        values = acf['rho'][(radius >= low) & (radius < high)]
        values = values[np.isfinite(values)]
        rings.append(dict(rmin=low, rmax=high,
                          median_rho=float(np.median(values)) if len(values) else None,
                          n_lags=len(values)))
    return rings


def variance_bootstrap(values, groups, rng, samples):
    """Cluster bootstrap from sufficient moments, bounded for million pixels."""
    _, inverse = np.unique(groups, return_inverse=True)
    count = np.bincount(inverse)
    total = np.bincount(inverse, weights=values)
    squared = np.bincount(inverse, weights=values**2)
    result = []
    for _ in range(samples):
        take = rng.integers(len(count), size=len(count))
        n, s, s2 = count[take].sum(), total[take].sum(), squared[take].sum()
        result.append(max(0., (s2 - s*s/n) / (n-1)))
    return np.percentile(result, [2.5, 50, 97.5]).tolist()


def aperture_statistics(image, mask, acf, rng, sides=(1, 2, 4, 8, 16, 32, 64), bootstrap=300):
    rows = []
    ny, nx = image.shape
    pixel_variance = float(np.var(image[mask], ddof=1))
    for side in sides:
        sy, sx = ny // side, nx // side
        patch = np.where(mask, image, 0.)[:sy*side, :sx*side].reshape(sy, side, sx, side)
        whole = mask[:sy*side, :sx*side].reshape(sy, side, sx, side).all(axis=(1, 3))
        sums = patch.sum(axis=(1, 3))[whole]
        yy, xx = np.indices((sy, sx))
        y, x = (yy[whole] + .5)*side, (xx[whole] + .5)*side
        if len(sums) < 20:
            rows.append(dict(side=side, n_apertures=len(sums), status='INSUFFICIENT_BLANK_APERTURES'))
            continue
        variance = float(np.var(sums, ddof=1))
        groups = (y // 128).astype(int)*10000 + (x // 128).astype(int)
        groups2 = (y // 256).astype(int)*10000 + (x // 256).astype(int)
        dx, dy = np.meshgrid(acf['dx'], acf['dy'])
        weight = np.maximum(side - np.abs(dx), 0) * np.maximum(side - np.abs(dy), 0)
        valid = weight > 0
        missing_lags = int(np.count_nonzero(~np.isfinite(acf['covariance'][valid])))
        support_complete = (acf['dx'][0] <= 1-side and acf['dx'][-1] >= side-1 and
                            acf['dy'][0] <= 1-side and acf['dy'][-1] >= side-1)
        if not support_complete:
            predicted, prediction_status = None, 'UNAVAILABLE_ACF_LAG_SUPPORT_TOO_SMALL'
        elif missing_lags:
            predicted = None
            prediction_status = 'UNAVAILABLE_REQUIRED_LAGS_HAVE_INSUFFICIENT_PAIRS'
        else:
            predicted = float(np.sum(weight[valid] * acf['covariance'][valid]))
            prediction_status = 'MEASURED_CONDITIONAL_ON_STATIONARITY'
        local = []
        for quadrant in range(4):
            use = ((x >= nx/2).astype(int) + 2*(y >= ny/2).astype(int)) == quadrant
            local.append(dict(quadrant=quadrant, n=int(use.sum()),
                              mean=float(sums[use].mean()) if use.any() else None,
                              variance=float(np.var(sums[use], ddof=1)) if use.sum() > 1 else None))
        rows.append(dict(side=side, area=side**2, n_apertures=int(len(sums)), status='DESCRIPTIVE',
                         variance=variance, mean=float(sums.mean()),
                         independent_pixel_prediction=side**2*pixel_variance,
                         variance_ratio_to_independent=variance/(side**2*pixel_variance),
                         covariance_sum_prediction_conditional_stationarity=predicted,
                         covariance_prediction_status=prediction_status,
                         covariance_prediction_missing_required_lags=missing_lags,
                         variance_bootstrap_block128_percentiles=variance_bootstrap(sums, groups, rng, bootstrap),
                         variance_bootstrap_block256_percentiles=variance_bootstrap(sums, groups2, rng, bootstrap),
                         quadrants=local,
                         interpretation='Nonoverlapping fully blank square apertures. Bootstrap blocks '
                         'may remain dependent. Covariance prediction requires stationarity and '
                         'matching sky support; no effective independent pixel count is certified.'))
    return rows


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def write_tsv(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def provenance(root):
    sources = []
    for path, excerpt in (
        (Path('/home/uqpdas/WAVES_BD_pipeline/measure_survey_gain.py'),
         'Chooses median of per-cutout profoundGainEst values from run.log; '
         'not original detector/coadd calibration.'),
        (Path('/home/uqpdas/WAVES_BD_pipeline/profuse_helpers.R'),
         'DEFAULT_MAGZERO is assumed 23.9. Header absence is explicitly acknowledged. '
         'robust_sky_gain accepts fixed override or estimates/clamps histogram-based gain.'),
        (root/'config/step0.yaml', 'Reported settings remain REPORTED, not measured calibration.'),
        (root/'R/sky_noise.R', 'Model variance convention: skyRMS^2 + max(model,0)/gain.'),
    ):
        item = dict(path=str(path), exists=path.is_file(), interpretation=excerpt)
        if path.is_file():
            item['sha256'] = sha256(path)
            lines = path.read_text().splitlines()
            if path.name == 'measure_survey_gain.py':
                item['relevant_lines'] = lines[:20] + [s for s in lines if 'median(' in s or 'options(bd.gain.fixed' in s]
            elif path.name == 'profuse_helpers.R':
                item['relevant_lines'] = lines[22:38] + lines[210:239]
        sources.append(item)
    return dict(local_source_identity=sources,
                producer_source_commit='UNAVAILABLE: downstream fit code does not identify image production',
                required_upstream=['Original WAVES 0.30arcsec tile identity/checksum and producing commit/config',
                                   'Photometric scale, extinction/DMAG convention and flux-per-pixel versus surface brightness',
                                   'Exposure counts, gains, flux scales, weights, coadd operator and resampling kernel',
                                   'Sky subtraction and variance/weight-map definitions; independent exposure splits if available'],
                gain_equations=['If y=a*x and Var_Poisson(x)=E[x]/G_x, then G_y=G_x/a.',
                                'For independent calibrated exposures E[f_e]=F, '
                                'y=sum(alpha_e*f_e), sum(alpha_e)=1: 1/G_y=sum(alpha_e^2/G_e).',
                                'General resampling y_p=sum_j A_pj*x_j gives '
                                'Cov(y_p,y_q)=sum_j A_pj*A_qj*Var(x_j). '
                                'Source-dependent variance can require a spatial/nonlocal operator rather than one scalar gain.'],
                primary_sources=[
                    dict(url='https://kids.strw.leidenuniv.nl/DR4/format.php',
                         applicability='Native KiDS conventions (0.2arcsec,ZP0,GAIN,DMAG) do not certify transformed WAVES cutouts.'),
                    dict(url='https://wavesurvey.org/input-catalogues/data-viewer/',
                         applicability='Official tile/cutout service entry; exact delivered-product provenance still required.'),
                    dict(url='https://github.com/asgr/ProFound/blob/1d56e383159606b7df9e81dcac2bf6e0051558ec/R/profoundUtility.R#L169-L192',
                         applicability='Pinned gain estimator uses sky histogram alone, not detector/exposure metadata.')],
                gain_status='UNVERIFIED; no scalar687 certification', units_zeropoint_status='UNVERIFIED',
                upstream_sky_status='UNVERIFIED_FROM_PROVENANCE; residuals measured separately',
                independent_pixel_count_status='NOT_CERTIFIED_FOR_BIC_OR_LIKELIHOOD')


def diagnostic_plot(path, image, mask, cells, plane, acf_raw, acf_detrended):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    good = [r for r in cells if r['accepted']]
    if good:
        sc = axes[0, 0].scatter([r['x'] for r in good], [r['y'] for r in good],
                                c=[r['median'] for r in good], marker='s', s=65)
        fig.colorbar(sc, ax=axes[0, 0], label='Native image units')
    else:
        axes[0, 0].text(.5, .5, 'No cells meet the original coverage criteria',
                        ha='center', transform=axes[0, 0].transAxes)
    axes[0, 0].set_title('128px sky-cell medians (diagnostic)')
    med, sigma = robust(image[mask])
    axes[0, 1].imshow(np.ma.masked_where(~mask, image), origin='lower',
                      vmin=med-2*sigma, vmax=med+2*sigma, interpolation='nearest')
    axes[0, 1].set_title('Sky pixels; source labels unchanged')
    for ax, acf, title in ((axes[1, 0], acf_raw, 'Raw sky pair-normalized correlation'),
                            (axes[1, 1], acf_detrended, 'After diagnostic plane removal')):
        if acf is None:
            ax.text(.5, .5, 'Unavailable: insufficient spatial coverage\nfor the 128px reference plane',
                    ha='center', transform=ax.transAxes)
            ax.set_axis_off()
            continue
        lag = int(acf['dx'][-1])
        im = ax.imshow(acf['rho'], origin='lower', extent=(-lag,lag,-lag,lag),
                       vmin=-.15, vmax=.3, cmap='coolwarm')
        ax.set_title(title); ax.set_xlabel('dx (pixels)'); ax.set_ylabel('dy (pixels)')
        fig.colorbar(im, ax=ax)
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)


def run(args):
    from astropy.io import fits
    import astropy
    import scipy
    if not __import__('os').environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Measurement requires a SLURM allocation')
    root, out = Path(args.repo).resolve(), Path(args.out).resolve()
    if not (root/'.hector_bd_root').is_file() or root/'outputs' not in out.parents or out.exists():
        raise ValueError('Require a NEW output directory below HectorBD outputs/')
    if args.lag < 64 or args.lag > 128 or not 100 <= args.bootstrap <= 2000:
        raise ValueError('lag must be64..128; bootstrap100..2000')
    with open(args.audit, newline='', encoding='utf-8-sig') as f:
        audit = list(csv.DictReader(f, delimiter='\t'))
    chosen = []
    for gid in IDS:
        matches = [r for r in audit if r['uberID'] == gid]
        if len(matches) != 1:
            raise ValueError(f'{gid}: audit must contain exactly one row')
        row = matches[0]
        if row['status'] != 'IDENTICAL' or float(row['wcs_max_sep_arcsec']) > .05:
            raise ValueError(f'{gid}: requires identical original labels and verified WCS')
        chosen.append(row)
    out.mkdir(parents=True)
    rng = np.random.default_rng(args.seed)
    report = dict(status=STATUS, ids=list(IDS), script_sha256=sha256(__file__),
                  audit_sha256=sha256(args.audit), provenance=provenance(root),
                  runtime=dict(python=sys.version, executable=sys.executable,
                               numpy=np.__version__, scipy=scipy.__version__, astropy=astropy.__version__,
                               platform=platform.platform(), slurm_job_id=__import__('os').environ['SLURM_JOB_ID']),
                  settings=vars(args), galaxies=[])
    report['repository_commit'] = subprocess.check_output(
        ['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    report['wrapper_sha256'] = sha256(root/'slurm/step2_spatial_calibration.sbatch')
    for row in chosen:
        gid = row['uberID']; folder = out/gid; folder.mkdir()
        ip, sp = Path(row['image_path']), Path(row['segmap_path'])
        before = {str(p): sha256(p) for p in (ip, sp)}
        with fits.open(ip, mode='readonly', memmap=False) as hdus:
            image = np.asarray(hdus[1].data, float)
            headers = [h.header.copy() for h in hdus]
        with fits.open(sp, mode='readonly', memmap=False) as hdus:
            labels = np.asarray(hdus[1].data)
            seg_headers = [h.header.copy() for h in hdus]
        if image.shape != (1333, 1333) or labels.shape != image.shape:
            raise ValueError(f'{gid}: unexpected original cutout dimensions')
        label_identity = array_sha256(labels)
        with (folder/'headers.txt').open('w') as f:
            for product, hs in (('science', headers), ('original_segmentation', seg_headers)):
                for index, header in enumerate(hs):
                    f.write(f'{product} HDU {index}\n'+header.tostring(sep='\n', padding=False)+'\n\n')
        metadata = {k: [str(header[k]) for header in headers if k in header] for k in PROVENANCE_KEYS}
        galaxy = dict(uberID=gid, image_path=str(ip), segmap_path=str(sp), hdus=dict(image=1, segmap=1),
                      source_sha256=before, labels_array_sha256=label_identity,
                      header_calibration_fields=metadata,
                      missing_calibration_fields=[k for k,v in metadata.items() if not v], masks=[])
        for radius in (10, 20, 40):
            mask = sky_mask(image, labels, radius)
            median, sigma = robust(image[mask])
            if not np.isfinite(sigma) or sigma <= 0:
                raise ValueError('invalid robust sky scale')
            print(gid, 'radius', radius, 'n_sky', int(mask.sum()), 'median', median, 'sigma_MAD', sigma, flush=True)
            cells = {side: cell_statistics(image, mask, side) for side in (64, 128, 256)}
            planes = {side: plane_diagnostic(cells[side], rng, args.bootstrap) for side in cells}
            for side in cells:
                write_tsv(folder/f'cells_d{radius}_s{side}.tsv', cells[side])
            plane_available = planes[128]['status'] == 'MEASURED'
            plane = plane_image(image.shape, planes[128]['coefficients']) if plane_available else None
            acfs, fields = {}, {}
            measurements = {'raw': image - median}
            if plane_available:
                measurements['plane_detrended'] = image - plane
            else:
                fields['plane_detrended'] = dict(status=planes[128]['status'],
                                                 reason=planes[128]['reason'])
            for name, field in measurements.items():
                acf = pairwise_acf(field, mask, args.lag)
                acfs[name] = acf
                np.savez_compressed(folder/f'acf_d{radius}_{name}.npz', **acf)
                local = []
                for q in range(4):
                    ny, nx = image.shape
                    sy, sx = slice((q//2)*ny//2, ((q//2)+1)*ny//2), slice((q%2)*nx//2, ((q%2)+1)*nx//2)
                    local.append(dict(quadrant=q, lags=[direct_lag(field[sy,sx], mask[sy,sx], dx, dy)
                         for dx,dy in ((1,0),(0,1),(2,0),(0,2),(5,0),(0,5),(16,0),(0,16),(64,0),(0,64),(64,64))]))
                fields[name] = dict(status='MEASURED', acf_rings=acf_summary(acf), local_quadrant_lags=local,
                                    pixel_mean=float(field[mask].mean()),
                                    pixel_variance=float(np.var(field[mask], ddof=1)),
                                    apertures=aperture_statistics(field, mask, acf, rng, bootstrap=args.bootstrap))
            saved = dict(sky_mask=mask, plane_available=np.asarray(plane_available),
                         plane_reference_cell_side=np.asarray(128))
            if plane_available:
                saved['plane'] = plane
            np.savez_compressed(folder/f'diagnostic_sky_d{radius}.npz', **saved)
            diagnostic_plot(folder/f'spatial_d{radius}.png', image, mask, cells[128], plane,
                            acfs['raw'], acfs.get('plane_detrended'))
            galaxy['masks'].append(dict(dilation=radius, sky_mask_sha256=array_sha256(mask),
                                       n_sky=int(mask.sum()), median=median, sigma_mad=sigma,
                                       plane_available=plane_available, plane_reference_cell_side=128,
                                       planes=planes, fields=fields))
        after = {str(p): sha256(p) for p in (ip, sp)}
        if before != after or label_identity != array_sha256(labels):
            raise RuntimeError('original science/segmentation changed during diagnostic')
        galaxy['source_unchanged'] = True
        write_json(folder/'spatial.json', galaxy)
        report['galaxies'].append(galaxy)
    write_json(out/'spatial_calibration.json', report)
    (out/'INTERPRETATION.txt').write_text(
        'DIAGNOSTIC ONLY. No original image/label or Step0 status was changed.\n'
        'Compare mask dilation, cell scale, held-out quadrants and raw/detrended covariance.\n'
        'A plane fitted to sky cells is a diagnostic comparison, not an approved subtraction.\n'
        'Bootstrap blocks can retain correlations; disagreement across block/region sizes is evidence of nonstationarity.\n'
        'Aperture factors describe blank square sums only, not universal Neff, BIC corrections or galaxy likelihoods.\n'
        'Units/zeropoint/upstream sky/gain remain unverified without exact image-production provenance.\n'
        'No source-Poisson gain, including687, was certified from a sky histogram.\n')
    paths = sorted(p for p in out.rglob('*') if p.is_file())
    (out/'SHA256SUMS.txt').write_text(''.join(f'{sha256(p)}  {p.relative_to(out)}\n' for p in paths))
    print('SPATIAL_CALIBRATION_DIAGNOSTICS_COMPLETE', 'GAIN_UNITS_NEFF_UNVERIFIED', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', default='.')
    parser.add_argument('--audit', default='outputs/step1_catalogue_1427444/step1_segmap_audit.tsv')
    parser.add_argument('--out', required=True)
    parser.add_argument('--lag', type=int, default=64)
    parser.add_argument('--bootstrap', type=int, default=300)
    parser.add_argument('--seed', type=int, default=20261108)
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
