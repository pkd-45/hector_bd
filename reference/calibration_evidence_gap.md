# WAVES delivered-image calibration evidence gap

Status: **UNVERIFIED; no observed-fit authorization.** This records the limits
of the available evidence, not a replacement calibration. Original images,
segmentations and historical products remain read-only.

## Evidence available on Friday

The three controls are `11690092907692`, `11690105411641` and
`41691060304604`. Their science images are
`/data/uqpdas/WAVES/<ID>/<ID>_r_1333.fits.gz`, HDU 1 (zero based).
The original delivered label maps are adjacent
`<ID>_waves_postprocessed_1333.fits.gz`, HDU 1. File hashes and the existing
measurements are in `outputs/step2_calibration_1437317/calibration.json`.

The supplied headers establish the grid, including 0.30 arcsec pixels, but
do not establish units, photometric zeropoint, exposure, effective gain,
producer version or upstream sky processing. Each inspected directory contains
only these two FITS products. The Mac location `/Volumes/waves/WAVES/tiles/`
is not mounted on Friday; `/data/uqpdas/WAVES/tiles/` is absent.

`/home/uqpdas/WAVES_BD_pipeline/profuse_helpers.R` explicitly makes its 23.9
default conditional on the cutouts already having the assumed flux units.
`measure_survey_gain.py` in that directory recommends a median of downstream
`profoundGainEst` estimates. Neither supplies original production calibration.
The three coarse estimates (1162.48, 783.88, 405.71) cannot certify 687.
The exact pinned ProFound implementation and its documentation describe an
approximate sky-histogram estimate, not an instrument/coadd calibration:
[implementation](https://github.com/asgr/ProFound/blob/1d56e383159606b7df9e81dcac2bf6e0051558ec/R/profoundUtility.R#L169-L192),
[documentation](https://github.com/asgr/ProFound/blob/1d56e383159606b7df9e81dcac2bf6e0051558ec/man/profoundGainEst.Rd).

[Native KiDS DR4 documentation](https://kids.strw.leidenuniv.nl/DR4/format.php)
describes a 0.2 arcsec product with its own zeropoint, gain and DMAG conventions.
Those conventions do not authenticate this later 0.30 arcsec WAVES product.

A historical source manifest records a possible independent Hector catalogue
at `/Users/uqpdas/Downloads/Hector_input_catalogue_v03.2.csv`, SHA256
`109ebe882c3e66c91773757f3dc3687ab581891d77c6ba31fccc79d29ec85dd1`.
It is not available at that path on Friday. The staged 832-object catalogue
export contains coordinates, mass and position angle, but no independent
photometry. A catalogue magnitude would also need its aperture, bandpass,
extinction and PSF conventions before a flux comparison could calibrate images.

## Information required from Sabine Bellstedt / WAVES production

No external message has been sent. The minimum production response should bind
the following information to the delivered files or their exact parent tiles:

1. Parent tile/exposure identifiers, release, producing scripts/configuration
   and version; the cutout/remapping operation and any subsequent scaling.
2. Pixel units (flux per pixel, surface brightness, counts or count rate), the
   exact AB magnitude conversion, zeropoint, extinction corrections and any
   spatial photometric corrections. Confirm how pixel area enters the scaling.
3. The meaning and units of gain, native detector gains/exposure times, coadd
   weights and their normalization, rejection/clipping policy and coverage maps.
   State whether effective gain is a scalar, tile value or spatial map.
4. Resampling kernel, flux-conservation convention, pixel scale transforms and
   associated variance/covariance propagation or noise-power calibration.
5. Sky subtraction before and after coaddition: method, mesh/smoothing scales,
   source masks, retained sky/variance maps and known residual-background terms.
6. Whether delivered weights/variance maps include source Poisson noise, sky,
   read noise and correlated/background uncertainty, and their exact units.

Independent repeat exposures/splits or variance products in final image units
can test the source-Poisson coefficient. A single coadd sky histogram cannot
separate it from background, resampling and other noise contributions.

## Conventions that must be checked, not assumed

For a linear flux scaling `y = a*x`, a gain expressed as electrons per image
unit transforms as `G_y = G_x/a`. For normalized exposure weights `alpha_e`
and exposures already in common flux units, independent source noise would give
`1/G_coadd = sum(alpha_e^2/G_e)`. Resampling additionally creates off-diagonal
source covariance. These equations require the actual exposure/scaling/weight
metadata and do not justify transferring a native KiDS gain to a WAVES cutout.

The pipeline's diagonal variance `skyRMS^2 + max(model,0)/gain` is a conditional
model. Its applicability requires validated source-noise scaling and evidence
about spatial background and covariance on the fitted support. Sky trends,
raw/detrended autocovariance, blank-aperture variance and spatial block
uncertainty must be examined before deriving any effective independent-pixel
count. Such a count is support/scale dependent and cannot be certified from a
five-pixel ACF window or treated as a universal BIC correction.

All corresponding Step 0 calibration and noise gates remain unverified until
these requirements are met. Synthetic recovery, historical P metadata and
pilot preparation can proceed independently without authorizing an observed fit.
