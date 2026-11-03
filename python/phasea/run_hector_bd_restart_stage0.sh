#!/usr/bin/env bash
# Evidence collection only: no fitting, package installation, or source edits.
# Run this on the Mac that already has the working SSH alias "friday".
set -euo pipefail

if [[ "${1:-}" == "--help" ]]; then
  printf '%s\n' 'Usage: bash run_hector_bd_restart_stage0.sh [SSH_ALIAS]' \
    'Default SSH alias: friday. Creates one new local evidence log in Downloads.'
  exit 0
fi
if [[ $# -gt 1 ]]; then
  printf '%s\n' 'ERROR: expected at most one SSH alias.' >&2
  exit 2
fi
bd_audit_host=${1:-friday}
case "$bd_audit_host" in
  -*|*[!a-zA-Z0-9_.@-]*|'') printf '%s\n' 'ERROR: invalid SSH alias.' >&2; exit 2 ;;
esac
command -v ssh >/dev/null 2>&1 || { printf '%s\n' 'ERROR: ssh is missing.' >&2; exit 2; }
bd_audit_parent=${BD_AUDIT_OUTPUT_PARENT:-"$HOME/Downloads"}
mkdir -p "$bd_audit_parent"
bd_audit_out=$(mktemp -d "$bd_audit_parent/HectorBD_RESTART_STAGE0.XXXXXXXX")
printf 'EVIDENCE_DIRECTORY=%s\n' "$bd_audit_out"
printf '%s\n' 'Collecting existing files and metadata on Friday. Historical directories are read only.'

if ssh "$bd_audit_host" 'bash -s' <<'BD_REMOTE' 2>&1 | tee "$bd_audit_out/FRIDAY_STAGE0_EVIDENCE.txt"
set -euo pipefail
bd_stage_tmp=$(mktemp -d "${TMPDIR:-/tmp}/hector_bd_stage0.XXXXXXXX")
trap 'rm -rf -- "$bd_stage_tmp"' EXIT
bd_stage_root=/data/scratch/uqpdas/WAVES_bd_v32_production842/baseline/FIXED4_seed666
bd_stage_code=/data/uqpdas/WAVES_v32
bd_stage_rlib=/data/uqpdas/WAVES_v38_profit272_pilot/Rlib
bd_stage_source_ok=0
printf 'COLLECTION_MODE=READ_ONLY_HISTORICAL_INPUTS\n'
printf 'COLLECTION_TIME_UTC=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
printf 'PRODUCTION_ROOT=%s\nCODE_ROOT=%s\n' "$bd_stage_root" "$bd_stage_code"

# This is source inspection, not source execution. Capture the runner and its
# immediate R/Python/shell/config neighbours, where the missing helpers may live.
printf '\nSOURCE_DIRECTORY_INVENTORY\n'
if [[ -d "$bd_stage_code" ]]; then
  ls -la "$bd_stage_code"
  if [[ -f "$bd_stage_code/4_run_production_v32.R" ]]; then
    bd_stage_source_ok=1
  else
    printf '%s\n' 'PRODUCTION_RUNNER=MISSING'
  fi
  shopt -s nullglob
  for bd_stage_file in "$bd_stage_code"/*.R "$bd_stage_code"/*.r "$bd_stage_code"/*.py "$bd_stage_code"/*.sh "$bd_stage_code"/*.json; do
    [[ -f "$bd_stage_file" ]] || continue
    bd_stage_bytes=$(wc -c < "$bd_stage_file")
    printf '\nSOURCE_BEGIN=%s\nSOURCE_BYTES=%s\n' "$bd_stage_file" "$bd_stage_bytes"
    if command -v sha256sum >/dev/null 2>&1; then
      sha256sum "$bd_stage_file"
    elif command -v shasum >/dev/null 2>&1; then
      shasum -a 256 "$bd_stage_file"
    else
      printf '%s\n' 'SOURCE_HASH=UNAVAILABLE'
    fi
    if [[ "$bd_stage_bytes" -le 2000000 ]]; then
      cat "$bd_stage_file"
    else
      printf '%s\n' 'SOURCE_CONTENT=SKIPPED_OVER_2MB'
    fi
    printf '\nSOURCE_END=%s\n' "$bd_stage_file"
  done
else
  printf '%s\n' 'SOURCE_DIRECTORY=MISSING'
fi

bd_stage_r=/home/uqpdas/.conda/envs/waves_r/bin/Rscript
if [[ ! -x "$bd_stage_r" ]]; then
  bd_stage_r=$(command -v Rscript || true)
fi
if [[ -z "$bd_stage_r" ]]; then
  printf '%s\n' 'R_METADATA=UNAVAILABLE_NO_RSCRIPT' 'COLLECTION_STATUS=PARTIAL'
  exit 3
fi
printf '\nRSCRIPT_EXECUTABLE=%s\n' "$bd_stage_r"

cat > "$bd_stage_tmp/inspect.R" <<'BD_R'
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
root <- args[[1L]]
oldlib <- args[[2L]]
paths_file <- args[[3L]]
cat("R_VERSION=", R.version.string, "\n", sep = "")
cat("DEFAULT_LIBRARY_PATHS=", paste(.libPaths(), collapse = ";"), "\n", sep = "")
packages <- c("ProFit", "ProFuse", "ProFound", "Highlander", "LaplacesDemon", "Rfits")
for (lib in unique(c(oldlib, .libPaths()))) {
  for (pkg in packages) {
    desc <- file.path(lib, pkg, "DESCRIPTION")
    cat("\nPACKAGE_DESCRIPTION=", desc, "\n", sep = "")
    if (file.exists(desc)) {
      d <- tryCatch(read.dcf(desc), error = function(e) NULL)
      if (!is.null(d)) {
        fields <- intersect(c("Package", "Version", "Date", "Built", "RemoteType", "RemoteRepo", "RemoteUsername", "RemoteRef", "RemoteSha", "GithubSHA1"), colnames(d))
        print(d[, fields, drop = FALSE])
      } else cat("DESCRIPTION_READ=FAILED\n")
    } else cat("INSTALLED_HERE=NO\n")
  }
}
cat("FINGERPRINT_SCOPE=INSTALLED_FILES_TODAY_NOT_PROOF_OF_HISTORICAL_LOADED_CODE\n")
if (requireNamespace("ProFound", quietly=TRUE)) {
  fx <- tryCatch(getExportedValue("ProFound", "profoundGainEst"), error=function(e) NULL)
  if (is.function(fx)) cat("PROFOUND_GAIN_EST_FORMALS=", paste(names(formals(fx)), collapse=","), "\n", sep="")
  else cat("PROFOUND_GAIN_EST_FORMALS=UNAVAILABLE_IN_ACTIVE_LIBRARY\n")
} else cat("PROFOUND_GAIN_EST_FORMALS=PROFOUND_NOT_AVAILABLE_IN_ACTIVE_LIBRARY\n")
cat("GAIN_EST_EXECUTED=NO_INPUTS_NOT_YET_VERIFIED\n")
files <- sort(Sys.glob(file.path(root, "*", "r_1333_v*", "inputs.rds")))
cat("\nINPUT_RDS_COUNT_ALL_VERSIONS=", length(files), "\n", sep = "")
if (length(files)) {
  id_groups <- split(files, basename(dirname(dirname(files))))
  ambiguous <- id_groups[lengths(id_groups) > 1L]
  cat("MULTIPLE_INPUT_VERSIONS_GALAXIES=", length(ambiguous), "\n", sep="")
  if (length(ambiguous)) for (g in names(ambiguous)) cat("MULTIPLE_INPUT_VERSIONS=", g, " ", paste(ambiguous[[g]], collapse=";"), "\n", sep="")
}
if (length(files)) cat(paste(files, collapse = "\n"), "\n")
ids <- if (length(files)) basename(dirname(dirname(files))) else character()
pilot <- c("11680009812466", "12690079611187", "176980350403709", "42710765002459", "41691060304604")
pilot_idx <- match(pilot, sub("^W", "", ids))
for (j in seq_along(pilot)) {
  if (is.na(pilot_idx[[j]])) cat("REQUESTED_PILOT_MISSING=", pilot[[j]], "\n", sep = "")
}
chosen <- unique(c(pilot_idx[!is.na(pilot_idx)], head(seq_along(files), 2L)))
path_rows <- list()
show_array <- function(v, label) {
  cat("ARRAY=", label, " TYPE=", typeof(v), " DIM=", paste(dim(v), collapse = "x"), " LENGTH=", length(v), "\n", sep = "")
  if (is.numeric(v) || is.logical(v)) {
    cat("FINITE=", sum(is.finite(v)), " NONFINITE=", sum(!is.finite(v)), "\n", sep = "")
    if (is.logical(v)) cat("TRUE=", sum(v, na.rm = TRUE), " NA=", sum(is.na(v)), "\n", sep = "")
    if (is.numeric(v) && any(is.finite(v))) print(range(v[is.finite(v)]))
  }
}
read_one <- function(f) tryCatch(readRDS(f), error = function(e) {
  cat("READ_RDS_FAILED=", f, " MESSAGE=", conditionMessage(e), "\n", sep = "")
  NULL
})
for (i in chosen) {
  f <- files[[i]]
  cat("\nINPUT_BEGIN=", f, "\n", sep = "")
  x <- read_one(f)
  if (!is.list(x)) { cat("INPUT_NOT_A_LIST_OR_UNREADABLE\n"); next }
  str(x, max.level = 2L, list.len = 50L, vec.len = 6L)
  for (key in intersect(c("image_path", "segmap_path"), names(x))) {
    val <- x[[key]]
    cat("EXACT_PATH_FIELD=", key, "\n", sep = "")
    dput(val)
    if (is.character(val) && length(val) == 1L && !is.na(val) && nzchar(val)) {
      path_rows[[length(path_rows) + 1L]] <- data.frame(galaxy = ids[[i]], field = key, path = val, stringsAsFactors = FALSE)
    }
  }
  for (key in intersect(c("target", "band", "side_fit", "magzero"), names(x))) {
    cat("METADATA_FIELD=", key, "\n", sep = "")
    str(x[[key]], max.level = 3L, list.len = 80L, vec.len = 20L)
  }
  for (key in intersect(c("image", "segim", "sigma", "mask", "region"), names(x))) show_array(x[[key]], key)
  for (name in c("fit_FS.rds", "fit_PD.rds", "fit_BD.rds", "psf.rds")) {
    ff <- file.path(dirname(f), name)
    cat("\nSAVED_OBJECT=", ff, "\n", sep = "")
    if (!file.exists(ff)) { cat("SAVED_OBJECT=MISSING\n"); next }
    obj <- read_one(ff)
    str(obj, max.level = 2L, list.len = 45L, vec.len = 4L)
    if (!is.list(obj)) next
    data <- if (is.list(obj$Data)) obj$Data else if (is.list(obj$fit) && is.list(obj$fit$Data)) obj$fit$Data else NULL
    if (is.list(data)) {
      cat("DATA_NAMES=", paste(names(data), collapse = ","), "\n", sep = "")
      for (key in intersect(c("image", "segim", "sigma", "mask", "region", "region_which", "calcregion", "psf"), names(data))) show_array(data[[key]], paste0("Data$", key))
    }
  }
  cat("INPUT_END=", f, "\n", sep = "")
}
paths <- if (length(path_rows)) do.call(rbind, path_rows) else data.frame(galaxy = character(), field = character(), path = character())
write.table(paths, paths_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = TRUE, qmethod = "double", na = "")
cat("\nPATHS_TSV_BEGIN\n")
cat(readLines(paths_file), sep = "\n")
cat("\nPATHS_TSV_END\n")
cat("\nPATH_EXPORT=COMPLETE\n")
cat("SEGMENT_PRODUCT=UNRESOLVED_PENDING_SOURCE_REVIEW\n")
cat("EXACT_SEGMENT_SUPPORT_COMPARISON=NOT_RUN\n")
BD_R

# Parse on the actual R runtime before any RDS is read. No -c quoting through SSH.
"$bd_stage_r" --vanilla -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[1L])); cat("R_SYNTAX_GATE=PASS\n")' "$bd_stage_tmp/inspect.R"
"$bd_stage_r" --vanilla "$bd_stage_tmp/inspect.R" "$bd_stage_root" "$bd_stage_rlib" "$bd_stage_tmp/paths.tsv"

bd_stage_python=/home/uqpdas/.conda/envs/waves_r/bin/python3
if [[ ! -x "$bd_stage_python" ]]; then
  bd_stage_python=$(command -v python3 || true)
fi
if [[ -z "$bd_stage_python" ]]; then
  printf '%s\n' 'FITS_HEADERS=UNAVAILABLE_NO_PYTHON3' 'COLLECTION_STATUS=PARTIAL'
  exit 3
fi
printf '\nPYTHON_EXECUTABLE=%s\n' "$bd_stage_python"
"$bd_stage_python" - "$bd_stage_tmp/paths.tsv" "$bd_stage_source_ok" <<'BD_PY'
import csv
import hashlib
import os
import sys

try:
    from astropy.io import fits
except ImportError as exc:
    print('FITS_HEADERS=UNAVAILABLE_ASTROPY_IMPORT:', exc)
    print('COLLECTION_STATUS=PARTIAL')
    raise SystemExit(3)

with open(sys.argv[1], newline='', encoding='utf-8') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
seen = set()
unavailable = 0
source_ok = sys.argv[2] == '1'
for row in rows:
    original = row['path']
    # Resolve only the recorded absolute path. Do not replace Mac prefixes or
    # choose a similarly named file on Friday as if it had the same lineage.
    print('\nRECORDED_PATH:', row['galaxy'], row['field'], repr(original))
    if not os.path.isabs(original) or not os.path.isfile(original):
        print('ORIGINAL_PATH_STATUS=MISSING_OR_NOT_ABSOLUTE')
        unavailable += 1
        continue
    if original in seen:
        print('HEADER_STATUS=ALREADY_REPORTED')
        continue
    seen.add(original)
    try:
        hasher = hashlib.sha256()
        with open(original, 'rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                hasher.update(block)
        print('FILE_SHA256:', hasher.hexdigest())
        with fits.open(original, mode='readonly', memmap=False) as hdus:
            hdus.info()
            for index, hdu in enumerate(hdus):
                print('HDU_HEADER_BEGIN:', index)
                print(repr(hdu.header))
                print('HDU_HEADER_END:', index)
    except Exception as exc:
        unavailable += 1
        print('HEADER_READ_FAILED:', type(exc).__name__, str(exc))
# These are targeted byte-evidence checks on the recorded original files.
# They do not guess HDUs when more than one plausible image extension exists.
print('\nSTEP0_REPORTED_CONFIGURATION=band:r;grid:1333x1333;scale_arcsec:0.30;unit:microJy;AB_zeropoint:23.9;gain:687')
print('STEP0_REPORTED_NOT_BYTES_VERIFIED=historical_noise_formula,upstream_sky_subtraction,ProFound_gain_687')
paired={}
for item in rows:
    paired.setdefault(item['galaxy'],{})[item['field']]=item['path']
for gid,pair in sorted(paired.items()):
    image_path=pair.get('image_path');seg_path=pair.get('segmap_path')
    print('\nSTEP0_GALAXY_BEGIN:',gid)
    if not image_path or not seg_path or not os.path.isfile(image_path) or not os.path.isfile(seg_path):
        print('STEP0_PAIR=BLOCKED_MISSING_RECORDED_FILE')
        continue
    try:
        import numpy as np
        with fits.open(image_path,mode='readonly',memmap=True) as ih, fits.open(seg_path,mode='readonly',memmap=True) as sh:
            candidates_i=[(k,h) for k,h in enumerate(ih) if getattr(h,'data',None) is not None and h.data.ndim==2]
            candidates_s=[(k,h) for k,h in enumerate(sh) if getattr(h,'data',None) is not None and h.data.ndim==2]
            print('STEP0_IMAGE_2D_HDUS=',[k for k,h in candidates_i])
            print('STEP0_SEG_2D_HDUS=',[k for k,h in candidates_s])
            if len(candidates_i)!=1 or len(candidates_s)!=1:
                print('STEP0_HDU=UNRESOLVED_MULTIPLE_OR_ZERO_2D_HDUS')
                continue
            ki,hi=candidates_i[0];ks,hs=candidates_s[0]
            im=np.asarray(hi.data);seg=np.asarray(hs.data)
            print('STEP0_HDUS=VERIFIED_UNIQUE_2D',ki,ks)
            print('STEP0_IMAGE_SHAPE=',str(im.shape),'SEG_SHAPE=',str(seg.shape))
            if im.shape!=seg.shape:
                print('STEP0_GRID=BLOCKED_SHAPE_MISMATCH')
                continue
            print('STEP0_GRID_1333=','VERIFIED_BYTES' if im.shape==(1333,1333) else 'DIFFERS_REPORTED_1333_NOT_DROPPED')
            print('STEP0_IMAGE_DTYPE=',im.dtype,'SEG_DTYPE=',seg.dtype)
            print('STEP0_BUNIT=',repr(hi.header.get('BUNIT','NOT_DECLARED')))
            print('STEP0_FILTER=',repr(hi.header.get('FILTER','NOT_DECLARED')))
            for key in ('MAGZERO','MAGZP','MAG_ZERO','MAGZERO_R','PHOTZP','ZP','ZEROPOINT','SKYSUB','SKY_SUB','BACKSUB','BKG_SUB','GAIN'):
                if key in hi.header:print('STEP0_IMAGE_HEADER_'+key+'=',repr(hi.header[key]))
            for key in ('PRODTYPE','SEGMAP','SEGTYPE','DETECT','DILATED','REGROUP','VERSION','ORIGIN'):
                if key in hs.header:print('STEP0_SEG_HEADER_'+key+'=',repr(hs.header[key]))
            for key in ('PC1_1','CD1_1','CDELT1','CDELT2'):
                if key in hi.header:print('STEP0_IMAGE_WCS_'+key+'=',repr(hi.header[key]))
                if key in hs.header:print('STEP0_SEG_WCS_'+key+'=',repr(hs.header[key]))
            try:
                from astropy.wcs import WCS
                from astropy.wcs.utils import proj_plane_pixel_scales
                if hi.header.get('CTYPE1') and hi.header.get('CTYPE2'):
                    pix=proj_plane_pixel_scales(WCS(hi.header).celestial)*3600.0
                    print('STEP0_PIXEL_SCALE_ARCSEC=',repr(pix.tolist()))
                    print('STEP0_SCALE_0p30_STATUS=', 'VERIFIED_HEADER' if np.all(np.isfinite(pix)) and np.allclose(pix,0.30,rtol=0,atol=0.002) else 'DIFFERS_FROM_REPORTED')
                else:print('STEP0_PIXEL_SCALE=UNRESOLVED_MISSING_CELESTIAL_WCS')
            except Exception as e:print('STEP0_PIXEL_SCALE=UNRESOLVED_WCS_ERROR',str(e)[:150])
            seg_numeric=np.issubdtype(seg.dtype,np.number)
            # Same shape alone is not proof of the same WCS pixel grid.
            try:
                from astropy.wcs import WCS
                if hi.header.get('CTYPE1') and hs.header.get('CTYPE1'):
                    iw=WCS(hi.header).celestial;sw=WCS(hs.header).celestial
                    rr=np.array([0,(im.shape[0]-1)/2,im.shape[0]-1],float)
                    cc=np.array([0,(im.shape[1]-1)/2,im.shape[1]-1],float)
                    xx,yy=np.meshgrid(cc,rr)
                    iworld=iw.pixel_to_world(xx.ravel(),yy.ravel())
                    sworld=sw.pixel_to_world(xx.ravel(),yy.ravel())
                    sep=iworld.separation(sworld).arcsec
                    print('STEP0_IMAGE_SEG_GRID_MAX_WCS_SEPARATION_ARCSEC=',float(np.max(sep)))
                    print('STEP0_IMAGE_SEG_WCS_ALIGNMENT=', 'VERIFIED_HEADER' if np.all(np.isfinite(sep)) and np.max(sep)<0.01 else 'BLOCKED_OR_DIFFERENT_WCS')
                else:print('STEP0_IMAGE_SEG_WCS_ALIGNMENT=UNRESOLVED_MISSING_WCS')
            except Exception as e:print('STEP0_IMAGE_SEG_WCS_ALIGNMENT=UNRESOLVED_ERROR',str(e)[:150])
            print('STEP0_SEG_INTEGER_LABELS=',bool(seg_numeric and np.all(np.isfinite(seg)) and np.all(seg==np.rint(seg))))
            if seg_numeric and np.all(np.isfinite(seg)):
                print('STEP0_SEG_NONZERO_LABELS=',int(np.unique(seg[seg!=0]).size))
                bg=(seg==0)&np.isfinite(im)
                nbg=int(np.count_nonzero(bg))
                print('STEP0_ORIGINAL_BACKGROUND_PIXEL_COUNT=',nbg)
                if nbg>=100:
                    # Diagnostic only; never treat as proof of upstream subtraction.
                    z=np.asarray(im[bg],float)
                    z=z[np.isfinite(z)]
                    median=float(np.median(z));mad=float(np.median(np.abs(z-median)))
                    print('STEP0_BACKGROUND_MEDIAN_DIAGNOSTIC=',repr(median))
                    print('STEP0_BACKGROUND_1p4826MAD_DIAGNOSTIC=',repr(1.4826*mad))
                    print('STEP0_SKY_SUBTRACTION=UNRESOLVED_PIXEL_MEDIAN_IS_NOT_PROOF_OF_PROCESSING_HISTORY')
                else:print('STEP0_SKY_BACKGROUND_DIAGNOSTIC=INSUFFICIENT_PIXELS')
            print('STEP0_SEGMENT_PRODUCT_SUBTYPE=REPORTED_DELIVERED_WAVES_MAP;UNRESOLVED_IF_NO_EXPLICIT_HEADER')
            print('STEP0_RADEC_LABEL_MAPPING=UNRESOLVED_REQUIRES_DOCUMENTED_CATALOGUE_COORDINATES_AND_SEG_WCS')
            print('STEP0_GAIN_687=REPORTED_PENDING_SOURCE_AND_REESTIMATION')
            print('STEP0_NO_PIXEL_SUPPORT_COMPARISON=TRUE')
    except Exception as e:
        print('STEP0_PAIR=BLOCKED_READ_EXCEPTION',type(e).__name__,str(e)[:240])
print('STEP0_G0_STATUS=EVIDENCE_ONLY_NOT_FIT_AUTHORISATION')
print('\nUNAVAILABLE_SOURCE_PATHS:', unavailable)
print('PRODUCTION_RUNNER_FOUND:', source_ok)
print('NEW_FITS=0 NEW_OPTIMISERS=0 SEGMENTATION_CHANGES=0')
print('EXACT_SEGMENT_SUPPORT_COMPARISON=NOT_RUN')
incomplete = unavailable or not rows or not source_ok
print('COLLECTION_STATUS=' + ('PARTIAL' if incomplete else 'EVIDENCE_COLLECTED'))
if incomplete:
    raise SystemExit(3)
BD_PY
BD_REMOTE
then
  awk '/^PATHS_TSV_BEGIN$/{printing=1;next} /^PATHS_TSV_END$/{printing=0} printing{print}' "$bd_audit_out/FRIDAY_STAGE0_EVIDENCE.txt" > "$bd_audit_out/paths.tsv"
  printf '\nEVIDENCE_LOG=%s/FRIDAY_STAGE0_EVIDENCE.txt\n' "$bd_audit_out"
  printf 'PATHS_TSV=%s/paths.tsv\n' "$bd_audit_out"
  printf '%s\n' 'Collection completed. This is metadata evidence, not a segment-coverage pass.'
else
  bd_audit_status=$?
  printf '\nEVIDENCE_LOG=%s/FRIDAY_STAGE0_EVIDENCE.txt\n' "$bd_audit_out" >&2
  printf '%s\n' 'Collection is incomplete. The partial log is preserved; do not rerun fitting.' >&2
  exit "$bd_audit_status"
fi
