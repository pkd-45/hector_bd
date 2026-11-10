# Friday source and scientific evidence handoff — 2026-11-10

Repository: `pkd-45/hector_bd`, branch `main`. The independently checked
starting commit was `7fc5ebf3dfea6096a119e87cb6608c89921fb3bc`.
The next step is Claude's independent review of the pushed source and actual
numerical evidence. This handoff does not authorize observed fitting.

## Intended source changes

- `python/phasea/step2_calibration.py` and its tests/SLURM wrapper: read-only
  three-control delivered-image calibration diagnostics.
- `python/phasea/step2_spatial_calibration.py` and its tests/SLURM wrapper:
  original-label sky-mask dilation comparisons, spatial cells, held-out plane
  diagnostics, raw/detrended covariance, blank apertures and block uncertainty.
  Insufficient coverage remains explicit; no universal effective pixel count
  or source-Poisson gain is certified.
- `python/phasea/step3_historical_p_review.py`,
  `python/phasea/src/export_historical_p_metadata.R` and the SLURM wrapper:
  metadata extraction for the 68 P endpoints from the exact R61/R53 bindings.
  Missing settings remain unknown; no historical fit is rerun.
- `scripts/freeze_pilot.py`, its tests and SLURM wrapper: a review-only
  candidate from the exact frozen 832-object tables and existing audit.
  It proposes replacing the WCS-ineligible anchor, records the proposal and
  never changes the frozen roster or approval configuration.
- `scripts/prepare_pilot_inputs.py`, `scripts/verify_pilot_fits.R`, their tests
  and `slurm/pilot_preparation.sbatch`: byte-verified decompressed copies and
  Astropy/Rfits orientation checks for W11690092907692 in a fresh output
  directory. The path index remains a blocked candidate.
- `scripts/build_pilot_psf_candidate.R`, its SLURM wrapper,
  `tests/test_psf_centering.R` and `slurm/psf_geometry.sbatch`: candidate PSF
  preparation source and a focused synthetic geometry diagnostic. Their
  existence is not E1–E7/G5 acceptance or authorization to run star fits.
- Synthetic recovery and notebook/CLI tests gain diagnostic evidence output;
  `tests/run_synthetic_validation.R` and its SLURM wrapper preserve numerical
  recovery, optimizer and runtime evidence without changing tolerances.
- `reference/calibration_evidence_gap.md` documents missing image-production
  provenance and the information needed to test calibration/noise assumptions.
- The supplied scientific review request and original Claude plan/location
  references are preserved under `reference/reviews/`. A copy of the supplied
  `hectorbd_make_claude_review_bundle.py` is tracked under `scripts/`, with an
  optional `--extra-evidence` glob argument to include detailed existing
  numerical evidence and report missing final outputs. The supplied original
  at `/home/uqpdas/hectorbd_make_claude_review_bundle.py` remains unchanged.

The audit included tracked diffs, every new source file and ignored paths.
`outputs/`, `data/`, `env/restart_Rlib/`, FITS, RDS, archives and logs stay
excluded by `.gitignore`; `config/verified_path_index.tsv` is also ignored.
No active global exclude file or additional local exclude rules were found.
Ignored output copies of test scripts duplicate the intended source; the older
`outputs/build_step1_claude_package.py` is historical Step 1 packaging code.
Generated libraries, raw science products and large logs are not source changes.

## Existing results collected without repeating science calculations

| Job | Result and evidence location |
| --- | --- |
| 1437317 | Completed coarse three-control calibration diagnostics: `outputs/step2_calibration_1437317/`. These do not verify calibration. |
| 1443677 | Completed pinned R environment/renderer run, including F1–F9 PASS: `outputs/r_environment_1443677.out`. These are synthetic renderer results, not per-object PSF acceptance. |
| 1453593 | Completed synthetic recovery and notebook/CLI parity; `RECOVERY_EXIT=0`, `PARITY_EXIT=0`: `outputs/synthetic_validation_1453593/`. |
| 1453639 | Completed all 10 spatial synthetic tests; stderr contains the unittest results. |
| 1453640 | Completed three-control spatial diagnostics: `outputs/step2_spatial_calibration_1453640/`. Units, gain and noise provenance remain unverified. |
| 1453616 | Failed after writing 68 historical-P metadata records: `outputs/step3_historical_p_1453616/`. No successful original `summary.json` exists. |

Synthetic recovery records all eight physical recovery checks passing.
The D_FIXED4 objective is approximately -1249.378 and BD_FREE_N approximately
-1249.424. Both `inputs.rds` and `g2.rds` have matching notebook/CLI byte hashes.
The three test-source hashes match the current source. Spatial script and
wrapper hashes likewise match their recorded evidence. Historical results
record the old base commit plus source hashes; a later source commit does not
turn them into executions of newly tracked files.

The 68 historical JSONL records have unique IDs matching all selected bindings.
They record accepted `Data$finesample=1` and `Data$psftype="empirical"` in all
68 records. Accepted/source kernels match numerically for 55 records and differ
for 13. These differences alone do not establish a rendering defect.
Accepted `Data$psffinesampled`, saved ProFit SHA and saved R version are absent
in all 68; source `Data$psffinesampled` and saved ProFit SHA are also absent.
Two records have no unique explicit PSF-cache log reference.

The failed R helper logged extraction completion, then a trailing parse error
(`unexpected ','` at `lags=I(unique(flags)),`). The current helper parses, but
has not been re-executed against the historical RDS files. Its executed SHA was
`733111386808caf30634bc694b7fa0b2daf0b556dd228182874ad3ac2dc9be14`;
the current helper SHA is
`5297a3c7925e19965646b2354bf3a2545b50d9fc519a015c5490773b900a5a93`.
The Python extraction driver matches its executed SHA. Preserve this
distinction; the original failed job is not a successful review completion.

## Already-submitted preparation jobs

Jobs 1453626 (synthetic PSF geometry), 1453629 (review-only pilot candidate)
and 1453634 (input unit tests) originally depended on successful historical
metadata job 1453616. They do not consume that job's results. Only that
scheduling dependency was removed to complete the already-submitted work.
No new job was submitted, no assertions were changed, and the original failure
remains recorded. Job 1453635 (input copies) retained its dependencies on
1453634 and 1453629. All four jobs reached terminal states:

| Job | Collected result |
| --- | --- |
| 1453626 | **FAILED**, exit 1:0. All six direct physical-centre controls pass; 15 of 21 total cases fail the unchanged centring/symmetry assertions. Current v14 cases show roughly 0.42–0.50 pixel offsets. Numerical metrics and runtime are in `outputs/psf_geometry_1453626/centering_metrics.tsv` and `centering_report.json`. No PSF source correction or acceptance change was made. |
| 1453629 | **COMPLETED**, exit 0:0. Seven selection tests pass. The candidate uses 799 eligible objects and contains S=7, P=6, D_FIXED4=7; `outputs/pilot_roster_1453629/pilot20.json`. Candidate IDs SHA256: `6319da64cf16215f140c2930ef837d502016bf99480053b97ff87eba8750cdc7`. The frozen 832 roster and approved configuration remain unchanged. |
| 1453634 | **COMPLETED**, exit 0:0. All four input-copy unit tests pass; `outputs/pilot_preparation_1453634.err`. |
| 1453635 | **FAILED**, exit 1:0. Gzip copies, `copy_manifest.json` and `Rfits_verification.json` were written to `outputs/pilot_inputs_1453635/`. The Python driver stopped at `image: explicit off-diagonal Rfits sample differs`, after checking full-array numeric hash/dimensions/finite counts. This is not successful input verification. The final originals-after record, candidate path index, candidate Step 0 record and completion checksums were not produced. No PSF candidate or observed fitting job was submitted. |

Independent comparison of the two existing JSON reports shows matching full
numeric array hashes, dimensions and finite counts for both image and labels;
the image's serialized sample values differ in their final decimal digits.
This observation does not change the original job's failure or its assertions.

## Verification and boundaries

The three coarse-calibration tests, seven pilot-selection tests and four
input-copy/orientation unit tests pass. The actual input-copy job still failed
as recorded above; unit-test success does not certify those real products.
All changed/new Python files pass AST parsing; all eight new SLURM files pass
`bash -n`. All six changed/new R files parse in the existing Friday R runtime.
No observed fitting, repeated extraction or repeated spatial/recovery run was
launched during this handoff.

`config/step0.yaml`, `config/run_pilot.yaml`, `config/psf_g4_pilot.tsv`, original
segmentation and the frozen 832-object roster are unchanged. Counts remain
S=352, P=68 and D_FIXED4=412, bound to roster SHA256
`874a350440b0361371db15b4f4e363bfc79eb924a31ff0ee5ae478ca41d275f3`.
Units/zeropoint, upstream sky, effective gain 687 and the noise protocol remain
unverified. Candidate preparation, synthetic success and historical metadata
do not promote Step 0 or per-object PSF/scientific acceptance gates.

## Review input identity and packaging

The user transferred the original review files during this handoff. Each
transferred file matches its member of
`outputs/review_inputs/HECTORBD_CLAUDE_REVIEW_STARTER_20261108.zip` (CRC checked):

| Supplied file | SHA256 |
| --- | --- |
| `HECTORBD_CLAUDE_SCIENTIFIC_REVIEW_20261108.md` | `c4cf7b128c1b5a276a2fa4637c0a13ae5a33349d1ecbe16639c8960d196dd205` |
| `CLAUDE_BD_RESTART_01_PLAN_AND_CHECKS.md` | `c00a70ae7d0808e63d7402198e1ded155313565ad463ac8dce7e58db09d9be87` |
| `CLAUDE_BD_RESTART_02_FILE_LOCATIONS.md` | `c6bba4d8eb7492a5dc750aaeebe00a9bd876551b99e10e1d8e770d839298de50` |
| Original packaging helper | `df513632f622d36a277e9b47557741d9ce61bfc5c9e129a583508c16c7d2b8b5` |

The supplied review Markdown is preserved byte-for-byte under ignored
`outputs/review_inputs/`; only four trailing Markdown hard-break spaces are
normalized in the tracked copy so `git diff --cached --check` passes. The
bundle's `REVIEW/CLAUDE_REVIEW_REQUEST.md` uses the original bytes. Its SHA placeholders
are resolved by the generated `REVIEW/PUSH_PROVENANCE.md` and
`MANIFEST_SHA256.json` after pushing; no self-referential commit SHA is inserted
into a tracked source file. The starter ZIP is retained under ignored outputs,
not committed. The root source checksum manifest is refreshed for this source
snapshot; the historical V03.1 manifest is preserved separately.

Packaging uses `git archive HEAD` only after a clean source tree and an exact
independent `git ls-remote origin refs/heads/main` match. Additional evidence
patterns include detailed synthetic recovery/parity records, historical JSONL,
spatial cells/covariance, WCS diagnostics, PSF geometry and preparation job logs.
Missing historical summary and input-copy completion products must remain
explicitly missing. No FITS/RDS, generated libraries or diagnostic sky-pixel
arrays are included. Failed assertions remain failed in the archived evidence.
