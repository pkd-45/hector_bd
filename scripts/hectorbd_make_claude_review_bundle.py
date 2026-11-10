#!/usr/bin/env python3
"""Read-only, fail-closed Friday HectorBD pushed-commit + evidence snapshot.

Does NOT git add/commit/push. Run only AFTER the coding agent has reviewed,
committed and pushed all intended source changes. Stdlib only.
"""
import argparse
import datetime as dt
import glob
import hashlib
import io
import json
import os
import pathlib
import subprocess
import sys
import tarfile
import zipfile

ROOT = pathlib.Path.cwd().resolve()
BASE = '7fc5ebf3dfea6096a119e87cb6608c89921fb3bc'
MAX_EVIDENCE_BYTES = 8 * 1024 * 1024
MAX_SOURCE_BYTES = 20 * 1024 * 1024


def run(*args, check=True):
    p = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    if check and p.returncode:
        raise SystemExit(f"ERROR running {' '.join(args)}:\n{p.stderr}\n{p.stdout}")
    return p.stdout.strip()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def add_record(store, meta, target, data, source):
    if target in store:
        raise SystemExit(f'Duplicate archive entry: {target}')
    store[target] = data
    meta.append({'archive_path': target, 'source': source, 'sha256': sha(data), 'size_bytes': len(data)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--review', type=pathlib.Path, required=True,
                    help='The Claude scientific review Markdown supplied with this script')
    ap.add_argument('--plan', type=pathlib.Path,
                    help='Original CLAUDE_BD_RESTART_01_PLAN_AND_CHECKS.md (recommended)')
    ap.add_argument('--locations', type=pathlib.Path,
                    help='Original CLAUDE_BD_RESTART_02_FILE_LOCATIONS.md (recommended)')
    ap.add_argument('--output-dir', type=pathlib.Path, default=pathlib.Path('outputs/claude_review'))
    ap.add_argument('--extra-evidence', action='append', default=[], metavar='GLOB',
                    help='Additional repository-relative evidence glob; repeat as needed. '
                         'Unmatched patterns are reported as missing, never as passing evidence.')
    args = ap.parse_args()

    if not (ROOT / '.hector_bd_root').is_file():
        raise SystemExit('Run from ~/hector_bd repository root; marker not found')
    if run('git', 'rev-parse', '--show-toplevel') != str(ROOT):
        raise SystemExit('Wrong git root')
    if run('git', 'branch', '--show-current') != 'main':
        raise SystemExit('Publish intended changes to main first (do not auto-merge)')
    dirty = run('git', 'status', '--porcelain=v1', '-uall')
    if dirty:
        raise SystemExit('Uncommitted or untracked files remain. Review/stage intended sources, or explicitly ignore generated files:\n' + dirty)

    local = run('git', 'rev-parse', 'HEAD')
    remote_ls = run('git', 'ls-remote', 'origin', 'refs/heads/main')
    if not remote_ls:
        raise SystemExit('Could not resolve origin/main; no online-verification claim possible')
    remote = remote_ls.split()[0]
    if remote != local:
        raise SystemExit(f'NOT ONLINE: local HEAD {local} != GitHub origin/main {remote}. Push, then rerun.')
    if subprocess.run(['git','merge-base','--is-ancestor',BASE,local],cwd=ROOT).returncode:
        raise SystemExit('Review base not ancestor of current main; investigate history')

    review = args.review.expanduser().resolve()
    if not review.is_file():
        raise SystemExit('Review Markdown missing: '+str(review))
    store, meta, missing = {}, [], []
    add_record(store, meta, 'REVIEW/CLAUDE_REVIEW_REQUEST.md', review.read_bytes(), str(review))
    for option, dest in [(args.plan, 'REVIEW/CLAUDE_BD_RESTART_01_PLAN_AND_CHECKS.md'),
                         (args.locations, 'REVIEW/CLAUDE_BD_RESTART_02_FILE_LOCATIONS.md')]:
        if option:
            path = option.expanduser().resolve()
            if not path.is_file():
                raise SystemExit('Explicit reference missing: '+str(path))
            add_record(store, meta, dest, path.read_bytes(), str(path))
        else:
            missing.append('OPTIONAL_REFERENCE_NOT_SUPPLIED: '+dest)

    # Source tree is EXACTLY the pushed git HEAD, not mutable working files.
    src_bytes = subprocess.run(['git','archive','--format=tar','HEAD'],cwd=ROOT,
                               capture_output=True,check=True).stdout
    with tarfile.open(fileobj=io.BytesIO(src_bytes), mode='r:') as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            if member.size > MAX_SOURCE_BYTES:
                missing.append(f'TRACKED_SOURCE_TOO_LARGE: {member.name} ({member.size})')
                continue
            data = tf.extractfile(member).read()
            add_record(store, meta, 'PUSHED_SOURCE/' + member.name, data,
                       f'git archive {local}:{member.name}')

    # Small existing evidence only. This does not generate or re-run any science.
    evidences = [
        'outputs/r_environment_1443187/package_manifest.tsv',
        'outputs/r_environment_1443677.out',
        'outputs/r_environment_1443677.err',
        'outputs/*1443186*.out', 'outputs/*1443186*.err',
        'outputs/*1453593*.out', 'outputs/*1453593*.err',
        'outputs/step2_calibration_1437317/calibration.tsv',
        'outputs/step2_calibration_1437317/calibration.json',
        'outputs/step2_calibration_1437317/SHA256SUMS.txt',
        'outputs/step1_catalogue_1427444/summary.json',
        'outputs/step1_catalogue_1427444/step1_segmap_audit.tsv',
        'outputs/*spatial*/*.json', 'outputs/*spatial*/*.tsv',
        'outputs/*pilot*/*.json', 'outputs/*pilot*/*.tsv',
        'outputs/*historical*/*.tsv', 'outputs/*historical*/*.json',
    ]
    for pattern in args.extra_evidence:
        if pathlib.Path(pattern).is_absolute() or '..' in pathlib.Path(pattern).parts:
            raise SystemExit('Extra evidence patterns must stay below the repository: '+pattern)
    evidences.extend(args.extra_evidence)
    covered = set()
    for pattern in evidences:
        matches = sorted(glob.glob(str(ROOT / pattern)))
        if not matches:
            missing.append('NOT_FOUND: '+pattern)
        for p in matches:
            path = pathlib.Path(p)
            if path in covered or not path.is_file():
                continue
            covered.add(path)
            rel = path.relative_to(ROOT).as_posix()
            if path.stat().st_size > MAX_EVIDENCE_BYTES:
                missing.append('EVIDENCE_TOO_LARGE_NOT_INCLUDED: '+rel)
                continue
            add_record(store, meta, 'FRIDAY_EVIDENCE/' + rel, path.read_bytes(), str(path))

    diff = run('git', 'diff', '--stat', f'{BASE}..{local}')
    names = run('git', 'diff', '--name-status', f'{BASE}..{local}')
    logs = run('git', 'log', '-10', '--format=%h %ad %s', '--date=iso-strict')
    prov = '\n'.join([
        '# Pushed GitHub snapshot — machine-checked on Friday',
        '',f'- Created (UTC): {dt.datetime.now(dt.timezone.utc).isoformat()}',
        f'- GitHub repository origin: `{run("git","remote","get-url","origin")}`',
        f'- Friday HEAD SHA: `{local}`', f'- GitHub origin/main SHA: `{remote}`',
        f'- Exact match: **YES**', '- Git working tree: **clean**',
        f'- Review diff base SHA: `{BASE}`', '', '## Changed paths since review base',
        '```', names, '```', '', '## Diff statistics', '```', diff, '```',
        '', '## Recent commits', '```', logs, '```',
        '', '## Missing/omitted evidence',
        'All entries below are explicitly **unverified** in this bundle. Check Friday source if necessary.',
        '```', '\n'.join(missing) if missing else 'None', '```',
        '', '## Dataset safety',
        'No FITS images, RDS fits, raw WAVES tiles or source data were included.',
        'SLURM results are archived only if pre-existing and below size limits.',
        '', '## Interpretation rule',
        'Archive presence does not by itself establish a passing test. Read and audit the logs.',
        'Missing result => NOT INSPECTED, not PASS.', ''
    ]).encode()
    add_record(store, meta, 'REVIEW/PUSH_PROVENANCE.md', prov, 'generated by review bundle helper')
    manifest = {'utc_created':dt.datetime.now(dt.timezone.utc).isoformat(),
                'local_head':local,'remote_main':remote,'base':BASE,'worktree_clean':True,
                'missing_or_omitted':missing,'files':meta}
    store['MANIFEST_SHA256.json'] = (json.dumps(manifest,indent=2,sort_keys=True)+'\n').encode()

    out = args.output_dir.expanduser().resolve()
    out.mkdir(parents=True,exist_ok=True)
    target = out / ('HECTORBD_CLAUDE_REVIEW_' + local[:12] + '.zip')
    if target.exists():
        raise SystemExit('Archive already exists; refusing overwrite: '+str(target))
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as zf:
        for name in sorted(store):
            zf.writestr(name,store[name])
    with zipfile.ZipFile(target) as zf:
        bad = zf.testzip()
        if bad:
            raise SystemExit('ZIP CRC failure: '+bad)
    digest = sha(target.read_bytes())
    check = target.with_suffix('.zip.sha256')
    check.write_text(f'{digest}  {target.name}\n')
    print('PUSH_MATCH_VERIFIED',local)
    print('BUNDLE',target)
    print('SHA256',digest)
    print('ENTRIES',len(store))
    print('MISSING_OR_OMITTED',len(missing))
    print('REVIEW_BASE',BASE)
    print('NO_OBSERVED_FITTING_RUN')

if __name__ == '__main__':
    main()
