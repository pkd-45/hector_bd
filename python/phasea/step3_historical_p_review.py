#!/usr/bin/env python3
"""Read-only metadata review of the 68 P endpoints selected by frozen R61.

The selection binds through R61's exact task index and R53's exact run plan.
No model is rendered, no fit is performed, and absent settings remain absent.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

R61 = Path('/scratch/uqpdas/WAVES_BD_PROFIT_PROFUSE_R61_R53_NATIVE_SPD_FIXED4_POLICY_EXECUTION_20260930')
R53 = Path('/scratch/uqpdas/WAVES_BD_PROFIT_PROFUSE_R53R1_FULL832_PRODUCTION_EXECUTION_PATH_REPAIR_20260921')
FROZEN_SHA = '874a350440b0361371db15b4f4e363bfc79eb924a31ff0ee5ae478ca41d275f3'
METRICS_SHA = 'f71a7f1b8fe23d094c5e247591983ca3f40d2c8cf5d1522b6cd9a3ec2ae6ccd5'
INDEX_SHA = '11e36004738e03b1f490ee353ebafc4dbb22a2bf49678214a64ae40965524a5e'
PLAN_SHA = 'dc94fd640918fa1481789d4ea388f8718b2f2eb38bc6aa84f2315b8f7905d6ff'
POLICY = 'R53_NATIVE_S_P_D_FIXED4_MODEL_SELECTION_V1'
EXPECTED_COUNTS = {'S': 352, 'P': 68, 'D_FIXED4': 412}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def git_provenance(repo, files):
    """Bind executed sources to the committed tree (Claude review P1-3)."""
    def git(*args):
        return subprocess.run(['git', '-C', str(repo), *args], check=True, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.strip()
    head = git('rev-parse', 'HEAD')
    sources = {}
    for path in files:
        relative = Path(path).resolve().relative_to(repo).as_posix()
        committed = git('rev-parse', f'HEAD:{relative}')
        executed = git('hash-object', '--', relative)
        sources[relative] = {'committed_blob': committed, 'executed_blob': executed,
                             'matches_HEAD': committed == executed}
    return {'git_HEAD': head, 'git_tracked_changes': git('status', '--porcelain',
            '--untracked-files=no').splitlines(), 'executed_sources': sources}


def centroid_offset_range(records):
    """Range of accepted-kernel physical centroid offsets from the grid centre (x, y)."""
    offsets = [r['accepted_psf_geometry']['centroid_offset_from_grid_centre_xy']
               for r in records if r.get('accepted_psf_geometry')]
    if not offsets:
        return None
    return {'count': len(offsets),
            'x_min': min(o[0] for o in offsets), 'x_max': max(o[0] for o in offsets),
            'y_min': min(o[1] for o in offsets), 'y_max': max(o[1] for o in offsets)}


def read_tsv(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def unique_index(rows, key):
    result = {}
    for row in rows:
        value = row[key]
        if value in result:
            raise ValueError(f'Duplicate {key}: {value}')
        result[value] = row
    return result


def checked_table(path, expected_sha, count):
    actual = sha256(path)
    if actual != expected_sha:
        raise ValueError(f'Frozen source SHA mismatch: {path}: {actual}')
    rows = read_tsv(path)
    if len(rows) != count:
        raise ValueError(f'Frozen table row count mismatch: {path}: {len(rows)}')
    return rows


def bindings(r61, r53):
    selection_path = r61 / 'results/R61_NATIVE_MODEL_SELECTION_832.tsv'
    metrics_path = r61 / 'results/R61_ENDPOINT_COMMON_METRICS_7488.tsv'
    index_path = r61 / 'contract/R61_ENDPOINT_TASK_INDEX_7488.tsv'
    plan_path = r53 / 'contract/R52_FULL832_PLAN_7488.tsv'
    roster = checked_table(selection_path, FROZEN_SHA, 832)
    unique_index(roster, 'uberID')
    if Counter(row['global_candidate'] for row in roster) != EXPECTED_COUNTS:
        raise ValueError('Frozen model counts differ from S352/P68/D_FIXED4=412')
    if any(row['native_policy'] != POLICY for row in roster):
        raise ValueError('Unexpected selection policy')
    metrics = unique_index(checked_table(metrics_path, METRICS_SHA, 7488), 'run_index')
    index = unique_index(checked_table(index_path, INDEX_SHA, 7488), 'run_index')
    plan = unique_index(checked_table(plan_path, PLAN_SHA, 7488), 'run_index')
    selected = []
    for row in roster:
        if row['global_candidate'] != 'P':
            continue
        uid = row['uberID']
        if not re.fullmatch(r'[0-9]+', uid):
            raise ValueError('Non-numeric frozen uberID')
        run = row['global_candidate_run_index']
        start = row['global_candidate_start_id']
        task, metric, source = index[run], metrics[run], plan[run]
        for field in ('run_index', 'uberID', 'model', 'start_id'):
            if len({task[field], metric[field], source[field]}) != 1:
                raise ValueError(f'Binding identity mismatch {uid} {field}')
        if (task['uberID'], task['model'], task['start_id']) != (uid, 'P', start):
            raise ValueError(f'Selected endpoint identity mismatch {uid}')
        if row['P_within_model_best_run_index'] != run or row['P_within_model_best_start'] != start:
            raise ValueError(f'Within-model endpoint mismatch {uid}')
        if task['accepted_fit_rds'] != metric['accepted_fit_rds']:
            raise ValueError(f'R61 task/metric path mismatch {uid}')
        accepted = Path(task['accepted_fit_rds'])
        # The path is checked against the recorded task, never guessed from family.
        expected_path = r53 / 'runs' / uid / 'P' / start / 'accepted_fit.rds'
        if accepted != expected_path or not accepted.is_file():
            raise ValueError(f'Unexpected/missing accepted endpoint {uid}: {accepted}')
        source_path = Path(source['source_rds'])
        if not source_path.is_absolute() or not source_path.is_file():
            raise ValueError(f'Missing original source RDS {uid}: {source_path}')
        result_path = accepted.with_name('result.tsv')
        result = read_tsv(result_path)
        if len(result) != 1:
            raise ValueError(f'Unexpected R53 result rows {uid}')
        for field in ('run_index', 'uberID', 'model', 'start_id', 'source_rds'):
            if result[0][field] != source[field]:
                raise ValueError(f'R53 result/plan mismatch {uid} {field}')
        run_log = source_path.with_name('run.log')
        log_lines = run_log.read_text(errors='replace').splitlines() if run_log.is_file() else []
        references = sorted(set(match.group(1).strip() for line in log_lines
            if (match := re.search(r'reusing cached PSF:\s*(/[^\r\n]+)', line))))
        cache = references[0] if len(references) == 1 else ''
        selected.append({
            'uberID': uid, 'model': 'P', 'selection_layer': 'R61', 'fit_layer': 'R53',
            'run_index': run, 'start_id': start,
            'accepted_fit_rds': str(accepted), 'source_rds': str(source_path),
            'r53_result_tsv': str(result_path), 'r53_result_sha256': sha256(result_path),
            'source_run_log': str(run_log),
            'source_run_log_sha256': sha256(run_log) if run_log.is_file() else '',
            'psf_recorded_cache_path': cache,
            'psf_cache_path_status': ('ONE_EXPLICIT_LOG_REFERENCE' if len(references) == 1
                else 'ABSENT_LOG_REFERENCE' if not references else 'AMBIGUOUS_LOG_REFERENCES'),
            'psf_log_lines_json': json.dumps([
                {'line': n, 'text': text} for n, text in enumerate(log_lines, 1)
                if re.search(r'psf|finesample', text, re.I)], separators=(',', ':')),
            'r61_delta_BIC': row['delta_BIC'],
            'r61_native_selection_class': row['native_selection_class'],
        })
    if len(selected) != 68:
        raise ValueError('Selected P endpoint count is not 68')
    return selected, {'roster_sha256': FROZEN_SHA, 'metrics_sha256': METRICS_SHA,
        'task_index_sha256': INDEX_SHA, 'r53_plan_sha256': PLAN_SHA,
        'frozen_model_counts': EXPECTED_COUNTS, 'selected_p_count': 68,
        'selection_policy': POLICY, 'r74_selection_fields_used': False}


def summarise(evidence, records):
    """Counts and per-record kernel evidence written to summary.json."""
    summary = dict(evidence)
    summary.update({'record_count': len(records),
        'accepted_finesample_presence_counts': dict(Counter(r['accepted']['Data_finesample']['presence'] for r in records)),
        'accepted_psftype_presence_counts': dict(Counter(r['accepted']['Data_psftype']['presence'] for r in records)),
        'accepted_finesample_values': dict(Counter(json.dumps(r['accepted']['Data_finesample']['value'], sort_keys=True) for r in records)),
        'accepted_psftype_values': dict(Counter(json.dumps(r['accepted']['Data_psftype']['value'], sort_keys=True) for r in records)),
        'metadata_flags': dict(Counter(flag for r in records for flag in r['metadata_flags'])),
        'accepted_source_kernel_comparison': dict(Counter(r['accepted_source_psf_comparison'] for r in records)),
        'accepted_source_kernel_value_classification': dict(Counter(
            r['accepted_source_psf_value_comparison']['class'] for r in records)),
        'accepted_cache_kernel_value_classification': dict(Counter(
            f'{key}:{value["class"]}' for r in records
            for key, value in (r['accepted_cache_kernel_value_comparisons'] or {}).items())),
        'accepted_source_kernel_differences': [{
            'uberID': r['uberID'], 'hash_comparison': r['accepted_source_psf_comparison'],
            **{key: r['accepted_source_psf_value_comparison'].get(key) for key in (
                'class', 'same_dimensions', 'sum_a', 'sum_b', 'max_abs_raw_difference_common',
                'max_abs_normalised_difference', 'normalised_difference_relative_to_peak',
                'pixels_differing_common', 'pixels_common')}}
            for r in records if r['accepted_source_psf_comparison'] != 'EXACT_NUMERIC_KERNEL_MATCH'],
        'accepted_kernel_centroid_offset_from_grid_centre_px': centroid_offset_range(records),
        'fit_authorized': False, 'scientific_interpretation': 'Historical metadata review only; absent settings do not establish a rendering defect.'})
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--r61', type=Path, default=R61)
    ap.add_argument('--r53', type=Path, default=R53)
    ap.add_argument('--rscript', default='/home/uqpdas/.conda/envs/waves_r/bin/Rscript')
    args = ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        ap.error('Run through SLURM; bulk historical RDS reads are forbidden on login')
    repo = args.repo.resolve()
    if not (repo / '.hector_bd_root').is_file():
        ap.error('Missing repository marker')
    out = args.out.absolute()
    if out.exists() or out.is_symlink():
        ap.error('Output must be a new, nonexistent directory')
    if (repo / 'outputs').resolve() not in out.resolve().parents:
        ap.error('Output must be under this repository outputs directory')
    helper = repo / 'python/phasea/src/export_historical_p_metadata.R'
    provenance = git_provenance(repo, [Path(__file__), helper])
    if not all(s['matches_HEAD'] for s in provenance['executed_sources'].values()):
        ap.error('Driver or R helper differs from the committed HEAD; commit before running')
    selected, evidence = bindings(args.r61, args.r53)
    out.mkdir(parents=True, exist_ok=False)
    binding_path = out / 'selected_p_bindings.tsv'
    with binding_path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(selected[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(selected)
    evidence.update({'job_id': os.environ['SLURM_JOB_ID'], 'script_sha256': sha256(__file__),
        'helper_sha256': sha256(helper), 'bindings_sha256': sha256(binding_path),
        'source_provenance': provenance,
        'r61_root': str(args.r61), 'r53_root': str(args.r53),
        'scope': 'READ_ONLY_HISTORICAL_METADATA_NO_RENDERING_NO_FITTING',
        'absent_settings_policy': 'ABSENT_IS_UNKNOWN_NEVER_DEFAULTED_OR_INFERRED'})
    (out / 'binding_evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    result = subprocess.run([args.rscript, '--vanilla', str(helper),
        str(binding_path), str(out)], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)
    (out / 'r_metadata_export.log').write_text(result.stdout)
    print(result.stdout, end='', flush=True)
    if result.returncode:
        raise SystemExit(result.returncode)
    records = [json.loads(line) for line in (out / 'historical_p_metadata.jsonl').read_text().splitlines()]
    if len(records) != 68 or {r['uberID'] for r in records} != {r['uberID'] for r in selected}:
        raise ValueError('R metadata output membership/count mismatch')
    summary = summarise(evidence, records)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)
    print('HISTORICAL_SELECTED_P_METADATA_REVIEW_COMPLETE', flush=True)


if __name__ == '__main__':
    main()
