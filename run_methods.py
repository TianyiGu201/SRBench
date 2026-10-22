"""Run independent method scripts sequentially from an explicit command list."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid


def load_jobs(filename, methods=None):
    filename = Path(filename).resolve()
    config = json.loads(filename.read_text(encoding='utf-8'))
    root = (filename.parent / config.get('project_root', '..')).resolve()
    context = {'root': str(root), 'package': str(Path(__file__).resolve().parent)}

    def expand(value):
        return str(value).format_map(context)

    def path(value):
        return str((root / expand(value)).resolve())

    jobs, identifiers = [], set()
    for item in config.get('jobs', []):
        identifier = item['id']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', identifier) or identifier in identifiers:
            raise ValueError(f'Invalid or duplicate job id: {identifier}')
        identifiers.add(identifier)
        if methods and item['method'] not in methods:
            continue
        steps = item.get('steps', [item.get('command')])
        if not steps or any(not isinstance(step, list) or not step or not all(isinstance(arg, str) for arg in step) for step in steps):
            raise ValueError(f'{identifier}: supply a command argument list or a steps list of argument lists')
        jobs.append({'id': identifier, 'method': item['method'],
                     'cwd': path(item.get('cwd', '.')),
                     'steps': [[expand(arg) for arg in step] for step in steps],
                     'inputs': [path(value) for value in item.get('inputs', [])],
                     'outputs': [path(value) for value in item.get('outputs', [])],
                     'env': {key: expand(value) for key, value in item.get('env', {}).items()}})
    if not jobs:
        raise ValueError('No matching method jobs')
    if methods and set(methods) - {item['method'] for item in jobs}:
        raise ValueError('Requested methods are absent from the configuration')
    return jobs, Path(path(config.get('log_dir', 'RUN_LOGS')))


def check_job(job):
    if not Path(job['cwd']).is_dir():
        raise FileNotFoundError(f'Working directory not found: {job["cwd"]}')
    for value in job['inputs']:
        if not Path(value).exists():
            raise FileNotFoundError(f'Input not found: {value}')
    for value in job['outputs']:
        if Path(value).exists():
            raise FileExistsError(f'Output already exists; choose a fresh run directory: {value}')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--methods', nargs='+', help='Run only these method names')
    parser.add_argument('--dry-run', action='store_true', help='Print commands without running or checking external inputs')
    parser.add_argument('--check', action='store_true', help='Check working directories and declared input/output paths only')
    parser.add_argument('--keep-going', action='store_true', help='Continue independent jobs after a failure')
    args = parser.parse_args(argv)
    jobs, log_root = load_jobs(args.config, args.methods)
    if args.dry_run:
        print(json.dumps(jobs, ensure_ascii=False, indent=2))
        return 0
    if args.check:
        for job in jobs:
            check_job(job)
            print(f'Checked {job["id"]}')
        return 0

    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid.uuid4().hex[:8]
    log_dir = log_root / stamp
    log_dir.mkdir(parents=True, exist_ok=False)
    results = []
    for job in jobs:
        start = time.monotonic()
        result = {'id': job['id'], 'method': job['method'], 'cwd': job['cwd'], 'steps': []}
        try:
            check_job(job)
            print(f'Running {job["id"]}', flush=True)
            for number, command in enumerate(job['steps'], 1):
                log_file = log_dir / f'{job["id"]}_step{number}.log'
                with log_file.open('w', encoding='utf-8') as log:
                    log.write(json.dumps(command, ensure_ascii=False) + '\n')
                    log.flush()
                    completed = subprocess.run(command, cwd=job['cwd'],
                        env={**os.environ, 'PYTHONUTF8': '1', **job['env']},
                        stdout=log, stderr=subprocess.STDOUT, shell=False)
                result['steps'].append({'command': command, 'returncode': completed.returncode, 'log': log_file.name})
                if completed.returncode:
                    raise RuntimeError(f'Step {number} exited with code {completed.returncode}')
            missing = [value for value in job['outputs'] if not Path(value).exists()]
            if missing:
                raise FileNotFoundError(f'Command completed but declared outputs are missing: {missing}')
            result['status'] = 'completed'
        except Exception as error:
            result.update(status='failed', error=str(error))
            print(f'Failed {job["id"]}: {error}', file=sys.stderr, flush=True)
        result['elapsed_seconds'] = time.monotonic() - start
        results.append(result)
        (log_dir / 'summary.json').write_text(json.dumps(results, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        if result['status'] == 'failed' and not args.keep_going:
            break
    print(f'Run records: {log_dir}')
    return 1 if any(result['status'] == 'failed' for result in results) else 0


if __name__ == '__main__':
    raise SystemExit(main())
