#!/usr/bin/env python3
"""Run all finite diagnostics and record ordered logs and JSON reports.

Usage (from any directory, Python 3.10+, standard library only):
    python3 -I -B run_checks.py [--out DIR]

Outputs go to DIR (default: <bundle>/results), so the suite also runs from a
read-only copy of the bundle. The same interpreter flags (-I, -B, -O) are
passed to every check script. All checks are explicit (no `assert`), so the
suite is equally valid under `python -O`. Exit status is nonzero if any
script fails. The mathematical proofs are in manuscript.md, not in this runner.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
TASKS = [
    ('checks/exact_rational_checks.py', 'exact_rational_log.txt', 'exact_rational_checks.json'),
    ('checks/second_exact_checks.py', 'second_exact_log.txt', 'second_exact_checks.json'),
    ('checks/check_convex_polygons.py', 'convex_polygon_log.txt', 'convex_polygon_checks.json'),
    ('checks/extension_checks.py', 'extension_log.txt', 'extension_checks.json'),
]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default=str(ROOT / 'results'),
                        help='output directory for logs and JSON (default: <bundle>/results)')
    out = Path(parser.parse_args().out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    flags = (['-I'] if sys.flags.isolated else []) + (['-B'] if sys.dont_write_bytecode else []) \
        + ['-O'] * sys.flags.optimize
    records = []
    for source, log, report in TASKS:
        started = time.perf_counter()
        wall_start = time.time()
        result = subprocess.run([sys.executable, *flags, str(ROOT / source), '--out', str(out)],
                                cwd=str(out), text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        (out / log).write_text(result.stdout)
        summary = None
        if (out / report).is_file() and (out / report).stat().st_mtime >= wall_start - 1:  # ignore stale reports
            data = json.loads((out / report).read_text())
            summary = {k: data.get(k) for k in ('status', 'checks_evaluated', 'checks_failed')}
        ok = result.returncode == 0 and summary is not None and summary['status'] == 'PASS'
        records.append({'script': source, 'script_sha256': sha256(ROOT / source),
                        'exit_code': result.returncode, 'passed': ok, 'log': log,
                        'report': report, 'report_summary': summary,
                        'elapsed_seconds': round(time.perf_counter() - started, 2)})
        print(f"{source}: {'PASS' if ok else 'FAIL'} (exit {result.returncode}"
              + (f", {summary['checks_evaluated']} checks, {summary['checks_failed']} failed" if summary else '')
              + ')', flush=True)
        if not ok:
            print(result.stdout[-4000:])
    run = {
        'date': datetime.date.today().isoformat(),
        'python': platform.python_version(),
        'implementation': platform.python_implementation(),
        'platform': platform.platform(),
        'interpreter_flags': flags,
        'scope': 'Finite analytic-identity and numerical diagnostics; no formal verification.',
        'manuscript_sha256': sha256(ROOT / 'manuscript.md'),
        'checks': records,
        'total_checks_evaluated': sum((r['report_summary'] or {}).get('checks_evaluated') or 0 for r in records),
        'all_passed': all(r['passed'] for r in records),
    }
    (out / 'run_summary.json').write_text(json.dumps(run, indent=2) + '\n')
    print(f"all passed: {run['all_passed']}; {run['total_checks_evaluated']} explicit checks; outputs in {out}")
    if not run['all_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
