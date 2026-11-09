"""Run driver.py under ThreadSanitizer, one process per case, and count reports.

Usage: python tests/tsan/run.py BUILD_DIR OUT_DIR [options] [CASE ...]

BUILD_DIR comes from build.py, using the same interpreter as this script.  CASE
is a key of CASES in test_thread_safety.py, all of them by default.  On a
free-threaded interpreter, the GIL is turned off in the driver.

Exits with 1 if a case has a report, fails or times out.  With --no-locks,
reports are expected instead: see README.md.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import glob
import os
import pathlib
import re
import site
import subprocess
import sys
import sysconfig
import time

HERE = pathlib.Path(__file__).resolve().parent


def real_python():
    # TSan fails to load through the bin/python3.x launcher
    macos = pathlib.Path(sys.base_prefix, 'Resources/Python.app/Contents/MacOS')
    name = 'PythonT' if sysconfig.get_config_var('Py_GIL_DISABLED') else 'Python'
    if not (macos / name).exists():
        sys.exit('no %s in %s: a framework build of Python is needed' % (name, macos))
    return str(macos / name)


def child_env(build_dir, tsan_options='abort_on_error=0'):
    runtime = subprocess.run(['clang++', '-print-resource-dir'], capture_output=True,
                             text=True, check=True).stdout.strip()
    # the real binary does not see the packages of a virtual environment
    path = [os.path.join(build_dir, 'src')] + site.getsitepackages()
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(path),
               DYLD_INSERT_LIBRARIES=runtime + '/lib/darwin/libclang_rt.tsan_osx_dynamic.dylib',
               TSAN_OPTIONS=tsan_options)
    if sysconfig.get_config_var('Py_GIL_DISABLED'):
        env['PYTHON_GIL'] = '0'
    return env


def run_case(case, args):
    tag = re.sub(r'\W', '_', case)
    for f in glob.glob(os.path.join(args.out_dir, tag + '.*')):
        os.remove(f)
    # abort_on_error defaults to 1 on macOS: every run with a report would abort
    # at exit and open a crash window
    env = child_env(args.build_dir,
                    'abort_on_error=0 halt_on_error=0 history_size=4 '
                    'second_deadlock_stack=1 report_signal_unsafe=0 log_path='
                    + os.path.join(args.out_dir, tag + '.tsan'))
    if args.no_locks:
        env['TSAN_NO_LOCKS'] = '1'
    start = time.time()
    output = os.path.join(args.out_dir, tag + '.out')
    with open(output, 'w') as out:
        p = subprocess.Popen([real_python(), str(HERE / 'driver.py'), case],
                             stdout=out, stderr=subprocess.STDOUT, env=env)
        try:
            rc = p.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
            rc = 'TIMEOUT'
    reports = 0
    for f in glob.glob(os.path.join(args.out_dir, tag + '.tsan*')):
        with open(f, errors='replace') as log:
            reports += log.read().count('WARNING: ThreadSanitizer')
    with open(output, errors='replace') as out:
        lines = out.read().strip().splitlines()
    last = lines[-1][:70] if lines else ''
    # TSan exits with 66 when it has reported something
    failed = rc not in (0, 66) or last != 'errors: none'
    print('%-18s %6.0fs rc=%-7s reports=%-5d %s' % (case, time.time() - start, rc,
                                                    reports, last), flush=True)
    return reports, failed


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('build_dir', type=os.path.abspath)
    parser.add_argument('out_dir', type=os.path.abspath)
    parser.add_argument('cases', nargs='*', metavar='CASE')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=1800, help='per case, in seconds')
    parser.add_argument('--no-locks', action='store_true',
                        help='positive control: disable the Python locks')
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    cases = args.cases
    if not cases:
        cases = subprocess.run([real_python(), str(HERE / 'driver.py'), '--list'],
                               capture_output=True, text=True, check=True,
                               env=child_env(args.build_dir)).stdout.splitlines()
    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(lambda case: run_case(case, args), cases))

    reports = sum(r for r, _ in results)
    failures = sum(f for _, f in results)
    print('%d cases, %d reports, %d failed' % (len(cases), reports, failures))
    if args.no_locks:
        return 1 if failures or not reports else 0
    return 1 if failures or reports else 0


if __name__ == '__main__':
    sys.exit(main())
