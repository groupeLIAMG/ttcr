"""Exercise the wrappers from Python threads, for ThreadSanitizer to watch.

Started by run.py, one case per process: python driver.py (CASE | --list)

Thread numbers stay distinct unless the operation's name says otherwise, since
calls sharing one are serialized by the Python locks.  TSan cannot see those
locks with the GIL off, and would report them as races, so the operations that
rely on them only run with the GIL on.
"""
import os
import pathlib
import sys
import threading

import numpy as np

TESTS_DIR = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TESTS_DIR))
os.chdir(TESTS_DIR)
import test_thread_safety as t

if sys.argv[1] == '--list':
    print('\n'.join(t.CASES))
    sys.exit()

case = sys.argv[1]
N = 4
make, kwargs, ndim = t.CASES[case]
errors = []
gil_on = not hasattr(sys, '_is_gil_enabled') or sys._is_gil_enabled()

if os.environ.get('TSAN_NO_LOCKS'):
    # positive control: without the locks, shared_slot is a known race
    from contextlib import nullcontext
    from ttcrpy._locks import GridLock
    GridLock.one = lambda self, thread_no=0: nullcontext()
    GridLock.all = lambda self: nullcontext()


# solves must last long enough to overlap: where they do not, the GIL orders
# them and hides a race from TSan.  Same extent as the tests, finer sampling.
def big3d(n_threads, **kw):
    x = np.arange(0.0, 8.05, 0.25)
    return t.rg.Grid3d(x, x.copy(), x.copy(), n_threads, **kw)


def big2d(n_threads, **kw):
    # the 2-D FSM timed out at 0.02 under TSan
    x = np.arange(0.0, 8.005, 0.05 if kw.get('method') == 'FSM' else 0.02)
    return t.rg.Grid2d(x, x.copy(), n_threads, **kw)


def big_points(dim):
    step = {3: 0.5, 2: 0.1}[dim]
    axes = [np.arange(0.0, 8.05, step)] * dim
    return np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, dim)


t.mesh_points = big_points
make = {t.grid3d: big3d, t.grid2d: big2d}.get(make, make)
kwargs = dict(kwargs, cell_slowness=kwargs.get('cell_slowness', True))


def threads(work, n=N):
    barrier = threading.Barrier(n)

    def target(i):
        try:
            barrier.wait()
            work(i)
        except Exception as e:
            errors.append('%s: %r' % (work.__name__, e))

    ts = [threading.Thread(target=target, args=(i,)) for i in range(n)]
    for th in ts:
        th.start()
    for th in ts:
        th.join()


def new_grid(n_threads=N):
    g = make(n_threads, **kwargs)
    g.set_slowness(t.as_float(1.0, g))
    return g


src, rcv = t.sources_and_receivers(ndim)
g = new_grid()


def distinct_thread_no(i):
    for _ in range(3):
        g.raytrace(src[i], rcv, thread_no=i)


def distinct_with_rays(i):
    for _ in range(2):
        g.raytrace(src[i], rcv, thread_no=i, return_rays=True)


def distinct_compute_h(i):
    if ndim == 3:
        for _ in range(2):
            g.compute_H(src[i], rcv, thread_no=i)


def shared_slot(i):
    for _ in range(2):
        g.raytrace(src[i], rcv, thread_no=0)


def reading_niter(i):
    for _ in range(3):
        g.raytrace(src[i], rcv, thread_no=i)
        g.get_niter()
        g.get_niterw()


def grid_traveltimes(i):
    g.raytrace(src[i], rcv, thread_no=i)
    g.get_grid_traveltimes(thread_no=i)


def default_thread_no(i):
    for _ in range(2):
        g.raytrace(src[i], rcv)


def changing_slowness(i):
    for _ in range(2):
        if i == 0:
            g.set_slowness(t.as_float(1.0, g))
        else:
            g.raytrace(src[i], rcv, thread_no=i)


def set_verbose_during_solve(i):
    module = t.rg if case.startswith('Grid') else t.tm
    for _ in range(2):
        if i == 0:
            module.set_verbose(0)
        else:
            g.raytrace(src[i], rcv, thread_no=i)


# Without the GIL, only the control with the locks disabled may share a thread
# number: it looks for exactly the race that TSan would report otherwise.
OPS = [distinct_thread_no, distinct_with_rays, distinct_compute_h, reading_niter,
       grid_traveltimes, set_verbose_during_solve]
if gil_on or os.environ.get('TSAN_NO_LOCKS'):
    OPS += [shared_slot]
if gil_on:
    OPS += [default_thread_no, changing_slowness]
for op in OPS:
    print('op', op.__name__, flush=True)
    threads(op)

# the C++ workers of a single call, without then with the thread pool
print('op workers', flush=True)
many_src = np.vstack(src[:N])
many_rcv = np.repeat(rcv[:1], N, axis=0)
g2 = new_grid()
g2.raytrace(many_src, many_rcv)
g2.set_use_thread_pool(True)
g2.raytrace(many_src, many_rcv)
g2.raytrace(many_src, many_rcv, return_rays=True)
# the FSM refuses L on some grids and has its own tests; node meshes have no L
if 'FSM' not in case and not case.startswith('Mesh3d node'):
    try:
        g2.raytrace(many_src, many_rcv, compute_L=True)
    except (NotImplementedError, ValueError) as e:
        print('no compute_L:', e)

print('errors:', errors if errors else 'none', flush=True)
