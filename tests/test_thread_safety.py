# -*- coding: utf-8 -*-
"""Calls sharing a grid from many Python threads must not corrupt each other.

The wrappers release the GIL around the solvers, so calls overlap in C++ on any
interpreter.  On a free-threaded one, run with PYTHON_GIL=0: importing the
modules would otherwise turn the GIL back on, and the Python code in between
would then be serialized.
"""
import threading
import time
import unittest

import numpy as np
from scipy.spatial import Delaunay

from ttcrpy._locks import GridLock, by_thread_no, exclusive, reader
import ttcrpy.rgrid as rg
import ttcrpy.tmesh as tm

N_THREADS = 8
N_REPEATS = 10
TIMEOUT = 300   # seconds, a deadlock fails the test instead of hanging


def run_threads(work, n_threads=N_THREADS):
    """Run work(i) in n_threads threads started together, and re-raise errors."""
    barrier = threading.Barrier(n_threads)
    errors = []

    def target(i):
        try:
            barrier.wait()
            work(i)
        except Exception as e:
            errors.append(e)

    threads = [threading.Thread(target=target, args=(i,), daemon=True)
               for i in range(n_threads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(TIMEOUT)
    if any(t.is_alive() for t in threads):
        raise AssertionError('threads are still running: possible deadlock')
    if errors:
        raise errors[0]


def holds(lock, acquire):
    """Tell whether another thread would have to wait for this lock."""
    result = []

    def target():
        got = acquire(blocking=False)
        result.append(not got)
        if got:
            lock.release()

    t = threading.Thread(target=target)
    t.start()
    t.join(TIMEOUT)
    return result[0]


def held(lock):
    """Which thread numbers of a GridLock are held by the calling thread."""
    return [i for i, l in enumerate(lock._locks) if holds(l, l.acquire)]


class Recorder:
    """Stand-in for a grid, recording which locks are held during a call."""

    def __init__(self, n_threads):
        self._lock = GridLock(n_threads)

    @exclusive
    def write(self):
        """write doc"""
        return held(self._lock)

    @reader
    def read(self):
        return held(self._lock)

    @by_thread_no('slowness', 'xi')
    def trace(self, source, slowness=None, xi=None, thread_no=None):
        """trace doc"""
        return held(self._lock)

    @by_thread_no()
    def get_tt(self, thread_no=0):
        return held(self._lock)


class TestGridLock(unittest.TestCase):

    def test_one_holds_only_that_thread_number(self):
        lock = GridLock(3)
        with lock.one(1):
            self.assertEqual(held(lock), [1])
        self.assertEqual(held(lock), [])

    def test_all_holds_every_thread_number(self):
        lock = GridLock(3)
        with lock.all():
            self.assertEqual(held(lock), [0, 1, 2])
        self.assertEqual(held(lock), [])

    def test_all_is_reentrant(self):
        lock = GridLock(3)
        with lock.all():
            with lock.all():
                with lock.one(2):
                    self.assertEqual(held(lock), [0, 1, 2])
            self.assertEqual(held(lock), [0, 1, 2])
        self.assertEqual(held(lock), [])

    def test_released_after_exception(self):
        lock = GridLock(3)
        with self.assertRaises(KeyError):
            with lock.all():
                raise KeyError
        self.assertEqual(held(lock), [])

    def test_at_least_one_lock(self):
        self.assertEqual(GridLock(0).n_threads, 1)

    def test_no_deadlock_under_contention(self):
        lock = GridLock(4)
        count = [0]

        def work(i):
            for _ in range(500):
                if i % 2:
                    with lock.all():
                        count[0] += 1
                else:
                    with lock.one(i % 4):
                        pass
        run_threads(work)
        self.assertEqual(count[0], 500 * (N_THREADS // 2))


class TestDecorators(unittest.TestCase):

    def test_exclusive(self):
        self.assertEqual(Recorder(3).write(), [0, 1, 2])

    def test_reader(self):
        self.assertEqual(Recorder(3).read(), [0])

    def test_thread_no_positional_and_keyword(self):
        r = Recorder(3)
        self.assertEqual(r.trace('src', None, None, 2), [2])
        self.assertEqual(r.trace('src', thread_no=1), [1])

    def test_thread_no_numpy_integer(self):
        self.assertEqual(Recorder(3).trace('src', thread_no=np.int64(1)), [1])

    def test_no_thread_no_uses_all_threads(self):
        self.assertEqual(Recorder(3).trace('src'), [0, 1, 2])

    def test_default_thread_no(self):
        self.assertEqual(Recorder(3).get_tt(), [0])
        self.assertEqual(Recorder(3).get_tt(2), [2])

    def test_changing_the_model_holds_everything(self):
        r = Recorder(3)
        self.assertEqual(r.trace('src', 1.0, thread_no=1), [0, 1, 2])
        self.assertEqual(r.trace('src', slowness=1.0, thread_no=1), [0, 1, 2])
        self.assertEqual(r.trace('src', xi=1.0, thread_no=1), [0, 1, 2])

    def test_invalid_thread_no(self):
        r = Recorder(3)
        for thread_no in (3, -1):
            with self.assertRaises(ValueError):
                r.trace('src', thread_no=thread_no)

    def test_single_thread_ignores_thread_no(self):
        self.assertEqual(Recorder(1).trace('src', thread_no=5), [0])

    def test_docstring_kept(self):
        self.assertEqual(Recorder.write.__doc__, 'write doc')
        self.assertEqual(Recorder.trace.__doc__, 'trace doc')


def grid3d(n_threads, **kwargs):
    x = np.arange(0.0, 8.5, 1.0)
    return rg.Grid3d(x, x.copy(), x.copy(), n_threads, cell_slowness=True,
                     **kwargs)


def grid2d(n_threads, **kwargs):
    x = np.arange(0.0, 8.5, 0.5)
    return rg.Grid2d(x, x.copy(), n_threads, cell_slowness=True, **kwargs)


def mesh_points(ndim):
    axes = [np.arange(0.0, 8.5, 2.0)] * ndim
    return np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, ndim)


def mesh3d(n_threads, **kwargs):
    nodes = mesh_points(3)
    simplices = Delaunay(nodes).simplices.astype(np.int64)
    return tm.Mesh3d(nodes, simplices, n_threads, cell_slowness=True, **kwargs)


def mesh2d(n_threads, **kwargs):
    nodes = mesh_points(2)
    simplices = Delaunay(nodes).simplices.astype(np.int64)
    return tm.Mesh2d(nodes, simplices, n_threads, cell_slowness=True, **kwargs)


# name -> (constructor, keywords, dimension)
CASES = {
    'Grid3d_d SPM': (grid3d, dict(method='SPM', nsnx=2, nsny=2, nsnz=2), 3),
    'Grid3d_f SPM': (grid3d, dict(method='SPM', nsnx=2, nsny=2, nsnz=2,
                                  dtype=np.float32), 3),
    'Grid3d_d FSM': (grid3d, dict(method='FSM'), 3),
    'Grid3d_f FSM': (grid3d, dict(method='FSM', dtype=np.float32), 3),
    'Grid3d_d DSPM': (grid3d, dict(method='DSPM', n_secondary=2,
                                   n_tertiary=2), 3),
    'Grid3d_f FSM GPU': (grid3d, dict(method='FSM', fsm_gpu=True,
                                      dtype=np.float32), 3),
    'Grid2d_d SPM': (grid2d, dict(method='SPM', nsnx=3, nsnz=3), 2),
    'Grid2d_d FSM': (grid2d, dict(method='FSM'), 2),
    'Grid2d_f SPM': (grid2d, dict(method='SPM', nsnx=3, nsnz=3,
                                  dtype=np.float32), 2),
    'Grid2d_f FSM': (grid2d, dict(method='FSM', dtype=np.float32), 2),
    'Mesh3d SPM': (mesh3d, dict(method='SPM', n_secondary=2), 3),
    # tt_from_rp: on this coarse mesh these two did not return from
    # getTraveltimeFromRaypath within a minute, even from a single thread
    'Mesh3d FSM': (mesh3d, dict(method='FSM', tt_from_rp=False), 3),
    'Mesh3d DSPM': (mesh3d, dict(method='DSPM', n_secondary=2, n_tertiary=2,
                                 tt_from_rp=False), 3),
    'Mesh2d SPM': (mesh2d, dict(method='SPM', n_secondary=2), 2),
    'Mesh2d FSM': (mesh2d, dict(method='FSM'), 2),
    'Mesh2d DSPM': (mesh2d, dict(method='DSPM', n_secondary=2,
                                 n_tertiary=2), 2),
}


def sources_and_receivers(ndim):
    src = [np.array([[1.0 + 0.5 * i] + [1.0] * (ndim - 1)])
           for i in range(N_THREADS)]
    rcv = np.array([[7.0] * ndim, [2.0] * (ndim - 1) + [7.0],
                    [7.0] + [2.0] * (ndim - 1)])
    return src, rcv


def as_float(value, grid):
    """Slowness array of the type the grid expects."""
    return value * np.ones(grid.get_number_of_cells())


class TestSharedGrid(unittest.TestCase):

    def assert_tt_equal(self, tt, expected):
        np.testing.assert_allclose(np.ravel(tt), np.ravel(expected), rtol=1e-5)

    def for_each_case(self, check):
        for name, (make, kwargs, ndim) in CASES.items():
            with self.subTest(name):
                check(make, kwargs, ndim)

    def test_default_thread_no(self):
        def check(make, kwargs, ndim):
            g = make(1, **kwargs)
            g.set_slowness(as_float(1.0, g))
            src, rcv = sources_and_receivers(ndim)
            ref = [g.raytrace(s, rcv) for s in src]

            def work(i):
                for _ in range(N_REPEATS):
                    self.assert_tt_equal(g.raytrace(src[i], rcv), ref[i])
            run_threads(work)
        self.for_each_case(check)

    def test_distinct_thread_no(self):
        def check(make, kwargs, ndim):
            g = make(N_THREADS, **kwargs)
            g.set_slowness(as_float(1.0, g))
            src, rcv = sources_and_receivers(ndim)
            ref = [g.raytrace(s, rcv, thread_no=i) for i, s in enumerate(src)]

            def work(i):
                for _ in range(N_REPEATS):
                    self.assert_tt_equal(g.raytrace(src[i], rcv, thread_no=i),
                                         ref[i])
            run_threads(work)
        self.for_each_case(check)

    def test_shared_thread_no(self):
        def check(make, kwargs, ndim):
            g = make(2, **kwargs)
            g.set_slowness(as_float(1.0, g))
            src, rcv = sources_and_receivers(ndim)
            ref = [g.raytrace(s, rcv, thread_no=0) for s in src]

            def work(i):
                for _ in range(N_REPEATS):
                    self.assert_tt_equal(
                        g.raytrace(src[i], rcv, thread_no=i % 2), ref[i])
            run_threads(work)
        self.for_each_case(check)

    def test_all_threads_at_once(self):
        def check(make, kwargs, ndim):
            g = make(2, **kwargs)
            g.set_slowness(as_float(1.0, g))
            src, rcv = sources_and_receivers(ndim)
            ref = g.raytrace(src[0], rcv)

            def work(i):
                for _ in range(N_REPEATS):
                    self.assert_tt_equal(g.raytrace(src[0], rcv), ref)
            run_threads(work)
        self.for_each_case(check)

    def test_slowness_passed_to_raytrace(self):
        # the slowness must stay in place for the whole of the call using it
        def check(make, kwargs, ndim):
            g = make(N_THREADS, **kwargs)
            src, rcv = sources_and_receivers(ndim)
            ref = {f: g.raytrace(src[0], rcv, slowness=as_float(f, g),
                                 thread_no=0)
                   for f in (1.0, 2.0)}

            def work(i):
                f = 1.0 + i % 2
                s = as_float(f, g)
                for _ in range(N_REPEATS):
                    self.assert_tt_equal(
                        g.raytrace(src[0], rcv, slowness=s, thread_no=i), ref[f])
            run_threads(work)
        self.for_each_case(check)

    def test_set_slowness_during_raytrace(self):
        # the result must come from one slowness model, not from a mix of both
        def check(make, kwargs, ndim):
            g = make(N_THREADS, **kwargs)
            src, rcv = sources_and_receivers(ndim)
            models = [as_float(1.0, g), as_float(2.0, g)]
            ref = [g.raytrace(src[0], rcv, slowness=m, thread_no=0)
                   for m in models]

            def work(i):
                for n in range(N_REPEATS):
                    if i == 0:
                        g.set_slowness(models[n % 2])
                        continue
                    tt = np.ravel(g.raytrace(src[0], rcv, thread_no=i))
                    if not any(np.allclose(tt, np.ravel(r), rtol=1e-5)
                               for r in ref):
                        raise AssertionError(f'mixed slowness models: {tt}')
            run_threads(work)
        self.for_each_case(check)

    def test_get_slowness_is_never_mixed(self):
        g = grid3d(2, method='SPM', nsnx=2, nsny=2, nsnz=2)
        models = [as_float(1.0, g), as_float(2.0, g)]
        g.set_slowness(models[0])

        def work(i):
            for n in range(200):
                if i == 0:
                    g.set_slowness(models[n % 2])
                else:
                    self.assertEqual(np.unique(g.get_slowness()).size, 1)
        run_threads(work)

    def test_methods_keep_their_documentation(self):
        for cls in (rg.Grid3d_d, rg.Grid3d_f, rg.Grid2d_d, rg.Grid2d_f,
                    tm.Mesh3d, tm.Mesh2d):
            with self.subTest(cls.__name__):
                self.assertIn('raytrace(', cls.raytrace.__doc__)
                self.assertIn('set_slowness(', cls.set_slowness.__doc__)
                self.assertEqual(cls.raytrace.__name__, 'raytrace')


class TestDecoratedMethods(unittest.TestCase):
    """Check each kind of method ends up with the locks it should hold."""

    def test_thread_no_out_of_range(self):
        g = grid3d(2, method='SPM', nsnx=2, nsny=2, nsnz=2)
        g.set_slowness(as_float(1.0, g))
        src, rcv = sources_and_receivers(3)
        with self.assertRaises(ValueError):
            g.raytrace(src[0], rcv, thread_no=2)

    def test_single_thread_still_ignores_thread_no(self):
        g = grid3d(1, method='SPM', nsnx=2, nsny=2, nsnz=2)
        g.set_slowness(as_float(1.0, g))
        src, rcv = sources_and_receivers(3)
        np.testing.assert_array_equal(g.raytrace(src[0], rcv, thread_no=3),
                                      g.raytrace(src[0], rcv))

    def test_nested_calls_do_not_deadlock(self):
        # raytrace(slowness=...) calls set_slowness with the locks already held
        g = grid3d(2, method='SPM', nsnx=2, nsny=2, nsnz=2)
        src, rcv = sources_and_receivers(3)
        results = []

        def target():
            results.append(g.raytrace(src[0], rcv, slowness=as_float(1.0, g)))

        t = threading.Thread(target=target, daemon=True)
        t.start()
        t.join(TIMEOUT)
        self.assertFalse(t.is_alive(), 'possible deadlock')
        self.assertEqual(len(results), 1)

    def test_pickling_gives_a_new_lock(self):
        import pickle
        g = grid3d(2, method='SPM', nsnx=2, nsny=2, nsnz=2)
        g.set_slowness(as_float(1.0, g))
        g2 = pickle.loads(pickle.dumps(g))
        self.assertIsNot(g2._lock, g._lock)
        self.assertEqual(g2._lock.n_threads, 2)


def longest_stall(call):
    """Run call() and time it, with another Python thread ticking meanwhile.

    Returns the duration of the call, and the longest time the other thread went
    without running.  A call holding the GIL stops it for as long as it lasts.
    """
    stop = threading.Event()
    stalls = [0.0]

    def ticker():
        last = time.perf_counter()
        while not stop.is_set():
            now = time.perf_counter()
            stalls[0] = max(stalls[0], now - last)
            last = now

    t = threading.Thread(target=ticker, daemon=True)
    t.start()
    start = time.perf_counter()
    try:
        call()
        duration = time.perf_counter() - start
    finally:
        stop.set()
        t.join(TIMEOUT)
    return duration, stalls[0]


class TestGilReleased(unittest.TestCase):
    """Other Python threads keep running while a solver is at work."""

    MIN_DURATION = 0.1  # seconds, shorter calls would prove nothing

    def check(self, call):
        duration, stall = longest_stall(call)
        self.assertGreater(duration, self.MIN_DURATION,
                           'call too short for this test to mean anything')
        self.assertLess(stall, 0.5 * duration,
                        'the GIL looks held: stalled %.3f s in a %.3f s call'
                        % (stall, duration))

    def big_grid3d(self, n_threads, **kwargs):
        x = np.arange(0.0, 60.5, 1.0)
        g = rg.Grid3d(x, x.copy(), x.copy(), n_threads, cell_slowness=True,
                      method='SPM', nsnx=2, nsny=2, nsnz=2, **kwargs)
        g.set_slowness(np.ones(g.get_number_of_cells()))
        return g

    def big_grid2d(self, n_threads, **kwargs):
        x = np.arange(0.0, 500.5, 1.0)
        g = rg.Grid2d(x, x.copy(), n_threads, cell_slowness=True,
                      method='SPM', nsnx=5, nsnz=5, **kwargs)
        g.set_slowness(np.ones(g.get_number_of_cells()))
        return g

    def big_mesh(self, ndim, n_threads):
        n = {3: 26, 2: 201}[ndim]
        axes = [np.linspace(0.0, 60.0, n)] * ndim
        nodes = np.stack(np.meshgrid(*axes, indexing='ij'),
                         axis=-1).reshape(-1, ndim)
        simplices = Delaunay(nodes).simplices.astype(np.int64)
        make = tm.Mesh3d if ndim == 3 else tm.Mesh2d
        g = make(nodes, simplices, n_threads, cell_slowness=True, method='SPM',
                 n_secondary=3)
        g.set_slowness(np.ones(g.get_number_of_cells()))
        return g

    def points(self, ndim):
        src = np.array([[5.0] * ndim])
        rcv = np.array([[55.0] * ndim, [55.0] + [5.0] * (ndim - 1)])
        return src, rcv

    def test_raytrace_single_thread(self):
        src, rcv = self.points(3)
        for dtype in (np.float64, np.float32):
            with self.subTest(dtype=dtype.__name__):
                g = self.big_grid3d(1, dtype=dtype)
                self.check(lambda: g.raytrace(src, rcv))

    def test_raytrace_in_one_thread_number(self):
        src, rcv = self.points(3)
        g = self.big_grid3d(2)
        self.check(lambda: g.raytrace(src, rcv, thread_no=1))
        self.check(lambda: g.raytrace(src, rcv, thread_no=1, return_rays=True))

    def test_raytrace_in_all_threads(self):
        src, rcv = self.points(3)
        g = self.big_grid3d(2)
        src2 = np.vstack((src, src + 1.0))
        rcv2 = rcv + 0.0
        self.check(lambda: g.raytrace(src2, rcv2))

    def test_raytrace_computing_l(self):
        src, rcv = self.points(3)
        g = self.big_grid3d(1)
        self.check(lambda: g.raytrace(src, rcv, compute_L=True))

    def test_compute_h(self):
        src, rcv = self.points(3)
        g = self.big_grid3d(1)
        self.check(lambda: g.compute_H(src, rcv))

    def test_grid2d(self):
        src, rcv = self.points(2)
        for dtype in (np.float64, np.float32):
            with self.subTest(dtype=dtype.__name__):
                g = self.big_grid2d(1, dtype=dtype)
                self.check(lambda: g.raytrace(src, rcv))

    def test_meshes(self):
        for ndim in (3, 2):
            with self.subTest(ndim=ndim):
                src, rcv = self.points(ndim)
                g = self.big_mesh(ndim, 1)
                self.check(lambda: g.raytrace(src, rcv))

    def test_mesh_compute_h(self):
        src, rcv = self.points(3)
        g = self.big_mesh(3, 1)
        self.check(lambda: g.compute_H(src, rcv))


class TestErrorsWithoutTheGil(unittest.TestCase):
    """A C++ exception gets back to Python and leaves no lock behind."""

    def check(self, g, ndim):
        g.set_slowness(as_float(1.0, g))
        src, rcv = sources_and_receivers(ndim)
        with self.assertRaises(RuntimeError):
            g.compute_H(src[0], rcv, radius_factor=-1.0)
        self.assertEqual(held(g._lock), [])
        g.compute_H(src[0], rcv)

    def test_grid3d(self):
        for dtype in (np.float64, np.float32):
            with self.subTest(dtype=dtype.__name__):
                self.check(grid3d(2, method='SPM', nsnx=2, nsny=2, nsnz=2,
                                  dtype=dtype), 3)

    def test_mesh3d(self):
        self.check(mesh3d(2, method='SPM', n_secondary=2), 3)


if __name__ == '__main__':
    unittest.main()
