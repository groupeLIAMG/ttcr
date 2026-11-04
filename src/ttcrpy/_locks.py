# -*- coding: utf-8 -*-
"""Locks protecting the C++ state of a single grid or mesh object.

A grid keeps its working data (traveltimes at the nodes, temporary nodes, GPU
solvers) in one set per thread number, while slowness and the other model
parameters are shared by all of them.  This gives two rules:

- calls using the same thread number must not overlap, as they would overwrite
  each other's traveltimes;
- a call changing the model must not overlap with anything else.

Each thread number therefore has its own lock.  A call working in one thread
number holds that lock, and a call changing the model holds all of them, so
calls with different thread numbers still run in parallel.

These are Python locks because a thread waiting for one releases the GIL.  A
thread waiting for a C++ mutex would keep it, and could block the thread holding
the mutex for good.

Locks are always taken in ascending order.  A call holding a single lock must
therefore never call a method that needs all of them: to change the model, take
all the locks from the start.
"""

import functools
import inspect
import operator
import threading
from contextlib import contextmanager


class GridLock:
    """One reentrant lock per thread number of a grid.

    Parameters
    ----------
    n_threads : int
        number of threads of the grid
    """

    def __init__(self, n_threads):
        self._locks = tuple(threading.RLock() for _ in range(max(1, n_threads)))

    @property
    def n_threads(self):
        """int: number of thread numbers guarded"""
        return len(self._locks)

    @contextmanager
    def one(self, thread_no=0):
        """Hold the lock of one thread number."""
        with self._locks[thread_no]:
            yield

    @contextmanager
    def all(self):
        """Hold the locks of all thread numbers."""
        n_held = 0
        try:
            for lock in self._locks:
                lock.acquire()
                n_held += 1
            yield
        finally:
            for lock in reversed(self._locks[:n_held]):
                lock.release()


def exclusive(method):
    """Decorator for a method that must not overlap with any other call.

    The object must hold its GridLock in the attribute ``_lock``, as must the
    objects of the other decorators in this module.
    """
    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        with self._lock.all():
            return method(self, *args, **kwargs)
    return wrapper


def reader(method):
    """Decorator for a method that only reads the model.

    Holding one lock is enough to keep out the calls changing the model, since
    these need all of them.
    """
    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        with self._lock.one():
            return method(self, *args, **kwargs)
    return wrapper


def by_thread_no(*writes):
    """Decorator for a method that has a ``thread_no`` argument.

    The method holds the lock of that thread number.  It holds all the locks if
    ``thread_no`` is None, because the work is then spread over all threads, or
    if it was given an argument changing the model.

    With a single thread, ``thread_no`` is passed on without being checked,
    because the only lock available is the right one whatever its value.

    Parameters
    ----------
    writes : str
        names of the arguments that change the model when they are not None

    Raises
    ------
    ValueError
        if the grid has many threads and ``thread_no`` is not one of them
    """
    def decorator(method):
        params = inspect.signature(method).parameters
        # position in *args, which does not hold self
        position = {name: i - 1 for i, name in enumerate(params)}
        default = params['thread_no'].default

        def get(name, args, kwargs, missing):
            if name in kwargs:
                return kwargs[name]
            if position[name] < len(args):
                return args[position[name]]
            return missing

        @functools.wraps(method)
        def wrapper(self, *args, **kwargs):
            lock = self._lock
            thread_no = get('thread_no', args, kwargs, default)
            if (thread_no is None
                    or any(get(w, args, kwargs, None) is not None
                           for w in writes)):
                held = lock.all()
            elif lock.n_threads == 1:
                held = lock.one()
            else:
                index = operator.index(thread_no)
                if not 0 <= index < lock.n_threads:
                    raise ValueError(
                        f'thread_no should be between 0 and '
                        f'{lock.n_threads - 1}, got {thread_no}')
                held = lock.one(index)
            with held:
                return method(self, *args, **kwargs)
        return wrapper
    return decorator
