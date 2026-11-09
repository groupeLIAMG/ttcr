# Thread tests under ThreadSanitizer

These scripts look for data races in the C++ code when Python threads share a
grid or a mesh.  They take a long time, so they are run by hand and are not
collected by pytest.  They were written for macOS on arm64, with a framework
build of Python, and have not been tried elsewhere.

```bash
python3.14 tests/tsan/build.py /path/to/tsan-build
python3.14 tests/tsan/run.py /path/to/tsan-build /path/to/out-314

python3.15t tests/tsan/build.py /path/to/tsan-build
python3.15t tests/tsan/run.py /path/to/tsan-build /path/to/out-315t
```

`run.py` prints one line per case and exits with 1 if a case has a report, fails
or times out.  The reports are in `OUT_DIR/<case>.tsan.<pid>`, and the driver's
output in `OUT_DIR/<case>.out`.

Check that TSan can still see a race by disabling the Python locks, so that
threads share a thread number.  This run must produce reports:

```bash
python3.14 tests/tsan/run.py /path/to/tsan-build /path/to/out-control \
    --no-locks "Grid3d_d SPM"
```

## Things that can go wrong

- `DYLD_INSERT_LIBRARIES` must reach the real binary inside the framework, since
  TSan fails to load through the `bin/python3.x` launcher.  `run.py` handles this.
- `abort_on_error` defaults to 1 on macOS, so each run with a report would abort
  at exit and open a crash window.  `run.py` sets it to 0.
- `PYTHON_GIL=0` is a fatal error on a GIL build, so `run.py` sets it only for a
  free-threaded driver.
- With the GIL on, the GIL can order solves that do not overlap and hide a race.
  This is why the driver enlarges the models.
- With the GIL off, TSan cannot see the Python locks, so the driver only runs
  distinct thread numbers there.  Calls relying on the locks are covered with the
  GIL on.

## What this does not cover

TSan only sees the paths that run.  Anisotropic media are covered by cases added
in `driver.py` (SPM only, mostly double precision), not by the thread tests.
There are no 2-D GPU grids, and mesh `compute_D`, `compute_K` and
`data_kernel_straight_rays` are not called.  A GPU failure is not triggered, so the fallback that clears `use_gpu`
does not run.
