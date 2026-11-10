.. _getting_started:


###############
Getting started
###############

.. _installing-ttcrpy:

*****************
Installing ttcrpy
*****************

You can use pip to install the package by doing::

  pip install ttcrpy

Python versions
===============

Wheels are built for Python 3.11 to 3.15 on Linux, macOS and Windows,
including the free-threaded builds of Python 3.14 and 3.15 (3.14t and 3.15t),
for which ttcrpy runs without the GIL (see :ref:`threads`).

Requirements
============

ttcrpy needs the following packages, installed along with it by ``pip install ttcrpy``:
  - numpy (https://numpy.org), 2.0 or later
  - scipy (https://www.scipy.org), 1.14 or later

VTK (https://www.vtk.org) is optional.  It is needed only to read and write
VTK files -- ``to_vtk``, the raypaths it can save, and ``builder`` -- and not
for raytracing.  Install it with ttcrpy by doing::

  pip install 'ttcrpy[vtk]'

Without it, those methods raise an ImportError saying so, and everything else
works as usual.

VTK does not yet publish wheels for every free-threaded build: as of VTK
9.7.2, there are none for macOS on Apple silicon (3.14t and 3.15t) nor for
Windows (3.15t).  There, ``pip install 'ttcrpy[vtk]'`` fails to find VTK, but
``pip install ttcrpy`` works, and so does raytracing.

Upgrading to ttcrpy 1.6
=======================

- Python 3.10 is no longer supported.
- ``Mesh3d`` raises a ValueError for a mesh holding tetrahedra without volume.
  Solvers cannot cross such cells: they used to hang, or to return wrong
  traveltimes.  ``scipy.spatial.Delaunay`` produces them for coplanar nodes,
  on a regular lattice or on the flat boundary of a model.  Those on the
  boundary can be removed; inside the mesh, move the nodes or use a mesh
  generator.
- ``raytrace`` and ``compute_H`` raise a ValueError for a ``thread_no`` outside
  0 to ``n_threads``-1, where they could crash or return wrong results
  (``raytrace`` still ignores ``thread_no`` when ``n_threads=1``).

Sparse arrays returned by ttcrpy
==================================

The methods that return sensitivity or derivative matrices -- ``raytrace``
with ``compute_L`` or ``compute_M``, ``compute_D``, ``compute_K`` and
``data_kernel_straight_rays`` -- hand back scipy **sparse arrays**
(``csr_array``), and no longer sparse matrices (``csr_matrix``).

The difference that matters is what the operators mean.  On a sparse array
``*`` multiplies elementwise and ``@`` is the matrix product, where on a
sparse matrix ``*`` was itself the matrix product.  Code written against the
older releases keeps running and quietly computes something else::

    tt = L * s      # a matrix product before, elementwise now
    tt = L @ s      # the matrix product, on either

Nothing else about the returned objects changes: they are indexed, stacked
and converted with ``toarray()`` as before.  scipy's own guide to the
difference is at
https://docs.scipy.org/doc/scipy/reference/sparse.migration_to_sparray.html

.. _threads:

********************************
Running calculations in parallel
********************************

A grid or a mesh can compute in parallel in two ways.

- Give it several threads with ``n_threads`` when creating it, and pass many
  sources to one call: ``raytrace`` spreads them over the threads.
- Share it between Python threads, each calling it with its own
  ``thread_no``, from 0 to ``n_threads``-1, and one source per call.  Each
  thread number has its own traveltimes, so calls with different thread
  numbers run at the same time.

The second way needs no care beyond choosing the thread numbers::

  from concurrent.futures import ThreadPoolExecutor
  from ttcrpy.rgrid import Grid3d

  g = Grid3d(x, y, z, n_threads=4)
  g.set_slowness(slowness)

  def trace(k):
      # thread k computes sources k, k+4, k+8, ...
      return [g.raytrace(src[i:i+1], rcv, thread_no=k)
              for i in range(k, len(src), 4)]

  with ThreadPoolExecutor(4) as pool:
      results = list(pool.map(trace, range(4)))

The grid takes care of the rest:

- calls using the same thread number wait for each other, and a call with
  ``thread_no=None`` uses all the threads, so it waits for the other calls;
- calls that change the model -- the ``set_*`` methods, or ``raytrace`` given
  ``slowness`` or another parameter of the model -- wait until no other call
  is running, and hold the others back meanwhile;
- an out-of-range ``thread_no`` raises a ValueError.  With ``n_threads=1``,
  ``raytrace`` ignores ``thread_no``, but ``compute_H`` still requires 0 or
  None.

The solvers run without the GIL, so Python threads compute in parallel on any
build of Python.  On a free-threaded build (3.14t, 3.15t), the Python code
between the calls runs in parallel as well.

Changing several parameters of an anisotropic model takes one call per
parameter, and a calculation in another thread can run between two of them.
Avoid changing the model while other threads are computing.  On a 2D
rectilinear grid, ``raytrace`` also accepts ``slowness``, ``xi``, ``theta``,
``Vp0``, ``Vs0``, ``delta``, ``epsilon`` and ``gamma``, and sets them and
computes in a single call.

***************
Simple examples
***************

An example showing how easy it is to use the code can be found at
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example_Grid3d.ipynb

A second example illustrating how to run jobs in parallel is given at
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example_tmesh_parallel.ipynb

An example illustrating how to use gmsh to build models with specific geometries is
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example4_circular_shapes.ipynb

Raytracing in anisotropic media is shown in
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example5_anisotropy_elliptical.ipynb

What a first-arrival solver returns for a triplicated qSV wavefront is discussed in
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example6_vti_triplication.ipynb

Examples with a layered model, showing GPU and float/double usage
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example1_3d_layered_model.ipynb
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example1_float32.ipynb
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example1_opencl.ipynb

