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

Requirements
============

ttcrpy needs the following packages:
  - numpy (https://numpy.org), 2.0 or later
  - scipy (https://www.scipy.org), 1.14 or later
  - vtk (https://www.vtk.org)

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

***************
Simple examples
***************

An example showing how easy it is to use the code can be found at
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example_Grid3d.ipynb

A second example illustrating how to run jobs in parallel is given at
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example_tmesh_parallel.ipynb

An example illutrating how to use gmsh to build models with specific geometries is
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example4_circular_shapes.ipynb

Raytracing in anisotropic media is shown in
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example5_anisotropy_elliptical.ipynb

What a first-arrival solver returns for a triplicated qSV wavefront is discussed in
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example6_vti_triplication.ipynb

Examples with a layered model, showing GPU and float/double usage
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example1_3d_layered_model.ipynb
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example1_float32.ipynb
https://github.com/groupeLIAMG/ttcr/blob/master/examples/example1_opencl.ipynb

