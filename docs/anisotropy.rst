********************
Anisotropy
********************

``ttcrpy`` can compute traveltimes in anisotropic media, and return the
sensitivity of those traveltimes to the parameters describing the medium.  Both
are available on 2D rectilinear grids (:class:`ttcrpy.rgrid.Grid2d`), on 2D
triangular meshes (:class:`ttcrpy.tmesh.Mesh2d`) and, for the models with a
vertical symmetry axis, on 3D rectilinear grids
(:class:`ttcrpy.rgrid.Grid3d`) and tetrahedral meshes
(:class:`ttcrpy.tmesh.Mesh3d`).

Anisotropy is described cell by cell, so it requires ``cell_slowness=True``, and
it is implemented for the Shortest-Path Method only, ``method='SPM'``.

Models
======

The model is chosen with the ``aniso`` argument of the constructor, and its
parameters are given with the setters listed below.  The order of the setters is
also the order in which the blocks of the sensitivity matrix appear, so it is
worth reading the table with :ref:`the section on the Jacobian <jacobian>`.

In two dimensions
-----------------

=======================  ============================================================================
``aniso``                parameters, in the order the setters take them
=======================  ============================================================================
``iso``                  ``set_slowness``
``elliptical``           ``set_slowness``, ``set_xi``
``vti_sh``               ``set_Vs0``, ``set_gamma``
``tilted_elliptical``    ``set_slowness``, ``set_xi``, ``set_tilt_angle``
``tti_sh``               ``set_Vs0``, ``set_gamma``, ``set_tilt_angle``
``weakly_anelliptical``  ``set_slowness``, ``set_s2``, ``set_s4``
``vti_psv``              ``set_Vp0``, ``set_Vs0``, ``set_epsilon``, ``set_delta``
``tti_psv``              ``set_Vp0``, ``set_Vs0``, ``set_epsilon``, ``set_delta``, ``set_tilt_angle``
=======================  ============================================================================

Note that ``set_slowness`` means different things from one model to the next.
For an elliptical medium it is the slowness along the horizontal axis, while for
a weakly anelliptical one it is the slowness along the vertical axis.

In three dimensions
-------------------

:class:`ttcrpy.rgrid.Grid3d` and :class:`ttcrpy.tmesh.Mesh3d` take the four
models whose symmetry axis is vertical, with the same ``aniso`` names and the
same setters.  On a mesh each setter takes one value per tetrahedron.  The
tilted models have no 3D counterpart yet: a tilted axis in three dimensions
needs a dip *and* an azimuth, which is a different model rather than the same
one with an extra parameter.

=======================  ============================================================================
``aniso``                parameters, in the order the setters take them
=======================  ============================================================================
``iso``                  ``set_slowness``
``vti_sh``               ``set_Vs0``, ``set_gamma``
``elliptical``           ``set_slowness``, ``set_chi``, ``set_psi``
``weakly_anelliptical``  ``set_slowness``, ``set_s2``, ``set_s4``
``vti_psv``              ``set_Vp0``, ``set_Vs0``, ``set_epsilon``, ``set_delta``
=======================  ============================================================================

``vti_sh``, ``vti_psv`` and ``weakly_anelliptical`` describe the same media as
in 2D, the polar angle being measured from the vertical axis and computed from
the horizontal magnitude :math:`\sqrt{l_x^2 + l_y^2}` of the ray segment.
Being transversely isotropic they are symmetric in azimuth, so the traveltime
depends on the direction of a segment only through that polar angle.

``elliptical`` is the exception: in 3D it is an **ellipsoid with two
independent ratios**, more general than a transversely isotropic medium.  It is
described by the *vertical* slowness :math:`s_z` and by

.. math::

   \chi = \frac{s_x}{s_z}, \qquad \psi = \frac{s_y}{s_z}

so that the traveltime along a segment of components :math:`(l_x, l_y, l_z)` is

.. math::

   dt = s_z \sqrt{\chi^2 l_x^2 + \psi^2 l_y^2 + l_z^2}

the axes of the ellipsoid being aligned with the global axes.  Beware that
``set_slowness`` therefore takes the **vertical** slowness here, where the 2D
``elliptical`` takes the horizontal one, and that :math:`\psi` is the ratio
along :math:`y`, not the ray angle it denotes elsewhere on this page.

Elliptical
----------

The medium is described by the horizontal slowness :math:`s_x` and by the
anisotropy ratio :math:`\xi = s_z / s_x`, the traveltime along a segment of
components :math:`(l_x, l_z)` being

.. math::

   dt = s_x \sqrt{l_x^2 + \xi^2 l_z^2}

Transversely isotropic, SH wave
-------------------------------

The medium is described by the vertical S-wave velocity :math:`V_{S0}` and by
Thomsen's :math:`\gamma` (Thomsen, 1986).  The phase velocity, for a phase angle :math:`\theta`
measured from the symmetry axis, is

.. math::

   v(\theta) = V_{S0} \sqrt{1 + 2 \gamma \sin^2\theta}

which is an ellipse, so the traveltime again has a closed form.

Transversely isotropic, qP and qSV waves
----------------------------------------

The medium is described by the two vertical velocities :math:`V_{P0}` and
:math:`V_{S0}` and by Thomsen's :math:`\epsilon` and :math:`\delta`.  The
compressional and shear waves are coupled, their phase velocities being the two
roots of one quartic.  The expression used is the exact one, so the medium is
not restricted to weak anisotropy (Thomsen, 1986; Tsvankin, 2012).

The wave to model is chosen with ``set_phase``:

.. code-block:: python

   g.set_phase('qSV')     # or 'qP', which is what a medium gives until asked

The quasi-compressional wave is the one modelled until ``set_phase`` is called.
The two other arguments the C++ ``setPhase()`` takes, ``1`` for qP and anything
else for qSV, are accepted as well.

Weakly anelliptical
-------------------

The formulation of Rommel (2024), in which the energy velocity is expanded in
powers of :math:`\sin^2\theta`,

.. math::

   v(\theta) = v_0 \left[ 1 + \left( s_2 + s_4 \sin^2\theta \right)
                          \sin^2\theta \right]

:math:`v_0` being the vertical velocity and :math:`s_2` and :math:`s_4` the
second- and fourth-order coefficients.

Tilted media
------------

The tilted models are the same media with the symmetry axis rotated by the angle
given to ``set_tilt_angle``, in radians and one value per cell.  A segment of
components :math:`(l_x, l_z)` is expressed in the frame of the symmetry axis,

.. math::

   t_1 = l_x \cos\theta_t + l_z \sin\theta_t, \qquad
   t_2 = l_z \cos\theta_t - l_x \sin\theta_t

:math:`t_2` being the component along the symmetry axis.  With this convention a
ray making an angle :math:`\psi` with the vertical makes an angle
:math:`\psi + \theta_t` with the symmetry axis, which therefore lies at
:math:`-\theta_t` from the vertical.  A tilt of zero reproduces the
corresponding vertical model exactly.

Phase and group velocity
========================

In an anisotropic medium the traveltime along a ray segment is its length
divided by the **group** (or energy) velocity taken in the direction of the
segment, not by the phase velocity.  The two differ everywhere except along the
symmetry directions (Červený, 2001), the group slowness being the support
function of the phase slowness surface,

.. math::

   \frac{1}{V_g(\psi)} = \max_{\theta} \frac{\cos(\psi-\theta)}{v(\theta)}
   \;\ge\; \frac{1}{v(\psi)}

so that using the phase velocity in the direction of the segment would make the
medium appear faster than it is.

For the elliptical, SH and weakly anelliptical media the group velocity has a
closed form.  For the coupled qP and qSV media it does not, and it is tabulated
against the ray angle when the medium parameters are set.  Where the
group-velocity surface of the qSV wave is multivalued, at a triplication, the
tabulation keeps the fastest branch, which is the first arrival.

.. _jacobian:

The Jacobian
============

Passing ``compute_L=True`` to ``raytrace`` returns, besides the traveltimes, the
matrix of partial derivatives of the traveltimes with respect to the parameters
of the medium.  It has one row per source-receiver pair and **one block of**
``ncells`` **columns per parameter**, the blocks appearing in the order the
setters are listed in the table above:

.. math::

   L = \left[\; \frac{\partial t}{\partial p_1}
        \;\middle|\; \frac{\partial t}{\partial p_2}
        \;\middle|\; \cdots \;\right]

So for a ``tti_psv`` medium on a grid of 1600 cells, ``L`` has
:math:`5 \times 1600 = 8000` columns, holding in turn the derivatives with
respect to :math:`V_{P0}`, :math:`V_{S0}`, :math:`\epsilon`, :math:`\delta` and
the tilt angle.

.. code-block:: python

   import numpy as np
   import ttcrpy.rgrid as rg

   x = np.arange(0., 2.05, 0.05)
   z = np.arange(0., 2.05, 0.05)
   ncells = (x.size-1) * (z.size-1)

   g = rg.Grid2d(x, z, method='SPM', nsnx=10, nsnz=10, aniso='vti_psv')
   g.set_Vp0(np.full(ncells, 3.094))
   g.set_Vs0(np.full(ncells, 1.51))
   g.set_epsilon(np.full(ncells, 0.256))
   g.set_delta(np.full(ncells, -0.0505))

   src = np.array([[0.3, 0.3]])
   rcv = np.array([[1.7, 1.6], [1.8, 0.6]])

   tt, L = g.raytrace(src, rcv, compute_L=True)

   dt_dVp0 = L[:, :ncells]                 # first block
   dt_dVs0 = L[:, ncells:2*ncells]         # second block

How it is obtained
------------------

Under Fermat's principle the raypath is stationary, so only the explicit
dependence of the traveltime on the medium matters.  Writing the group velocity
through the phase velocity of the branch carrying the arrival,
:math:`V_g = v / \cos(\psi - \theta_s)`, and noting that :math:`\theta_s` is
stationary for :math:`\cos(\psi-\theta)/v(\theta)`, the envelope theorem removes
the derivative of :math:`\theta_s` and leaves, for a segment of length
:math:`\ell`,

.. math::

   \frac{\partial}{\partial p} \frac{\ell}{V_g}
     = -\frac{\ell}{V_g\, v} \frac{\partial v}{\partial p}

evaluated on that branch.  This holds whichever branch carries the arrival and
is therefore valid through the triplications of the qSV wave.  The derivative
with respect to the tilt angle follows from the same relation, the ray angle
entering it as the angle to the symmetry axis:

.. math::

   \frac{\partial}{\partial \psi} \frac{\ell}{V_g}
     = -\frac{\ell\, v'}{v\, V_g}

Checking a Jacobian
-------------------

The traveltime is homogeneous of degree one in a slowness and of degree minus
one in a velocity, and the group velocity of a coupled medium scales with both
of its vertical velocities together.  Euler's identities therefore give exact
checks, which need no finite difference:

.. code-block:: python

   # a medium described by a slowness
   assert np.allclose(L[:, :ncells] @ slowness, tt)

   # a medium described by a velocity
   assert np.allclose(L[:, :ncells] @ Vs0, -tt)

   # a coupled medium, whose group velocity scales with both velocities
   assert np.allclose(L[:, :ncells] @ Vp0 + L[:, ncells:2*ncells] @ Vs0, -tt)

Restrictions
============

``compute_L`` is available for the Shortest-Path and Dynamic Shortest-Path
methods, with the slowness defined at the cells.  It is refused for the Fast
Sweeping Method, and for a mesh whose slowness is defined at the nodes.

Pickling a grid or a mesh keeps its geometry, the options it was built with, the
anisotropy model and the phase, but **not the parameters of the medium**.  A
grid or mesh rebuilt from a pickle has to be given them again:

.. code-block:: python

   import pickle

   g2 = pickle.loads(pickle.dumps(g))
   g2.set_Vp0(np.full(ncells, 3.094))      # and the rest

References
==========

- Thomsen, L., 1986. Weak elastic anisotropy. Geophysics, 51(10), 1954-1966.
  DOI : 10.1190/1.1442051
  https://library.seg.org/doi/10.1190/1.1442051

  The parameters :math:`\epsilon`, :math:`\delta` and :math:`\gamma` used here
  to describe transversely isotropic media, and the exact phase velocities from
  which the weak-anisotropy approximations are drawn.

- Tsvankin, I., 2012. Seismic Signatures and Analysis of Reflection Data in
  Anisotropic Media, 3rd edition. Society of Exploration Geophysicists,
  Geophysical References Series No. 19.

  Phase and group velocities of transversely isotropic media, and the coupling
  of the qP and qSV waves.

- Červený, V., 2001. Seismic Ray Theory. Cambridge University Press.

  Ray theory in anisotropic media, and the relation between the phase and the
  group velocity that the traveltime along a ray segment rests on.

- Rommel, B. E., 2024. Weakly anelliptical traveltime analysis: Ambiguity
  between subsurface and elasticity. Geophysics, 89(4), C171-C182.
  DOI : 10.1190/geo2023-0274.1
  https://library.seg.org/doi/10.1190/geo2023-0274.1

  The weakly anelliptical model, and the notation this package follows for it.
  The companion material is at
  https://github.com/bjornrommel/steinkauz

See also the :ref:`references` describing the algorithms themselves.
