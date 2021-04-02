********************
Performance
********************

3D Rectilinear Grids
====================

Models
------

Performance tests were conducted using two different slowness models: a
layer-cake model and a vertical gradient model.  Analytic solutions exist for
both models, which allows accuracy evaluation.  Besides, tests were done for
three level of discretization: coarse, medium, and fine.  The following figures
show the models.

.. image:: images/layers_coarse.*
   :width: 220
.. image:: images/layers_medium.*
   :width: 220
.. image:: images/layers_fine.*
   :width: 220

.. image:: images/gradient_coarse.*
   :width: 220
.. image:: images/gradient_medium.*
   :width: 220
.. image:: images/gradient_fine.*
   :width: 220

The following figures show the results of the tests.  In these figures, models
are labelled by two letters: "L" or "G" for layers or gradient, and "C", "M" or
"F" for coarse, medium or fine.

It is important to note that for the layers model, slowness values are assigned
to cells, whereas for the gradient model, slowness values are assigned to the
nodes of the grid.

Whole-grid accuracy
-------------------

In this section, the accuracy of the traveltimes computed over the grid nodes
(without using the option to update the traveltimes using the raypaths) is
evaluated.  Error is computed for nodes for which the coordinates are round numbers.

Fast-Sweeping Method
^^^^^^^^^^^^^^^^^^^^

The results are shown first for the FSM.  Accuracy is better for the gradient model,
except for the coarse models.  In the latter case, cells are too large (as thick as
the layers) for the solver to yield satisfying accuracy.

.. image:: images/accuracy_vs_cpu_fsm.*
   :width: 800

Shortest-Path Method
^^^^^^^^^^^^^^^^^^^^

Results for the SPM are shown next.  In the legend, the number next to the model
label is the number of secondary nodes employed.  Increasing this number obviously
has an impact on both accuracy and computation time.  Using 5 secondary nodes
appears to be a good compromise.

.. image:: images/accuracy_vs_cpu_spm.*
   :width: 800

Dynamic Shortest-Path Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Results for the DSPM are shown next, in a rather busy figure.  In the legend, the
first number next to the model label is the number of secondary nodes, the second
number is the number of tertiary nodes, and the last number is the radius of the
sphere containing the tertiary nodes around the source.

.. image:: images/accuracy_vs_cpu_dspm.*
   :width: 800

Results by model
^^^^^^^^^^^^^^^^

The next set of figures contains the accuracy achieved with the three methods for
each model.  In all cases, the lowest errors are obtained with the SPM with 15
secondary nodes (at the cost of very high computation time).  For the gradient
model, the FSM is very competitive for the medium and fine models.  Otherwise,
the DSPM often appears to offer a good compromise.

.. image:: images/accuracy_vs_cpu_lc.*
   :width: 800

.. image:: images/accuracy_vs_cpu_lm.*
   :width: 800

.. image:: images/accuracy_vs_cpu_lf.*
   :width: 800

.. image:: images/accuracy_vs_cpu_gc.*
   :width: 800

.. image:: images/accuracy_vs_cpu_gm.*
   :width: 800

.. image:: images/accuracy_vs_cpu_gf.*
   :width: 800
