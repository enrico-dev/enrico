Coupled Driver Input -- enrico.xml
==================================

Parameters and settings for each individual physics code are set normally in
their respective input files. Parameters related to the coupled simulation are
listed in a special ``enrico.xml`` file. This file accepts the following
elements:

``<driver_heatfluids>``
-----------------------

The physics driver for solving fluid and heat transfer equations. Valid options
are "nek5000" and "surrogate".

``<driver_transport>``
----------------------

The physics driver for solving particle transport. Valid options are "openmc",
"shift", and "surrogate".

``<power>``
-----------

The power of the reactor in units of [W].

``<max_timesteps>``
-------------------

The maximum number of timesteps.

``<max_picard_iter>``
---------------------

The maximum number of Picard iterations within a timestep.

``<epsilon>``
-------------

Convergence criterion, :math:`\epsilon`. If :math:`T_i` and :math:`T_{i+1}` are
the set of temperatures at iterations :math:`i` and :math:`i+1`, convergence is
reached if

.. math::
    \lvert T_{i+1} - T_i \rvert < \epsilon

*Default*: 1.0e-3

``<alpha>``
-----------

Underrelaxation parameter used on a heat source update. Let :math:`q_i` be the
heat source at iteration :math:`i` and :math:`\tilde{q}_{i+1}` be next estimate of
the heat source as determined by the neutronics solver. Then, the heat source
for iteration :math:`i + 1` is:

.. math::
    q_{i+1} = (1 - \alpha) q_i + \alpha \tilde{q}_{i+1}

Choosing :math:`\alpha = 1` corresponds to no underrelaxation.

*Default*: 1.0

``<temperature_ic>``
--------------------

The initial temperature distribution can be determined either from the
neutronics solver or the heat-fluids solver. A value of "neutronics" will use
the temperatures specified in the model for the neutronics solver whereas a
value of "heat" will use the temperatures specified in the model for the
heat-fluids solver.

*Default*: neutronics

``<pressure_bc>``
--------------

The pressure of the system in units of [MPa].

``<nek5000>``
-------------

This element holds settings specific to Nek5000 that are not contained in the
Nek5000 user files.

``<casename>``
~~~~~~~~~~~~~~

The Nek5000 casename.

``<heat_surrogate>``
--------------------

This element holds settings specific to the surrogate heat transfer solver.

``<clad_inner_radius>``
~~~~~~~~~~~~~~~~~~~~~~~

The cladding inner radius in units of [cm].

``<clad_outer_radius>``
~~~~~~~~~~~~~~~~~~~~~~~

The cladding outer radius in units of [cm].

``<pellet_radius>``
~~~~~~~~~~~~~~~~~~~

The fuel pellet radius in units of [cm].

``<fuel_rings>``
~~~~~~~~~~~~~~~~

The number of rings the fuel pellet should be subdivided into when solving the
heat equation.

``<clad_rings>``
~~~~~~~~~~~~~~~~

The number of rings in the cladding should be subdivided into when solving the
heat equation.

``<pin_centers>``
~~~~~~~~~~~~~~~~~

A list of the (x, y) locations of the centers of each fuel pin given in units of
[cm]. For example, if there were two pins at (0, 5) and (3, 2), the input should
be given as:

.. code-block:: xml

    <pin_centers>0.0 5.0 3.0 2.0</pin_centers>

``<z>``
~~~~~~~

Values along the z-axis that subdivide the fuel region in units of [cm].

``<tolerance>``
~~~~~~~~~~~~~~~

Tolerance on the heat equation solver.

``<viz>``
~~~~~~~~~

This element indicates visualization settings for the heat solver. It has the
following attributes:

- `filename`: File prefix for output VTK files

It also has the following subelements:

- ``<iterations>``: what iterations to write output at
- ``<resolution>``: resolution of the VTK objects
- ``<data>``: what data to write. Either "all", "source", or "temperature".
- ``<regions>``: what regions to write output for. Either "all", "fuel", or "cladding".
