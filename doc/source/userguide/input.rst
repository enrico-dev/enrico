.. _userguide_input:

ENRICO Runtime Settings 
=======================

Parameters and settings for each individual physics code are set normally in
their respective input files. Parameters related to the coupled simulation are
listed in a special ``enrico.xml`` file. This file accepts the following
elements:

``<heat_fluids>``
~~~~~~~~~~~~~~~~~

Base element for heat/fluids driver parameters

``<driver>``
------------

The physics driver for solving fluid and heat transfer equations. Valid options
are "nek5000", "nekrs", and "surrogate".

``<pressure_bc>``
-----------------

The pressure of the outlet boundary condition in units of [MPa].

Nek5000- and nekRS-specific Parameters
--------------------------------------

Under the ``<heat_fluids>`` element, these sub-elements are available Nek5000 and nekRS

* ``<casename>``: Required. The casename for the problem (i.e, the basename for .par or .rea files).
* ``<output_heat_source>``: Optional. Can be ``true`` or ``false`` (default ``false``).  If true, output the heat
  source to a Nek field file in units of [W/cm^3].  In an ENRICO simulation, the field file ``<casename>#.f#####`` always
  contains temperature. The heat source is additionally output as follows:
  - For Nek5000 runs, the heat source is output as the first passive scalar in ``<casename>#.f######``.
  - For nekRS runs, the heat source is output as the temperature field in a second field file, ``qsc<casename>#.f#####``.


Surrogate-specific Parameters
-----------------------------

Under the ``<heat_fluids>`` element, these surrogate-specific sub-elements are available:

* ``<clad_inner_radius>``: The cladding inner radius in units of [cm].
* ``<clad_outer_radius>``: The cladding outer radius in units of [cm].
* ``<pellet_radius>``: The fuel pellet radius in units of [cm].
* ``<fuel_rings>``: The number of rings the fuel pellet should be subdivided
  into when solving the heat equation.
* ``<clad_rings>``: The number of rings in the cladding should be subdivided
  into when solving the heat equation.
* ``<n_pins_x>``: Number of pins in the assembly in the x-direction.
* ``<n_pins_y>``: Number of pins in the assembly in the y-direction.
* ``<pin_pitch>``: Pitch, or distance between centers along the x- and y-axes,
  between pins. The pitch must be greater than the outer diameter of the pins,
  which would correspond to touching pins. This pitch is used to determine the
  pin-pin spacing and the pin- to assembly-edge spacing, which is taken to be
  half a pitch.
* ``<z>``: Values along the z-axis that subdivide the fuel region in units of [cm].
* ``<inlet_temperature>``: Fluid inlet temperature in [K].
* ``<mass_flowrate>``: Fluid mass flowrate in [kg/s].
* ``<max_subchannel_its>`` * Maximum number of iterations to perform in the
  solution of the subchannel equations. Convergence is based on the relative
  change measured in the 1-norm in enthalpy and pressure between two
  successive iterations. This defaults to 100.
* ``<subchannel_tol_h>``: Convergence tolerance to use for enthalpy between
  two successive iterations of the subchannel solver. This defaults to a
  value of 1e-2.
* ``<subchannel_tol_p>``: Convergence tolerance to use for pressure between
  two successive iterations of the subchannel solver. This defaults to a value
  of 1e-2.
* ``<heat_tol>``: Tolerance on the heat equation solver. This defaults to a value of 1e-4.
* ``<verbosity>``: Degree of output printing for diagnostic checking. This
  defaults to `none`, but may be set to `low` and `high`. Both `low` and `high`
  perform error checks such as ensuring conservation of mass and energy, while
  `high` prints some subchannel solution metrics for each channel.
* ``<viz>``: This element indicates visualization settings for the heat solver.

  - ``filename`` (attribute): File prefix for output VTK files
  - ``<iterations>``: what iterations to write output at
  - ``<resolution>``: resolution of the VTK objects. When fluid regions are
    included, the resolution must be divisible by the number of channels per rod
    (typically 4)
  - ``<data>``: what data to write. Either "all", "source", "temperature", or "density".
  - ``<regions>``: what regions to write output for. Either "all", "solid", or "fluid".

``<neutronics>``
~~~~~~~~~~~~~~~~

Base element for neutron transport driver parameters

``<driver>``
------------

The physics driver for solving particle transport. Valid options are "openmc",
"shift", and "surrogate".

Shift-specific Parameters
-------------------------

Under the ``<neutronics>`` element, these Shift-specific sub-elements are available:

* ``<filename>``: Path to the Shift XML input file

``<coupling>``
~~~~~~~~~~~~~~

Base node for coupling parameters

``<verbose>``
-------------

A boolean ("true" or "false", default "false").  If true, print detailed info, including
MPI communicator layouts for each rank and volume comparisons for each cell (summarized
volume statistics are always printed).

``<power>``
-----------

The power of the reactor in units of [W].

``<max_timesteps>``
-------------------

The maximum number of timesteps.

``<max_picard_iter>``
---------------------

The maximum number of Picard iterations within a timestep.

.. _epsilon:

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
heat source at iteration :math:`i` and :math:`\tilde{q}_{i+1}` be the next
estimate of the heat source as determined by the neutronics solver. Then, the
heat source for iteration :math:`i + 1` is:

.. math::
    q_{i+1} = (1 - \alpha) q_i + \alpha \tilde{q}_{i+1}

Choosing :math:`\alpha = 1` corresponds to no underrelaxation. A special value
of "robbins-monro" indicates that Robbins-Monro relaxation is to be used:

.. math::
    q_{i+1} = \frac{1}{i} q_i + \left (1 - \frac{1}{i} \right) \tilde{q}_{i+1}

*Default*: 1.0

``<alpha_T>``
-------------

Underrelaxation parameter used on a temperature update. Let :math:`T_i` be the
temperature at iteration :math:`i` and :math:`\tilde{T}_{i+1}` be the next
estimate of the temperature as determined by the thermal-fluids solver. Then,
the temperature for iteration :math:`i + 1` is:

.. math::
    T_{i+1} = (1 - \alpha_T) T_i + \alpha_T \tilde{T}_{i+1}

Choosing :math:`\alpha_T = 1` corresponds to no underrelaxation. A special value
of "robbins-monro" indicates that Robbins-Monro relaxation is to be used:

.. math::
    T_{i+1} = \frac{1}{i} T_i + \left (1 - \frac{1}{i} \right) \tilde{T}_{i+1}

*Default*: 1.0

``<alpha_rho>``
---------------

Underrelaxation parameter used on a density update update. Let :math:`\rho_i` be
the density at iteration :math:`i` and :math:`\tilde{\rho}_{i+1}` be the next
estimate of the density as determined by the thermal-fluids solver. Then, the
density for iteration :math:`i + 1` is:

.. math::
    \rho_{i+1} = (1 - \alpha_\rho) \rho_i + \alpha_\rho \tilde{\rho}_{i+1}

Choosing :math:`\alpha_\rho = 1` corresponds to no underrelaxation. A special
value of "robbins-monro" indicates that Robbins-Monro relaxation is to be used:

.. math::
    \rho_{i+1} = \frac{1}{i} \rho_i + \left (1 - \frac{1}{i} \right) \tilde{\rho}_{i+1}

*Default*: 1.0

``<temperature_ic>``
--------------------

The initial temperature distribution can be determined either from the
neutronics solver or the heat-fluids solver. A value of "neutronics" will use
the temperatures specified in the model for the neutronics solver whereas a
value of "heat_fluids" will use the temperatures specified in the model for the
heat-fluids solver.

*Default*: neutronics

``<density_ic>``
----------------

The initial density distribution can be determined either from the
neutronics solver or the heat-fluids solver. A value of "neutronics" will use
the densities specified in the model for the neutronics solver whereas a
value of "heat_fluids" will use the densities specified in the model for the
heat-fluids solver. Note that this density initial condition strictly refers
to the fluid density - the solid density is constant throughout the simulation,
and is unchanged from the value used in the neutronics input.

*Default*: neutronics

``<convergence_norm>``
----------------------

This element indicates the type of norm to use for convergence checks. At each
Picard iteration, the norm of the difference between the temperature at the
previous and current iterations is compared to the value of :ref:`epsilon` in
order to determine convergence. Valid values for this element are "L1", "L2",
and "Linf".

*Default*: Linf
