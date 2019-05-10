ENRICO: Exascale Nuclear Reactor Investigative COde
===================================================

ENRICO is an application that automates the workflow for solving a coupled
particle transport, heat transfer, and fluid dynamics problem. Currently
supported solvers include the `OpenMC <https://openmc.readthedocs.io>`_ and
Shift Monte Carlo codes and the `Nek5000 <https://nek5000.mcs.anl.gov>`_
computational fluid dynamics code. In addition, several simplified surrogate
physics models are available for obtaining approximate solutions that can be
used either to accelerate convergence or for testing.

The code establishes a mapping between the geometry representation in the Monte
Carlo transport solver and the spectral element mesh used by Nek5000 so that the
output of each code can be used as the input of the other. Heat generation rates
determined by the transport solver are used as heat sources in the CFD
simulation, and the temperature/density fields from the CFD simulation are used
as inputs for the transport solve.

Documentation
-------------

.. toctree::
   :maxdepth: 1

   quickstart
   input
   cppapi
   devguide
   license

Acknowledgment
--------------

Development of ENRICO was supported by the Exascale Computing Project
(17-SC-20-SC), a joint project of the U.S. Department of Energy’s Office of
Science and National Nuclear Security Administration, responsible for delivering
a capable exascale ecosystem, including software, applications, and hardware
technology, to support the nation’s exascale computing imperative.

.. image:: ecp.png
    :height: 100px
    :align: center
