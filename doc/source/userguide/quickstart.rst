.. _userguide_quickstart:

=================
Quick Start Guide
=================

This describes usage for a generic UNIX-based platform.  See ":ref:`userguide_platform_specific`"
for instructions for several specific platforms.  

Requirements
------------

ENRICO itself requires the following:

* C, C++, and Fortran compilers
* MPI libraries for C and Fortran
* CMake v3.3 or higher

ENRICO's component solvers have additional dependencies:

* OpenMC: Requires HDF5 v1.8, v1.10, or v1.12
* nekRS: Requires GNU compilers
* SHIFT: See SHIFT docs for full requirements


Configuring, Building, and Installing
-------------------------------------

The general workflow is:

1. Create a build directory and enter it::

    $ mkdir build
    $ cd build

2. Configure with CMake, specifying the desired heat/fluids solver.     
   
  - With Nek5000.  The path to Nek5000's case-specific compile-time input files must be specified by the `USR_LOC`
    CMake variable::

    $ CC=mpicc CXX=mpicxx FC=mpifort cmake -DNEK_DIST=nek5000 -DUSR_LOC=../tests/singlerod/short/nek5000 ..

  - With nekRS::

    $ CC=mpicc CXX=mpicxx FC=mpifort cmake -DNEK_DIST=nekrs ..

  - With heat surrogate solver only::

    $ CC=mpicc CXX=mpicxx FC=mpifort cmake ..
    
  Note that ENRICO always installs OpenMC and the heat surrogate solver; and that SHIFT is installed if 
  the SHIFT source is available in the source tree.  Any installed solver can be selected at
  runtime (see ":ref:`userguide_input`" for details).

  Other commonly-used CMake variables (such as ``CMAKE_BUILD_TYPE``, ``CMAKE_INSTALL_PREFIX``,
  etc.) are also supported. 


3. Run make and install.  By default (i.e., if ``CMAKE_INSTALL_PREFIX`` was not specified), the 
   executables will be in ``install/bin``::

    $ make -j4 enrico install

4. Set environment variables

  - **Required for nekRS:** Set ``NEKRS_HOME`` to the absolute path of the installation directory::

    $ export NEKRS_HOME=$(realpath install)

  - Optional for all builds: Add the installation directory to the current path::

    $ export PATH=$(realpath install/bin):$PATH

Running a Case
--------------

ENRICO must be run from the directory containing the case's input files.  This includes the input
files for the selected single-physics solvers; and the ENRICO-specific ``enrico.xml`` input file.  
See ":ref:`userguide_input`" for details about ``enrico.xml``.

For example, the full-length fuel rod test case can be run as follows. These commands assume you
have added ``build/install/bin`` to the ``PATH`` environment variable; if not, you must refer to 
the full path to the ``enrico`` executable.

- For OpenMC + Nek5000.  ``rod_l.run03.tgz`` contains the restart solution for Nek5000::

    $ cd tests/singlerod/long/openmc_nek5000
    $ tar -xzf rod_l.run03.tgz    
    $ mpirun -np 64 enrico

- For OpenMC + NekRS. ``rod_l.run03.tgz`` contains the restart solution for nekRS::

    $ cd tests/singlerod/long/openmc_nekrs
    $ tar -xzf rod_l.run03.tgz    
    $ mpirun -np 64 enrico

- For OpenMC + heat surrogate::

    $ cd tests/singlerod/long/openmc_heat_surrogate
    $ mpirun -np 64 enrico
