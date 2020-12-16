Quick Start Guide
=================

Configuring
-----------

ENRICO can be compiled to use Nek5000, NekRS, or neither.  This is controlled by the ``-DNEK_DIST``
CMake option.  The following values are allowed:

- ``-DNEK_DIST=nek5000``: (Default) If compiling with Nek5000, you must also specify the location of the input deck using the ``-DUSR_LOC`` option at configure time
- ``-DNEK_DIST=nekrs``
- ``-DNEK_DIST=none``

With any of the these options, the heat surrogate with also be available.

Building and Installing
-----------------------

To be usable, NekRS **must** be installed (e.g., via ``make install``) in addition to being compiled.
While this is not necessary for other builds, ``make install`` is still fully supported.  The default
installation directory is the ``install/``  subdirectory under the build directory.  It can be changed
using the usual ``-DCMAKE_INSTALL_PREFIX`` option.

The general workflow is:

1. Create a build directory in an arbitrary location and enter it::

    $ mkdir build
    $ cd build

2. Run CMake using the desired Nek distribution

  - For Nek5000::

    $ CC=mpicc CXX=mpicxx FC=mpifort cmake -DNEK_DIST=nek5000 -DUSR_LOC=../tests/singlerod/short/nek5000 ..

  - For NekRS::

    $ CC=mpicc CXX=mpicxx FC=mpifort cmake -DNEK_DIST=nekrs ..

  - Without Nek:

    $ CC=mpicc CXX=mpicxx FC=mpifort cmake ..

3. Run make and install::

    $ make -j4 enrico install

4. Set environment variables

  - **Required for NekRS:** Set ``NEKRS_HOME`` to the absolute path of the installation directory::

    $ export NEKRS_HOME=$(realpath install)

  - Optional for all builds: Add the installation directory to the current path::

    $ export PATH=$(realpath install/bin):$PATH

Running a Case
--------------

ENRICO must be run from the directory containing the case's input files.  This includes the input
files for both physics libraries; and the ENRICO-specific ``enrico.xml`` input file.

For the include short singlerod test case, you can run the simulations as follows. (These assume you
have added ``build/install/bin`` to your ``PATH``; if not, you must refer to the full path to ``enrico``)

- For OpenMC + Nek5000::

    $ cd tests/singlerod/short/openmc_nek5000
    $ mpirun -np 32 enrico

- For OpenMC + NekRS::

    $ cd tests/singlerod/short/openmc_nekrs
    $ mpirun -np 32 enrico

- For OpenMC + heat surrogate::

    $ cd tests/singlerod/short/openmc_heat_surrogate
    $ mpirun -np 32 enrico
