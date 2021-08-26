.. _userguide_platform_specific:

=======================
Platform-Specific Usage 
=======================

OLCF Summit
-----------

OpenMC + nekRS
~~~~~~~~~~~~~~

*Last updated: 2020-06-24*

Building ENRICO on Summit with OpenMC and nekRS requires these modules::

  $ module load gcc cmake cuda hdf5 python/3.7.0-anaconda3-5.3.0 openblas

We recommend running CMake with these variables defined::

  $ CC=mpicc CXX=mpicxx FC=mpifort cmake -DNEK_DIST=nekrs -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_LIBDIR=lib -DOCCA_CXX="g++" -DOCCA_CXXFLAGS="-O2 -ftree-vectorize -funroll-loops -mcpu=native -mtune=native" ..

Then compile and set the environment variables, as described in ":ref:`userguide_quickstart`"::

  $ make -j8 enrico install
  $ export NEKRS_HOME=$(realpath install)
  $ export PATH=$(realpath install/bin):$PATH


To run, you should use the provided script ``summit_bsub_openmc_nekrs.sh``, which is installed in `install/bin`.  
Its usage is::

  summit_bsub_openmc_nekrs.sh <casename> <number of compute nodes> <hh:mm>
