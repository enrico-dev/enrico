# STREAM: Solution Transfers for REActor Multiphysics

# General Compilation

STREAM must be compiled in a directory that contains the Nek5000 input deck (SIZE and \<casename\>.usr files).  We 
refer to this directory as the **case directory**.  Compilation follows this general pattern:

  1. `cd` into the case directory

  2. Run CMake for the given Nek5000 case:

     `$ CC=mpicc CXX=mpicxx F90=mpif90 cmake -DCASENAME=<casename> $HOME/stream`
   
  3. Run `make` for a desired driver:

     `make <driver>`
     
The availabe targets are described below.

# API Tests

The tests in `tests/api_tests/` run OpenMC and Nek5000 via the single-physics and coupled C++ drivers.  They do not test
solution transfer.  The API tests have input decks from the Nek5000 "ethier" case.

* `tests_openmc_api`: Instantiates and runs an `OpenmcDriver`.
* `tests_nek5000_api`: Instantiates and runs a `NekDriver`.
* `tests_coupled_api`: Instantiates an `OpenmcNekDriver`.  Runs the `OpenmcDriver` and `NekDriver` member objects.

## Compiling

The tests can be compiled by entering the `tests/api_tests/` directory, configuring with CMake, and running `make` for
one or more of the desired targets:

``` Console
$ cd tests/api_tests/
$ CC=mpicc CXX=mpicxx F90=mpif90 cmake -DCASENAME=ethier ../../
$ make test_openmc_api
$ make test_nek5000_api
$ make test_coupled_api
```

## Running

The tests can be run using `mpirun` with a given number of processes (`-np`) and OpenMP threads (`OMP_NUM_THREADS`).
The number of OpenMP threads does not affect `test_nek5000_api`, since it is run with MPI only.  

``` Console
$ OMP_NUM_THREADS=4 mpirun -np 4 ./test_openmc_api
$ mpirun -np 4 ./test_nek5000_api
$ OMP_NUM_THREADS=4 mpirun -np 4 ./test_coupled_api
```

# Library Targets
 
In a given case directory, you may also compile coupled and single-physics libraries without a driver.  The available 
targets are:
* `stream`: The STREAM library
* `libopenmc`: The OpenMC library
* `libnek5000`: The Nek5000 library

# Nek5000 Utilities

Nek5000 provides several utilities for working with mesh data.  These may be compiled using the following
targets:
* `genmap`: The Nek5000 tool for mesh partitioning (generates a .map file)
* `genbox`: The Nek5000 tool for generating a simple box mesh (generates a .box file)
