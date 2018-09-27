# STREAM: Solution Transfers for REActor Multiphysics

## General Compilation

STREAM must be compiled in a directory that contains the Nek5000 input deck (SIZE and \<casename\>.usr files).  We
refer to this directory as the **case directory**.  Compilation follows this general pattern:

  1. `cd` into the case directory

  2. Run CMake for the given Nek5000 case:

     `$ CC=mpicc CXX=mpicxx FC=mpif90 cmake -DCASENAME=<casename> $HOME/stream`

  3. Run `make` for a desired driver:

     `make <driver>`

The availabe targets are described below.

## Tests

The tests can be compiled by entering the `tests/singlerod/short/` directory, configuring with CMake, and running `make` for
one or more of the desired targets:

``` Console
$ cd tests/singlerod/short/
$ CC=mpicc CXX=mpicxx FC=mpif90 cmake ../../../
$ make test_openmc_singlerod
$ make test_nek5000_singlerod
$ make test_coupled_singlerod
```

## Running

The tests can be run using `mpirun` with a given number of processes (`-np`) and OpenMP threads (`OMP_NUM_THREADS`).
The number of OpenMP threads does not affect `test_nek5000_singlerod`, since it is run with MPI only.

``` Console
$ OMP_NUM_THREADS=4 mpirun -np 4 ./test_openmc_singlerod
$ mpirun -np 4 ./test_nek5000_singlerod
$ OMP_NUM_THREADS=4 mpirun -np 4 ./test_coupled_singlerod
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
