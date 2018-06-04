# STREAM: Solution Transfers for REActor Multiphysics

# General Compilation

STREAM must be compiled in the directory that contains the Nek5000 input deck (SIZE and \<casename\>.usr files).  We 
refer to this directory as the **case directory**.  Compilation follows this general pattern:

  1. `cd` into the case directory

  2. Run CMake for the given Nek5000 case:

     `$ CC=mpicc CXX=mpicxx F90=mpif90 cmake -DCASENAME=<casename> $HOME/stream`
   
  3. Run `make` for a desired driver:

     `make <driver>`
     
The availabe targets are described below.

# Tests

## Available API tests

The available tests run OpenMC and Nek5000 via the single-physics and coupled C++ drivers.  They do not implement
solution transfer.  
* `tests/openmc_api`: Instantiates and runs an `OpenmcDriver`.
* `tests/nek5000_api`: Instantiates and runs a `NekDriver`.
* `tests/coupled_driver`: Instantiates an `OpenmcNekDriver`.  Runs the `OpenmcDriver` and `NekDriver` member objects.

## Compiling and Running API Tests

Currently, the API tests have input decks from the Nek5000 "ethier" case.  Since the API tests simply run the
C++ drivers without coupling, the particular input deck is not significant.

### OpenmcDriver Test
``` Console
$ cd tests/openmc_driver
$ CC=mpicc CXX=mpicxx F90=mpif90 cmake -DCASENAME=ethier ../../
$ make test_openmc_driver
$ mpirun -np 4 ./test_openmc_driver
```

### NekDriver Test
``` Console
$ cd tests/nek_driver
$ CC=mpicc CXX=mpicxx F90=mpif90 cmake -DCASENAME=ethier ../../
$ make test_nek_driver
$ mpirun -np 4 ./test_nek_driver
```

### CoupledDriver Test
``` Console
$ cd tests/coupled_driver
$ CC=mpicc CXX=mpicxx F90=mpif90 cmake -DCASENAME=ethier ../../
$ make test_coupled_driver
$ mpirun -np 4 ./test_coupled_driver
```

# Library Targets
 
In a given case directory, you may also compile coupled and single-physics libraries without a driver.  The available 
targets are:
* `stream`: The STREAM library
* `libopenmc`: The OpenMC library
* `libnek5000`: The Nek5000 library

# Nek5000 targets

Nek5000 provides several utilities for working with mesh data.  These may be compiled using the following
targets:
* `genmap`: The Nek5000 tool for mesh partitioning (generates a .map file)
* `genbox`: The Nek5000 tool for generating a simple box mesh (generates a .box file)
