# STREAM: Solution Transfers for REActor Multiphysics

## Compilation

STREAM must be compiled in the directory that contains your Nek5000 input deck (SIZE and <casename>.usr files).  For
this example, we will assume you are compiling in `$HOME/stream/src/nek5000/examples/eddy_uv/`.

The first step is to run cmake.  You must point to the CMakeLists in the top-level STREAM directory.  You must also
provide the name of your Nek5000 case (i.e., the basename for <casename>.usr).  For example:
``` Console
$ CC=mpicc CXX=mpicxx F90=mpif90 cmake -DCASENAME=eddy_uv $HOME/stream
```

The second step is to run `make`.  The available targets are:
* `stream`: The STREAM driver (default target)
* `openmc`: The standalone OpenMC driver
* `libopenmc`: The OpenMC library
* `nek5000`: The standalone OpenMC library
* `libnek5000`: The Nek5000 library
* `genmap`: The Nek5000 tool for mesh partitioning (generates a .map file)
* `genbox`: The Nek5000 tool for generating a simple box mesh (generates a .box file)
