# ENRICO: Exascale Nuclear Reactor Investigative COde

[![License](https://img.shields.io/github/license/enrico-dev/enrico.svg)](http://enrico-docs.readthedocs.io/en/latest/license.html)
[![Travis CI build status (Linux)](https://travis-ci.org/enrico-dev/enrico.svg?branch=master)](https://travis-ci.org/enrico-dev/enrico)

## General Compilation

ENRICO can be configured and compiled in an arbitrary build directory.  However, at configure time, you must point to
the location of the Nek5000 input deck using the `-DUSR_LOC` flag.  Here is a sample workflow, for the compiling the
`tests/singlerod/short/` case:

  1. Create and enter a build directory anywhere.  For example:
  ``` Console
    $ mkdir tests/singlerod/short/build
    $ cd tests/singlerod/short/build
  ```
  2. Enter the build directory and run CMake, using `USR_LOC` to specify the location of the Nek5000 input deck:
  ``` Console
    $ CC=mpicc CXX=mpicxx FC=mpif90 cmake -DUSR_LOC=../nek5000 ../../../../
  ```
  3. Make it:
  ``` Console
    $ make -j4 enrico
  ```

# Library Targets

In a given case directory, you may also compile coupled and single-physics libraries without a driver.  The available
targets are:
* `enrico`: The ENRICO executable
* `libopenmc`: The OpenMC library
* `libnek5000`: The Nek5000 library

# Nek5000 Utilities

Nek5000 provides several utilities for working with mesh data.  These may be compiled using the following
targets:
* `genmap`: The Nek5000 tool for mesh partitioning (generates a .map file)
* `genbox`: The Nek5000 tool for generating a simple box mesh (generates a .box file)
