# ENRICO: Exascale Nuclear Reactor Investigative COde

ENRICO is an application that automates the workflow for solving a coupled particle transport, heat transfer, and fluid dynamics problem. Individual-physics solvers for particle transport and thermal-fluids are chosen at runtime. Currently supported solvers include the [OpenMC](https://docs.openmc.org/en/stable/) and [Shift](https://www.casl.gov/sites/default/files/docs/CASL-U-2015-0170-000.pdf) Monte Carlo particle transport codes and the [Nek5000](https://nek5000.mcs.anl.gov/) and [nekRS](https://github.com/Nek5000/nekRS) computational fluid dynamics codes. A simple surrogate thermal-fluids solver is also available for testing purposes.

[![License](https://img.shields.io/github/license/enrico-dev/enrico.svg)](http://enrico-docs.readthedocs.io/en/latest/license.html)
[![Travis CI build status (Linux)](https://travis-ci.org/enrico-dev/enrico.svg?branch=master)](https://travis-ci.org/enrico-dev/enrico)

## Configuring

For the thermal-fluids solver, ENRICO can be compiled to use Nek5000, nekRS, or neither. This is controlled by the `-DNEK_DIST`
CMake option.  The following values are allowed:

  * `-DNEK_DIST=nek5000`: (Default) If compiling with Nek5000, you must also specify the location of the input
  deck using the `-DUSR_LOC` option at configure time
  * `-DNEK_DIST=nekrs`
  * `-DNEK_DIST=none`

With any of the these options, the heat surrogate will still be available.

## Building and Installing

To obtain the necessary source files for building, first clone the repository:
``` Console
$ git clone git@github.com:enrico-dev/enrico.git
```
Next, fetch and checkout the submodules containing the various single-physics
applications:
``` Console
$ cd enrico
$ git submodule update --init --recursive
```

Next, the general workflow is for building and installing is:

  1. Create a build directory in an arbitrary location and enter it:
  ``` Console
    $ mkdir build
    $ cd build
  ```
  2. Run CMake using the desired Nek distribution
     * For Nek5000:
     ``` Console
       $ CC=mpicc CXX=mpicxx FC=mpifort cmake -DNEK_DIST=nek5000 -DUSR_LOC=../tests/singlerod/short/nek5000 ..
     ```
     * For nekRS:
     ``` Console
       $ CC=mpicc CXX=mpicxx FC=mpifort cmake -DNEK_DIST=nekrs ..
     ```
     * Without Nek5000 or nekRS:
     ``` Console
       $ CC=mpicc CXX=mpicxx FC=mpifort cmake -DNEK_DIST=none ..
     ```
  3. Run make and install:
     ``` Console
       $ make -j4 enrico install
     ```

### Important Information

To be usable, nekRS **must** be installed (e.g., via `make install`) in addition to being compiled.
While this is not necessary for other builds, `make install` is still fully supported.  The default
installation directory is the `install/` subdirectory under the build directory.  It can be changed
using the usual `-DCMAKE_INSTALL_PREFIX` option.


## Other Dependencies

* For OpenMC:

You must supply a nuclear data library to run OpenMC. See the
[OpenMC documentation](https://docs.openmc.org/en/stable/usersguide/cross_sections.html)
for instructions.

## Running a Case

ENRICO must be run from the directory containing the case's input files.  This includes the input
files for the physics applications; and the ENRICO-specific `enrico.xml` input file.  See the
[documentation](https://enrico-docs.readthedocs.io/en/latest/input.html) for a description of the
`enrico.xml` file.

For the included short singlerod test case, you can run the simulations as follows. (These assume you
have added `build/` to your `PATH`; if not, you must refer to the full path to `enrico`)

  * For OpenMC + Nek5000:
  ``` Console
  $ cd tests/singlerod/short/openmc_nek5000
  $ gunzip -f ../nek5000/rodcht.run01.gz
  $ mpirun -np 32 enrico
  (This particular example requires a restart file, `rodcht.run01`.)
  ```
  * For OpenMC + nekRS
  ``` Console
  $ cd tests/singlerod/short/openmc_nekrs
  $ mpirun -np 32 enrico
  ```
  * For OpenMC + heat surrogate
  ``` Console
  $ cd tests/singlerod/short/openmc_heat_surrogate
  $ mpirun -np 32 enrico
  ```
