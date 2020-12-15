# ENRICO: Exascale Nuclear Reactor Investigative COde

ENRICO is an application that automates the workflow for solving a coupled particle transport, heat transfer, and fluid dynamics problem. Individual-physics solvers for particle transport and thermal-fluids are chosen at runtime. Currently supported solvers include the [OpenMC](https://docs.openmc.org/en/stable/) and [Shift](https://www.casl.gov/sites/default/files/docs/CASL-U-2015-0170-000.pdf) Monte Carlo particle transport codes and the [Nek5000](https://nek5000.mcs.anl.gov/) and [nekRS](https://github.com/Nek5000/nekRS) computational fluid dynamics codes. A simple surrogate thermal-fluids solver is also available for testing purposes.

[![License](https://img.shields.io/github/license/enrico-dev/enrico.svg)](http://enrico-docs.readthedocs.io/en/latest/license.html)
[![Travis CI build status (Linux)](https://travis-ci.org/enrico-dev/enrico.svg?branch=master)](https://travis-ci.org/enrico-dev/enrico)

## Configuring

### On General Compute Environments

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
     * With NekRS on OLCF Summit:
     ``` Console
       $ CC=mpicc CXX=mpicxx FC=mpifort cmake -DNEK_DIST=nekrs -DCMAKE_INSTALL_LIBDIR=lib -DOCCA_CXX="g++" -DOCCA_CXXFLAGS="-O2 -ftree-vectorize -funrool-loops -mcpu=native -mtune=native" ..
     ```

  3. Run make and install:
     ``` Console
       $ make -j4 enrico install
     ```
  4. Set environment variables
     * **Required for nekRS:** Set `NEKRS_HOME` to the absolute path of the installation directory:
     ``` Console
       $ export NEKRS_HOME=$(realpath install)
     ```
     * Optional for all builds: Add the installation directory to the current path
     ``` Console
       $ export PATH=$(realpath install/bin):$PATH
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

### General Compute Environments

ENRICO must be run from the directory containing the case's input files.  This includes the input
files for the physics applications; and the ENRICO-specific `enrico.xml` input file.  See the
[documentation](https://enrico-docs.readthedocs.io/en/latest/input.html) for a description of the
`enrico.xml` file.

For the included short singlerod test case, you can run the simulations as follows. (These assume you
have added `build/install/bin` to your `PATH`; if not, you must refer to the full path to `enrico`)

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

### OLCF Summmit with NekRS

To run ENRICO with optimal runtime parameters on OLCF Summit, use the script
`scripts/bsub_enrico_summit.py`.  This is based on the script
`vendor/nekRS/scripts/nrsqsub_summit`, but `bsub_enrico_summit.py` handles
input parameters slightly differently.  In particular, the ENRICO script writes
parameters to the NekRS .par file, rather than passing it on the command line
as in the NekRS script.  Because of the way that ENRICO initializes NekRS, it
is more robust to rely on the .par file with ENRICO runs.  

The script has the following command-line options:

```
usage: bsub_enrico_summit.py [-h] -n NODES -t TIME [-b {CUDA,SERIAL}]
                             [--no-precompile]

optional arguments:
  -h, --help            show this help message and exit
  -n NODES, --nodes NODES
                        Specify the number of NODES (not processes) to use
  -t TIME, --time TIME  Sets the runtime limit of the job
  -b {CUDA,SERIAL}, --backend {CUDA,SERIAL}
                        Sets the OCCA kernel backend [default: CUDA]
  --no-precompile       Skips pre-compile step for NekRS kernels
```

It should be run from the directory containing ENRICO's and the physics
drivers' setup files.  It infers the NekRS casename from the `<casename>`
specified in `enrico.xml`.  It requires the `NEKRS_HOME` environment variable
to be set, as described above.  

