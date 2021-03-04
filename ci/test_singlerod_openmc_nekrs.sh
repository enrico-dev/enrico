#!/usr/bin/env bash
set -ex

cd tests/singlerod/short/openmc_nekrs
export NEKRS_HOME=$(realpath ../build/install)
../build/install/bin/nrspre rod_short 2
mpirun -np 2 ../build/install/bin/enrico
