#!/usr/bin/env bash
set -ex

cd tests/singlerod/short/openmc_nekrs
export NEKRS_HOME=$(realpath ../build/install)
cat rod_short.par
../build/install/bin/nrspre rod_short 2
mpirun -np 2 ../build/install/bin/enrico
