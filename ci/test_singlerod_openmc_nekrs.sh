#!/usr/bin/env bash
set -ex

cd tests/singlerod/short/openmc_nekrs
export NEKRS_HOME=$(realpath ../build/install)
mpirun -np 2 ../build/install/bin/enrico
