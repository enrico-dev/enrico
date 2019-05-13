#!/usr/bin/env bash
set -ex

cd tests/singlerod/short/openmc_heat_surrogate
mpirun -np 2 ../build/enrico
