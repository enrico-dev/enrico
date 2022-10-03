#!/usr/bin/env bash
set -ex

cd tests/core/partial/openmc_heat_surrogate
mpirun -np 2 ../build/install/bin/enrico
