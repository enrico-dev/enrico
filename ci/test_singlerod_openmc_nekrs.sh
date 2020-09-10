#!/usr/bin/env bash
set -ex

cd tests/singlerod/short/openmc_nekrs
mpirun -np 2 ../build/install/bin/enrico
