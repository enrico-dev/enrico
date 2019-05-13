#!/usr/bin/env bash
set -ex

cd tests/singlerod/short/openmc_nek5000
mpirun -np 8 ../build/enrico
