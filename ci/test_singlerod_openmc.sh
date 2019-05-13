#!/usr/bin/env bash
set -ex

cd tests/singlerod/short/openmc_nek5000
mpirun -np 2 ../build/test_openmc_singlerod
