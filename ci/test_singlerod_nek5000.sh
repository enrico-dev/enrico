#!/usr/bin/env bash
set -ex

cd tests/singlerod/short/nek5000
mpirun -np 8 ../build/test_nek5000_singlerod
