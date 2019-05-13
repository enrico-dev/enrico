#!/usr/bin/env bash
set -ex

cd tests/singlerod/short/openmc_nek5000
if [ -a rodcht.run01.gz ]; then gunzip -f rodcht.run01.gz ; fi
mpirun -np 8 ../build/test_nek5000_singlerod
