#!/usr/bin/env bash

set -ex

cd tests/singlerod/short/build
make -j -l4 enrico
if [ "$MODE" = "openmc_nek5000" ]; then
  make -j -l4 test_openmc_singlerod
  make -j -l4 test_nek5000_singlerod
fi
