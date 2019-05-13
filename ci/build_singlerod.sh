#!/usr/bin/env bash

set -ex

cd tests/singlerod/short/build
make -j -l4 enrico
make -j -l4 test_nek5000_singlerod
make -j -l4 test_openmc_singlerod
