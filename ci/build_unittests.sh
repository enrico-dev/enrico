#!/usr/bin/env bash

set -ex

cd tests/singlerod/short/build
make -j4 unittests install
