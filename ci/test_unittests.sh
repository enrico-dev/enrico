#!/usr/bin/env bash
set -ex

cd tests/unit
../singlerod/short/build/install/bin/unittests
