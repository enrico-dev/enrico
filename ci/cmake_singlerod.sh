#!/usr/bin/env bash
set -ex

mkdir -p tests/singlerod/short/build
cd tests/singlerod/short/build
cmake -DUSR_LOC=../ ../../../..
