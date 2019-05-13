#!/usr/bin/env bash
set -ex

cd tests/singlerod/short
gunzip -f rodcht.run01.gz
cp rodcht.run01 nek5000/
cp rodcht.run01 openmc_nek5000/
