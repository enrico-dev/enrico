#!/usr/bin/env bash
set -ex

if [ "$MODE" = "openmc_nek5000" ]; then
  cd tests/singlerod/short/nek5000
  gunzip -f rodcht.run01.gz
  cp rodcht.run01 "../$MODE/"
fi
