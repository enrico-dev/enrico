#!/usr/bin/env bash
set -e

if [[ -z "$MODE" ]] || [[ "$MODE" == "openmc_nek5000" ]] ; then
  cd tests/singlerod/short/nek5000
  patch -N rodcht.rea rodcht_rea.diff || true
  patch -N SIZE SIZE.diff || true
fi
