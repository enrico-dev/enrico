#!/usr/bin/env bash
set -e

if [ "$MODE" = "openmc_nek5000" ] ; then
  cd tests/singlerod/short/nek5000
  patch -N rodcht.rea rodcht_rea.diff || true
  patch -N SIZE SIZE.diff || true
elif [ "$MODE" = "openmc_nekrs" ] ; then
  cd tests/singlerod/short/nekrs
  patch -N rod_short.par rod_short_par.diff || true
fi
