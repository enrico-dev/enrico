#!/usr/bin/env bash
set -e

if [ "$MODE" = "openmc_nek5000" ]; then
  curdir="$(pwd)"
  source ci/test_singlerod_openmc.sh
  cd "$curdir"
  source ci/test_singlerod_nek5000.sh
  cd "$curdir"
  source ci/test_singlerod_openmc_nek5000.sh
elif [ "$MODE" = "openmc_nekrs" ]; then
  source ci/test_singlerod_openmc_nekrs.sh
elif [ "$MODE" = "openmc_heat_surrogate" ]; then
  source ci/test_singlerod_openmc_heat_surrogate.sh
else
  echo "Invalid test mode (provided MODE=\"$MODE\""
  exit 1
fi
