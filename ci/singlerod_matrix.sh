#!/usr/bin/env bash
set -e

if [ "$MODE" = "openmc" ]; then
  source ci/test_singlerod_openmc.sh
elif [ "$MODE" = "nek5000" ]; then
  source ci/test_singlerod_nek5000.sh
elif [ "$MODE" = "openmc_nek5000" ]; then
  source ci/test_singlerod_openmc_nek5000.sh
elif [ "$MODE" = "openmc_heat_surrogate" ]; then
  source ci/test_singlerod_openmc_heat_surrogate.sh
else
  echo "Invalid test mode (provided MODE=\"$MODE\""
  exit 1
fi
