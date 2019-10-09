#!/usr/bin/env bash
set -e

if [[ -z "$MODE" ]] || [[ "$MODE" == "openmc_nek5000" ]] ; then
  cd tests/singlerod/short
  patch -N rodcht.rea ci_config/rodcht_rea.diff || true
  patch -N SIZE ci_config/SIZE.diff || true
fi
