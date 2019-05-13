#!/usr/bin/env bash
set -ex

cd tests/singlerod/short
patch -N rodcht.rea ci_config/rodcht_rea.diff || true
patch -N SIZE ci_config/SIZE.diff || true
