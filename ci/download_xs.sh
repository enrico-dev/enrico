#!/usr/bin/env bash
set -ex

if [ "$MODE" = "openmc" ] || [ "$MODE" = "openmc_nek5000" ] || [ "$MODE" = "openmc_heat_surrogate" ]; then
  if [ ! -d $HOME/endf71_multitemp/ ]; then
    cd $HOME
    wget https://anl.box.com/shared/static/46osfq56h4bd68r6e6pbhupsk4gbnvj2.xz -O - | tar -xvJ
  fi
fi
