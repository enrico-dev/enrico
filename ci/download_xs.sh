#!/usr/bin/env bash
set -ex

if [ ! -e $HOME/endf71_multitemp/cross_sections.xml ]; then
  cd $HOME
  wget https://anl.box.com/shared/static/46osfq56h4bd68r6e6pbhupsk4gbnvj2.xz -O - | tar -xvJ
fi
