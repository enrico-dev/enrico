#!/usr/bin/env bash
set -ex

if [ ! -e $HOME/endfb-vii.1-hdf5/cross_sections.xml ]; then
  cd $HOME
  wget https://anl.box.com/shared/static/9igk353zpy8fn9ttvtrqgzvw1vtejoz6.xz -O - | tar -xvJ
fi
