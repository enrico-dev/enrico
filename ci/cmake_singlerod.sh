#!/usr/bin/env bash
set -ex

mkdir -p tests/singlerod/short/build
cd tests/singlerod/short/build
cmake -DBUILD_SHARED_LIBS=on -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpif90 \
  -DUSR_LOC=../ ../../../..
