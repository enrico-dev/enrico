#!/usr/bin/env bash
set -ex

mkdir -p tests/singlerod/short/build
cd tests/singlerod/short/build
if [ "$MODE" = "openmc_nek5000" ]; then
  cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpif90 \
    -DNEK_DIST=nek5000 -DUSR_LOC=../nek5000 ../../../..
elif [ "$MODE" = "openmc_nekrs" ]; then
  cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpif90 \
    -DNEK_DIST=nekrs ../../../..
else
  cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpif90 \
    -DNEK_DIST=none ../../../..
fi
