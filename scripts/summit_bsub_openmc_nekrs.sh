#!/usr/bin/env bash
set -e

if [ $# -ne 3 ]; then
  echo "usage: $0 <casename> <number of compute nodes> <hh:mm>"
  exit 1
fi

if [ ! -d "$NEKRS_HOME" ]; then
  echo "ERROR: The env variable NEKRS_HOME does not point to an existing directory"
  exit 1
fi

if [ -z "$PROJ_ID" ]; then
  echo "ERROR: The env variable PROJ_ID is not set"
  exit 1
fi

: ${OCCA_CACHE_DIR:="$PWD/.cache/occa"}
NVME_HOME="/mnt/bb/$USER/"
XL_HOME="/sw/summit/xl/16.1.1-3/xlC/16.1.1"

export NEKRS_HOME
export OCCA_CACHE_DIR
export NEKRS_HYPRE_NUM_THREADS=1
export OGS_MPI_SUPPORT=1
export OCCA_CXX="$XL_HOME/bin/xlc" 
export OCCA_CXXFLAGS="-O3 -qarch=pwr9 -qhot -DUSE_OCCA_MEM_BYTE_ALIGN=64" 
export OCCA_LDFLAGS="$XL_HOME/lib/libibmc++.a"
export OMPI_MCA_io=romio321
export ROMIO_HINTS="$(pwd)/.romio_hint"
if [ ! -f "$ROMIO_HINTS" ]; then
  echo "romio_no_indep_rw true"   >$ROMIO_HINTS
  echo "romio_cb_write enable"   >>$ROMIO_HINTS
  echo "romio_ds_write enable"   >>$ROMIO_HINTS
  echo "romio_cb_read enable"    >>$ROMIO_HINTS
  echo "romio_ds_read enable"    >>$ROMIO_HINTS
  echo "cb_buffer_size 16777216" >>$ROMIO_HINTS
  echo "cb_config_list *:1"      >>$ROMIO_HINTS
fi

module unload darshan-runtime
module load gcc cmake cuda hdf5 openblas

case=$1
nodes=$2
time=$3

gpu_per_node=6
cpu_per_rs=7
let nn=$nodes*$gpu_per_node
let ntasks=nn
backend=CUDA

if [ ! -f $case.par ]; then
  echo "Cannot find" $case.par
  exit 1
fi

if [ ! -f $case.co2 ]; then
  echo "Cannot find" $case.co2
  exit 1
fi

if [ ! -f $case.udf ]; then
  echo "Cannot find" $case.udf
  exit 1
fi

if [ ! -f $case.oudf ]; then
  echo "Cannot find" $case.oudf
  exit 1
fi

if [ ! -f $case.re2 ]; then
  echo "Cannot find" $case.re2
  exit 1
fi

mkdir -p $OCCA_CACHE_DIR 2>/dev/null

while true; do
  read -p "Do you want precompile (recommended)? [N]" yn
  case $yn in
    [Yy]* )
      echo $NEKRS_HOME
      mpirun -pami_noib -np 1 $NEKRS_HOME/bin/nekrs --setup $case --build-only $ntasks --backend $backend;
      if [ $? -ne 0 ]; then
        exit 1
      fi
      break ;;
    * )
      break ;;
  esac
done

jsrun="jsrun -X 1 -n$nodes -r1 -a1 -c1 -g0 -b packed:1 -d packed cp -a $OCCA_CACHE_DIR/* $NVME_HOME; export OCCA_CACHE_DIR=$NVME_HOME; jsrun --smpiargs='-gpu' -X 1 -n$nn -r$gpu_per_node -a1 -c$cpu_per_rs -g1 -b rs -d packed $NEKRS_HOME/bin/enrico --setup $case --backend $backend --device-id 0" 

cmd="bsub -nnodes $nodes -alloc_flags NVME -W $time -P $PROJ_ID -J nekRS_$case \"${jsrun}\""

echo $cmd
$cmd
