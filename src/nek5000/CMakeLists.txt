cmake_minimum_required(VERSION 3.1)
project(Nek5000 Fortran C)

# =============================================================================
# Configure options
# =============================================================================

option(CASENAME    "Name of problem case (.usr and .rea files)")
option(MPI         "Use MPI (default: ON)"                          ON)
option(PROFILING   "Use internal profiling (default: ON)"           ON)
option(NOMPIIO     "Deactivate MPI-IO support (default: OFF)"       OFF)
option(BGQ         "Use BGQ optimized mxm (default: OFF"            OFF)
option(XSMM        "Use libxsmm for mxm (default: OFF)"             OFF)
option(CVODE       "Use CVODE support (default: OFF)"               OFF)
option(VENDOR_BLAS "Use VENDOR BLAS/LAPACK (default: OFF)"          OFF)
option(EXTBAR      "Use overlapping mesh solver (default: OFF)"     OFF)
option(NEKNEK      "Add underscore to exit() (default: OFF)"        OFF)
option(NEKNEK      "USE overlapping mesh solver (default: OFF)"     OFF)
option(CMTNEK      "Use DG compressible-flow solver (default: OFF)" OFF)

# if(NOT CASENAME)
#   message(FATAL_ERROR "Must provide a case name via '-DCASENAME=<casename>'")
# endif()
# =============================================================================
# Compiler Options
# =============================================================================

# Detect name-mangling scheme
include(FortranCInterface)
FortranCInterface_VERIFY()
if(FortranCInterface_GLOBAL_SUFFIX STREQUAL "_")
  add_definitions(-DUNDERSCORE)
endif()

# Set fixed-format source
# set(CMAKE_Fortran_FORMAT "FIXED")

# Use preprocessor and 8-byte reals
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -cpp")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -Mpreprocess")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -fpconstant -fpp")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL XL)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qrealsize=8 -qdpc=e -qsuffix=cpp=f -qsuppress=cmpmsg")
else()
  message(FATAL_ERROR "Specified compiler (${CMAKE_Fortran_COMPILER_ID}) is not supported by Nek5000")
endif()

# =============================================================================
# Preprocessor Symbols
# =============================================================================

if(PROFILING)
  add_definitions(-DTIMER)
endif()

if(MPI)
  add_definitions(-DMPI)
endif()

if(NOT MPI OR NOMPIIO)
  add_definitions(-DNOMPIIO)
endif()

if(BGQ)
  add_definitions(-DBGQ)
endif()

if(XSMM)
  add_definitions(-DXSMM)
endif()

if(CVODE)
  add_definitions(-DCVODE)
endif()

if(VENDOR_BLAS)
  add_definitions(-DVENDOR_BLAS)
endif()

if(NEKNEK)
  add_definitions(-DNEKNEK)
endif()

if(CMTNEK)
  add_definitions(-DCMTNEK)
endif()

# =============================================================================
# Header files
# =============================================================================

include_directories(${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR}
    ${CMAKE_CURRENT_LIST_DIR}/core)

# TODO: Make mpif.h a proper dependency
if(MPI)
  file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/mpif.h)
else()
  configure_file( ${CMAKE_CURRENT_LIST_DIR}/core/mpi_dummy.h ${CMAKE_CURRENT_BINARY_DIR}/mpif.h )
endif()

if(CMT)
  include_directories(${CMAKE_CURRENT_LIST_DIR}/core/cmt/)
endif()

# =============================================================================
# SESSION File
# =============================================================================

file(WRITE ${CMAKE_BINARY_DIR}/SESSION.NAME "${CASENAME}\n${CMAKE_BINARY_DIR}\n")

# =============================================================================
# Nek5000 Sources and Targets
# =============================================================================

# Create <casename>.f from <casename>.usr
set(USR_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/${CASENAME}.f)
ADD_CUSTOM_COMMAND(
    OUTPUT ${USR_SOURCE}
    COMMAND ${CMAKE_COMMAND}
    -DCASENAME=${CASENAME} -DINFILE_DIR=${CMAKE_BINARY_DIR} -DOUTFILE_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -P ${CMAKE_CURRENT_LIST_DIR}/core/mkuserfile.cmake
    DEPENDS ${CMAKE_BINARY_DIR}/${CASENAME}.usr
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set(CORE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/dictionary.c
    ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/finiparser.c
    ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/iniparser.c
    ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/nek_in_situ.f
    ${CMAKE_CURRENT_LIST_DIR}/core/bdry.f
    ${CMAKE_CURRENT_LIST_DIR}/core/byte.c
    ${CMAKE_CURRENT_LIST_DIR}/core/byte_mpi.f
    ${CMAKE_CURRENT_LIST_DIR}/core/calcz.f
    ${CMAKE_CURRENT_LIST_DIR}/core/chelpers.c
    ${CMAKE_CURRENT_LIST_DIR}/core/coef.f
    ${CMAKE_CURRENT_LIST_DIR}/core/conduct.f
    ${CMAKE_CURRENT_LIST_DIR}/core/connect1.f
    ${CMAKE_CURRENT_LIST_DIR}/core/connect2.f
    ${CMAKE_CURRENT_LIST_DIR}/core/convect.f
    ${CMAKE_CURRENT_LIST_DIR}/core/cvode_driver.f
    ${CMAKE_CURRENT_LIST_DIR}/core/drive1.f
    ${CMAKE_CURRENT_LIST_DIR}/core/drive2.f
    ${CMAKE_CURRENT_LIST_DIR}/core/dssum.f
    ${CMAKE_CURRENT_LIST_DIR}/core/edgec.f
    ${CMAKE_CURRENT_LIST_DIR}/core/eigsolv.f
    ${CMAKE_CURRENT_LIST_DIR}/core/fast3d.f
    ${CMAKE_CURRENT_LIST_DIR}/core/fasts.f
    ${CMAKE_CURRENT_LIST_DIR}/core/gauss.f
    ${CMAKE_CURRENT_LIST_DIR}/core/genbox.f
    ${CMAKE_CURRENT_LIST_DIR}/core/genxyz.f
    ${CMAKE_CURRENT_LIST_DIR}/core/gfdm_op.f
    ${CMAKE_CURRENT_LIST_DIR}/core/gfdm_par.f
    ${CMAKE_CURRENT_LIST_DIR}/core/gfdm_solve.f
    ${CMAKE_CURRENT_LIST_DIR}/core/gfldr.f
    ${CMAKE_CURRENT_LIST_DIR}/core/gmres.f
    ${CMAKE_CURRENT_LIST_DIR}/core/hmholtz.f
    ${CMAKE_CURRENT_LIST_DIR}/core/hpf.f
    ${CMAKE_CURRENT_LIST_DIR}/core/hsmg.f
    ${CMAKE_CURRENT_LIST_DIR}/core/ic.f
    ${CMAKE_CURRENT_LIST_DIR}/core/induct.f
    ${CMAKE_CURRENT_LIST_DIR}/core/intp_usr.f
    ${CMAKE_CURRENT_LIST_DIR}/core/makeq.f
    ${CMAKE_CURRENT_LIST_DIR}/core/makeq_aux.f
    ${CMAKE_CURRENT_LIST_DIR}/core/map2.f
    ${CMAKE_CURRENT_LIST_DIR}/core/math.f
    ${CMAKE_CURRENT_LIST_DIR}/core/mvmesh.f
    ${CMAKE_CURRENT_LIST_DIR}/core/mxm_wrapper.f
    ${CMAKE_CURRENT_LIST_DIR}/core/navier0.f
    ${CMAKE_CURRENT_LIST_DIR}/core/navier1.f
    ${CMAKE_CURRENT_LIST_DIR}/core/navier2.f
    ${CMAKE_CURRENT_LIST_DIR}/core/navier3.f
    ${CMAKE_CURRENT_LIST_DIR}/core/navier4.f
    ${CMAKE_CURRENT_LIST_DIR}/core/navier5.f
    ${CMAKE_CURRENT_LIST_DIR}/core/navier6.f
    ${CMAKE_CURRENT_LIST_DIR}/core/navier7.f
    ${CMAKE_CURRENT_LIST_DIR}/core/navier8.f
    ${CMAKE_CURRENT_LIST_DIR}/core/nek_comm.c
    ${CMAKE_CURRENT_LIST_DIR}/core/nek_interface.f
    ${CMAKE_CURRENT_LIST_DIR}/core/papi.f
    ${CMAKE_CURRENT_LIST_DIR}/core/perturb.f
    ${CMAKE_CURRENT_LIST_DIR}/core/plan4.f
    ${CMAKE_CURRENT_LIST_DIR}/core/plan5.f
    ${CMAKE_CURRENT_LIST_DIR}/core/planx.f
    ${CMAKE_CURRENT_LIST_DIR}/core/postpro.f
    ${CMAKE_CURRENT_LIST_DIR}/core/prepost.f
    ${CMAKE_CURRENT_LIST_DIR}/core/reader_par.f
    ${CMAKE_CURRENT_LIST_DIR}/core/reader_re2.f
    ${CMAKE_CURRENT_LIST_DIR}/core/reader_rea.f
    ${CMAKE_CURRENT_LIST_DIR}/core/speclib.f
    ${CMAKE_CURRENT_LIST_DIR}/core/ssolv.f
    ${CMAKE_CURRENT_LIST_DIR}/core/subs1.f
    ${CMAKE_CURRENT_LIST_DIR}/core/subs2.f
    ${CMAKE_CURRENT_LIST_DIR}/core/vprops.f)

if(VENDOR_BLAS)
  unset(BLAS_SOURCES)
else()
  set(BLAS_SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/blas.f
      ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/dsygv.f)
endif()

if(BGQ)
  set(MXM_SOURCES ${CMAKE_CURRENT_LIST_DIR}/core/mxm_bgq.f)
else()
  set(MXM_SOURCES ${CMAKE_CURRENT_LIST_DIR}/core/mxm_std.f)
endif()

if(NEKNEK)
  set(NEKNEK_SOURCES ${CMAKE_CURRENT_LIST_DIR}/core/multimesh.f)
else()
  set(NEKNEK_SOURCES ${CMAKE_CURRENT_LIST_DIR}/core/singlmesh.f)
endif()

if(MPI)
  set(COMM_MPI_SOURCES ${CMAKE_CURRENT_LIST_DIR}/core/comm_mpi.f)
else()
  set(COMM_MPI_SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/core/comm_mpi.f
      ${CMAKE_CURRENT_LIST_DIR}/core/mpi_dummy.f)
endif()

if(CMT)
  set(CMT_SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/MixtPerf.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/artvisc.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/ausm.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/bc.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/diagnostics.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/diffusive_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/drive1_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/drive2_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/driver3_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/eqnsolver_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/face.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/filters_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/inflow_bc.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/intpdiff.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/outflow_bc.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/step.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/surface_fluxes.f
      ${CMAKE_CURRENT_LIST_DIR}/core/cmt/wall_bc.f )
else()
  unset(CMT_SOURCES)
endif()

if(VISIT)
  set(VISIT_SOURCES ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/visit.f)
else()
  unset(VISIT_SOURCES)
endif()

# Library
add_library(libnek5000
    ${USR_SOURCE} ${CORE_SOURCES} ${BLAS_SOURCES} ${MXM_SOURCES}
    ${NEKNEK_SOURCES} ${COMM_MPI_SOURCES} ${CMT_SOURCES} ${VISIT_SOURCES})
set_target_properties(libnek5000 PROPERTIES OUTPUT_NAME nek5000)

# Driver
add_executable(nek5000 ${CMAKE_CURRENT_LIST_DIR}/core/drive.f)
target_link_libraries(nek5000 PUBLIC libnek5000)

# =============================================================================
# Link to 3rd-party libraries
# =============================================================================

# gslib
# TODO: Check for gslib presence first.  Something like:
# find_library(GSLIB gs PATHS ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/gslib/lib/)
# target_link_libraries(nek5000 PUBLIC ${GSLIB})
if(NOT DEFINED USREXIT)
  set(USREXIT TRUE)
endif()
add_subdirectory(3rd_party/gslib/gslib-1.0.1/)
target_link_libraries(nek5000 PUBLIC gs)

# MPI
# ROR: 2018-03-14: This does not work as intended. It finds the first libraries
# available in the system path, rather than the libraries specified by the MPI
# compiler wrapper in CC and FC. This check would be nice, but it is not strictly
# necessary.

# if(MPI)
#   find_package(MPI REQUIRED)
# endif()

# CVODE
# This will find and link to a prebuild CVODE library.
# TODO: Build CVODE from scratch
if(CVODE)
  set(CVODE_LIBDIR ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/cvode/lib/)
  find_library(SUNDIALS_FCVODE sundials_fcvode PATHS ${CVODE_LIBDIR})
  find_library(SUNDIALS_CVODE sundials_cvode PATHS ${CVODE_LIBDIR})
  if(MPI)
    find_library(SUNDIALS_FNVEC sundials_fnvecparallel PATHS ${CVODE_LIBDIR})
    find_library(SUNDIALS_NVEC sundials_fcvodeparallel PATHS ${CVODE_LIBDIR})
  else()
    find_library(SUNDIALS_FNVEC sundials_fnvecserial PATHS ${CVODE_LIBDIR})
    find_library(SUNDIALS_NVEC sundials_fcvodeserial PATHS ${CVODE_LIBDIR})
  endif()
  target_link_libraries(nek5000 PUBLIC
      ${SUNDIALS_FCVODE} ${SUNDIALS_CVODE}
      ${SUNDIALS_FVEC} ${SUNDIALS_FVEC})
endif()

# BLAS
# TODO: Check for BLAS libraries
# if(VENDOR_BLAS)
# endif()

# VISIT
# TODO: Did visit libraries, but need to do visit headers
if(VISIT)
  target_include_directories(libnek5000 PUBLIC
      ${VISIT_INSTALL}/libsim/V2/include/
      ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/)
  target_include_directories(nek5000 PUBLIC
      ${VISIT_INSTALL}/libsim/V2/include/
      ${CMAKE_CURRENT_LIST_DIR}/core/3rd_party/)
endif()

# =============================================================================
# Tools targets and sources
# =============================================================================
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/tools)
