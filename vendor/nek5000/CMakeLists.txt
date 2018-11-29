cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
project(Nek5000 Fortran C)

# =============================================================================
# Configure options
# =============================================================================

option(USR_LOC     "Name of problem case (.usr and .rea files)")
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

if (NOT USR_LOC)
  set(USR_LOC "${CMAKE_BINARY_DIR}")
  message(STATUS "Nek5000 was not provided with -DUSR_LOC.  Attempting to find .usr files in ${USR_LOC}")
endif()

get_filename_component(USR_LOC "${USR_LOC}" ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}")

if (USR_LOC MATCHES ".*\.usr")
  if (EXISTS "${USR_LOC}")
    get_filename_component(CASENAME "${USR_LOC}" NAME_WE)
    get_filename_component(CASEDIR "${USR_LOC}" DIRECTORY)
  else  ()
    message(FATAL_ERROR "Could not locate a Nek5000 .usr file at '${USR_LOC}'. Please check that -DUSR_LOC is correctly set.")
  endif()
else()
  file(GLOB usr_files "${USR_LOC}/*.usr")
  if (usr_files)
    list(GET usr_files 0 first_usr_file)
    get_filename_component(CASENAME "${first_usr_file}" NAME_WE)
	  set(CASEDIR "${USR_LOC}")
    unset(first_usr_file)
  else ()
    message(FATAL_ERROR "Could not locate a Nek5000 .usr file at '${USR_LOC}'. Please check that -DUSR_LOC is correctly set.")
  endif()
endif()

message(STATUS "Nek5000 configured with '${CASENAME}.usr' found in ${CASEDIR}")

# =============================================================================
# Compiler Options
# =============================================================================

# Detect name-mangling scheme
include(FortranCInterface)
FortranCInterface_VERIFY()

# Set fixed-format source
# set(CMAKE_Fortran_FORMAT "FIXED")

# Use preprocessor and 8-byte reals
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -cpp")
  else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -cpp -mcmodel=large")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -Mpreprocess -mcmodel=medium")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -fpconstant -fpp -mcmodel=large")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL XL)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qrealsize=8 -qdpc=e -qsuffix=cpp=f -qsuppress=cmpmsg")
else()
  message(FATAL_ERROR "Specified compiler (${CMAKE_Fortran_COMPILER_ID}) is not supported by Nek5000")
endif()

# =============================================================================
# Header files
# =============================================================================

# Nek5000
include_directories(${CASEDIR} ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR} core)
# TODO: Make mpif.h a proper dependency
if(MPI)
  file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/mpif.h)
else()
  configure_file(core/mpi_dummy.h ${CMAKE_CURRENT_BINARY_DIR}/mpif.h)
endif()
if(CMT)
  include_directories(core/cmt/)
endif()

# =============================================================================
# gslib
# =============================================================================

add_subdirectory(3rd_party/gslib/gslib-1.0.1/)

# =============================================================================
# Nek5000 Sources and Targets
# =============================================================================

# Create <casename>.f from <casename>.usr
set(USR_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/${CASENAME}.f)
ADD_CUSTOM_COMMAND(
    OUTPUT ${USR_SOURCE}
    COMMAND ${CMAKE_COMMAND}
    -DCASENAME=${CASENAME} -DINFILE_DIR=${CASEDIR} -DOUTFILE_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -P ${CMAKE_CURRENT_LIST_DIR}/core/mkuserfile.cmake
    DEPENDS ${CASEDIR}/${CASENAME}.usr
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set(CORE_SOURCES
    core/3rd_party/dictionary.c
    core/3rd_party/finiparser.c
    core/3rd_party/iniparser.c
    core/3rd_party/nek_in_situ.f
    core/bdry.f
    core/byte.c
    core/byte_mpi.f
    core/calcz.f
    core/chelpers.c
    core/coef.f
    core/conduct.f
    core/connect1.f
    core/connect2.f
    core/convect.f
    core/cvode_driver.f
    core/drive1.f
    core/drive2.f
    core/dssum.f
    core/edgec.f
    core/eigsolv.f
    core/fast3d.f
    core/fasts.f
    core/gauss.f
    core/genbox.f
    core/genxyz.f
    core/gfdm_op.f
    core/gfdm_par.f
    core/gfdm_solve.f
    core/gfldr.f
    core/gmres.f
    core/hmholtz.f
    core/hpf.f
    core/hsmg.f
    core/ic.f
    core/induct.f
    core/intp_usr.f
    core/makeq.f
    core/makeq_aux.f
    core/map2.f
    core/math.f
    core/mvmesh.f
    core/mxm_wrapper.f
    core/navier0.f
    core/navier1.f
    core/navier2.f
    core/navier3.f
    core/navier4.f
    core/navier5.f
    core/navier6.f
    core/navier7.f
    core/navier8.f
    core/nek_comm.c
    core/nek_interface.f
    core/nek_interface_types.f
    core/papi.f
    core/perturb.f
    core/plan4.f
    core/plan5.f
    core/planx.f
    core/postpro.f
    core/prepost.f
    core/reader_par.f
    core/reader_re2.f
    core/reader_rea.f
    core/speclib.f
    core/ssolv.f
    core/subs1.f
    core/subs2.f
    core/vprops.f)

if(VENDOR_BLAS)
  unset(BLAS_SOURCES)
else()
  set(BLAS_SOURCES
      core/3rd_party/blas.f
      core/3rd_party/dsygv.f)
endif()

if(BGQ)
  set(MXM_SOURCES core/mxm_bgq.f)
else()
  set(MXM_SOURCES core/mxm_std.f)
endif()

if(NEKNEK)
  set(NEKNEK_SOURCES core/multimesh.f)
else()
  set(NEKNEK_SOURCES core/singlmesh.f)
endif()

if(MPI)
  set(COMM_MPI_SOURCES core/comm_mpi.f)
else()
  set(COMM_MPI_SOURCES
      core/comm_mpi.f
      core/mpi_dummy.f)
endif()

if(CMT)
  set(CMT_SOURCES
      core/cmt/MixtPerf.f
      core/cmt/artvisc.f
      core/cmt/ausm.f
      core/cmt/bc.f
      core/cmt/diagnostics.f
      core/cmt/diffusive_cmt.f
      core/cmt/drive1_cmt.f
      core/cmt/drive2_cmt.f
      core/cmt/driver3_cmt.f
      core/cmt/eqnsolver_cmt.f
      core/cmt/face.f
      core/cmt/filters_cmt.f
      core/cmt/inflow_bc.f
      core/cmt/intpdiff.f
      core/cmt/outflow_bc.f
      core/cmt/step.f
      core/cmt/surface_fluxes.f
      core/cmt/wall_bc.f)
else()
  unset(CMT_SOURCES)
endif()

if(VISIT)
  set(VISIT_SOURCES core/3rd_party/visit.f)
else()
  unset(VISIT_SOURCES)
endif()

# Library
add_library(libnek5000
    ${USR_SOURCE} ${CORE_SOURCES} ${BLAS_SOURCES} ${MXM_SOURCES}
    ${NEKNEK_SOURCES} ${COMM_MPI_SOURCES} ${CMT_SOURCES} ${VISIT_SOURCES})
target_link_libraries(libnek5000 PUBLIC libgs)
set_target_properties(libnek5000 PROPERTIES OUTPUT_NAME nek5000)

# Driver
add_executable(nek5000 core/drive.f)
target_link_libraries(nek5000 PUBLIC libnek5000)

# =============================================================================
# Preprocessor Symbols
# =============================================================================

if(PROFILING)
  target_compile_definitions(libnek5000 PRIVATE -DTIMER)
endif()

if(MPI OR NOT DEFINED MPI)
  set(MPI TRUE)
  target_compile_definitions(libnek5000 PRIVATE -DMPI)
endif()

if(NOT MPI OR NOMPIIO)
  target_compile_definitions(libnek5000 PRIVATE -DNOMPIIO)
endif()

if(FortranCInterface_GLOBAL_SUFFIX STREQUAL "_")
  target_compile_definitions(libnek5000 PRIVATE -DUNDERSCORE)
endif()

if(BGQ)
  target_compile_definitions(libnek5000 PRIVATE -DBGQ)
endif()

if(XSMM)
  target_compile_definitions(libnek5000 PRIVATE -DXSMM)
endif()

if(CVODE)
  target_compile_definitions(libnek5000 PRIVATE -DCVODE)
endif()

if(VENDOR_BLAS)
  target_compile_definitions(libnek5000 PRIVATE -DVENDOR_BLAS)
endif()

if(NEKNEK)
  target_compile_definitions(libnek5000 PRIVATE -DNEKNEK)
endif()

if(CMTNEK)
  target_compile_definitions(libnek5000 PRIVATE -DCMTNEK)
endif()

# =============================================================================
# Link to 3rd-party libraries
# =============================================================================

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
  set(CVODE_LIBDIR core/3rd_party/cvode/lib/)
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
      core/3rd_party/)
  target_include_directories(nek5000 PUBLIC
      ${VISIT_INSTALL}/libsim/V2/include/
      core/3rd_party/)
endif()

# =============================================================================
# Tools targets and sources
# =============================================================================
add_subdirectory(tools)
