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

include_directories(${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR})

# TODO: Make mpif.h a proper dependency
if(MPI)
  file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/mpif.h)
else()
  configure_file(
      ${CMAKE_CURRENT_LIST_DIR}/mpi_dummy.h
      ${CMAKE_CURRENT_BINARY_DIR}/mpif.h
  )
endif()

if(CMT)
  include_directories(${CMAKE_CURRENT_LIST_DIR}/cmt/)
endif()

# =============================================================================
# SESSION File
# =============================================================================

file(WRITE ${CMAKE_BINARY_DIR}/SESSION.NAME "${CASENAME}\n${CMAKE_BINARY_DIR}\n")

# =============================================================================
# Sources
# =============================================================================

set(USR_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/${CASENAME}.f)
ADD_CUSTOM_COMMAND(
    OUTPUT ${USR_SOURCE}
    COMMAND ${CMAKE_COMMAND}
            -DCASENAME=${CASENAME} -DINFILE_DIR=${CMAKE_BINARY_DIR} -DOUTFILE_DIR=${CMAKE_CURRENT_BINARY_DIR}
            -P ${CMAKE_CURRENT_LIST_DIR}/mkuserfile.cmake
    DEPENDS ${CMAKE_BINARY_DIR}/${CASENAME}.usr
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

set(CORE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/3rd_party/dictionary.c
    ${CMAKE_CURRENT_LIST_DIR}/3rd_party/finiparser.c
    ${CMAKE_CURRENT_LIST_DIR}/3rd_party/iniparser.c
    ${CMAKE_CURRENT_LIST_DIR}/3rd_party/nek_in_situ.f
    ${CMAKE_CURRENT_LIST_DIR}/bdry.f
    ${CMAKE_CURRENT_LIST_DIR}/byte.c
    ${CMAKE_CURRENT_LIST_DIR}/byte_mpi.f
    ${CMAKE_CURRENT_LIST_DIR}/calcz.f
    ${CMAKE_CURRENT_LIST_DIR}/chelpers.c
    ${CMAKE_CURRENT_LIST_DIR}/coef.f
    ${CMAKE_CURRENT_LIST_DIR}/conduct.f
    ${CMAKE_CURRENT_LIST_DIR}/connect1.f
    ${CMAKE_CURRENT_LIST_DIR}/connect2.f
    ${CMAKE_CURRENT_LIST_DIR}/convect.f
    ${CMAKE_CURRENT_LIST_DIR}/cvode_driver.f
    ${CMAKE_CURRENT_LIST_DIR}/drive1.f
    ${CMAKE_CURRENT_LIST_DIR}/drive2.f
    ${CMAKE_CURRENT_LIST_DIR}/dssum.f
    ${CMAKE_CURRENT_LIST_DIR}/edgec.f
    ${CMAKE_CURRENT_LIST_DIR}/eigsolv.f
    ${CMAKE_CURRENT_LIST_DIR}/fast3d.f
    ${CMAKE_CURRENT_LIST_DIR}/fasts.f
    ${CMAKE_CURRENT_LIST_DIR}/gauss.f
    ${CMAKE_CURRENT_LIST_DIR}/genbox.f
    ${CMAKE_CURRENT_LIST_DIR}/genxyz.f
    ${CMAKE_CURRENT_LIST_DIR}/gfdm_op.f
    ${CMAKE_CURRENT_LIST_DIR}/gfdm_par.f
    ${CMAKE_CURRENT_LIST_DIR}/gfdm_solve.f
    ${CMAKE_CURRENT_LIST_DIR}/gfldr.f
    ${CMAKE_CURRENT_LIST_DIR}/gmres.f
    ${CMAKE_CURRENT_LIST_DIR}/hmholtz.f
    ${CMAKE_CURRENT_LIST_DIR}/hpf.f
    ${CMAKE_CURRENT_LIST_DIR}/hsmg.f
    ${CMAKE_CURRENT_LIST_DIR}/ic.f
    ${CMAKE_CURRENT_LIST_DIR}/induct.f
    ${CMAKE_CURRENT_LIST_DIR}/intp_usr.f
    ${CMAKE_CURRENT_LIST_DIR}/makeq.f
    ${CMAKE_CURRENT_LIST_DIR}/makeq_aux.f
    ${CMAKE_CURRENT_LIST_DIR}/map2.f
    ${CMAKE_CURRENT_LIST_DIR}/math.f
    ${CMAKE_CURRENT_LIST_DIR}/mvmesh.f
    ${CMAKE_CURRENT_LIST_DIR}/mxm_wrapper.f
    ${CMAKE_CURRENT_LIST_DIR}/navier0.f
    ${CMAKE_CURRENT_LIST_DIR}/navier1.f
    ${CMAKE_CURRENT_LIST_DIR}/navier2.f
    ${CMAKE_CURRENT_LIST_DIR}/navier3.f
    ${CMAKE_CURRENT_LIST_DIR}/navier4.f
    ${CMAKE_CURRENT_LIST_DIR}/navier5.f
    ${CMAKE_CURRENT_LIST_DIR}/navier6.f
    ${CMAKE_CURRENT_LIST_DIR}/navier7.f
    ${CMAKE_CURRENT_LIST_DIR}/navier8.f
    ${CMAKE_CURRENT_LIST_DIR}/nek_comm.c
    ${CMAKE_CURRENT_LIST_DIR}/nek_interface.f
    ${CMAKE_CURRENT_LIST_DIR}/papi.f
    ${CMAKE_CURRENT_LIST_DIR}/perturb.f
    ${CMAKE_CURRENT_LIST_DIR}/plan4.f
    ${CMAKE_CURRENT_LIST_DIR}/plan5.f
    ${CMAKE_CURRENT_LIST_DIR}/planx.f
    ${CMAKE_CURRENT_LIST_DIR}/postpro.f
    ${CMAKE_CURRENT_LIST_DIR}/prepost.f
    ${CMAKE_CURRENT_LIST_DIR}/reader_par.f
    ${CMAKE_CURRENT_LIST_DIR}/reader_re2.f
    ${CMAKE_CURRENT_LIST_DIR}/reader_rea.f
    ${CMAKE_CURRENT_LIST_DIR}/speclib.f
    ${CMAKE_CURRENT_LIST_DIR}/ssolv.f
    ${CMAKE_CURRENT_LIST_DIR}/subs1.f
    ${CMAKE_CURRENT_LIST_DIR}/subs2.f
    ${CMAKE_CURRENT_LIST_DIR}/vprops.f)

if(VENDOR_BLAS)
  unset(BLAS_SOURCES)
else()
  set(BLAS_SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/3rd_party/blas.f
      ${CMAKE_CURRENT_LIST_DIR}/3rd_party/dsygv.f)
endif()

if(BGQ)
  set(MXM_SOURCES ${CMAKE_CURRENT_LIST_DIR}/mxm_bgq.f)
else()
  set(MXM_SOURCES ${CMAKE_CURRENT_LIST_DIR}/mxm_std.f)
endif()

if(NEKNEK)
  set(NEKNEK_SOURCES ${CMAKE_CURRENT_LIST_DIR}/multimesh.f)
else()
  set(NEKNEK_SOURCES ${CMAKE_CURRENT_LIST_DIR}/singlmesh.f)
endif()

if(MPI)
  set(COMM_MPI_SOURCES ${CMAKE_CURRENT_LIST_DIR}/comm_mpi.f)
else()
  set(COMM_MPI_SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/comm_mpi.f
      ${CMAKE_CURRENT_LIST_DIR}/mpi_dummy.f)
endif()

if(CMT)
  set(CMT_SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/cmt/MixtPerf.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/artvisc.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/ausm.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/bc.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/diagnostics.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/diffusive_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/drive1_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/drive2_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/driver3_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/eqnsolver_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/face.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/filters_cmt.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/inflow_bc.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/intpdiff.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/outflow_bc.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/step.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/surface_fluxes.f
      ${CMAKE_CURRENT_LIST_DIR}/cmt/wall_bc.f )
else()
  unset(CMT_SOURCES)
endif()

if(VISIT)
  set(VISIT_SOURCES ${CMAKE_CURRENT_LIST_DIR}/3rd_party/visit.f)
else()
  unset(VISIT_SOURCES)
endif()

# =============================================================================
# Targets
# =============================================================================

# Targets for "libnek" (the library) and "nek5000" (the driver)
# ROR: 2018-03-14: Building both "libnek5000" and "nek5000" would involve creating
# two targets named "nek5000", and CMake won't allow that.  Hence, we build libnek
# and nek5000.

# Library
add_library(libnek5000
    ${USR_SOURCE} ${CORE_SOURCES} ${BLAS_SOURCES} ${MXM_SOURCES}
    ${NEKNEK_SOURCES} ${COMM_MPI_SOURCES} ${CMT_SOURCES} ${VISIT_SOURCES})
set_target_properties(libnek5000 PROPERTIES OUTPUT_NAME nek5000)

# Driver
add_executable(nek5000 ${CMAKE_CURRENT_LIST_DIR}/drive.f)
target_link_libraries(nek5000 PUBLIC libnek5000)

# =============================================================================
# Link to 3rd-party libraries
# =============================================================================

# -----------------------------------------------------------------------------
# gslib
# -----------------------------------------------------------------------------

# TODO: Check for gslib presence first.  Something like:
# find_library(GSLIB gs PATHS ${CMAKE_CURRENT_LIST_DIR}/3rd_party/gslib/lib/)
# target_link_libraries(nek5000 PUBLIC ${GSLIB})

target_link_libraries(nek5000 PUBLIC gs)

# -----------------------------------------------------------------------------
# MPI
# -----------------------------------------------------------------------------

# ROR: 2018-03-14: This does not work as intended. It finds the first libraries
# available in the system path, rather than the libraries specified by the MPI
# compiler wrapper in CC and FC. This check would be nice, but it is not strictly
# necessary.

# if(MPI)
#   find_package(MPI REQUIRED)
# endif()

# -----------------------------------------------------------------------------
# CVODE
# -----------------------------------------------------------------------------

# This will find and link to a prebuild CVODE library.
# TODO: Build CVODE from scratch

if(CVODE)
  set(CVODE_LIBDIR ${CMAKE_CURRENT_LIST_DIR}/3rd_party/cvode/lib/)
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

# -----------------------------------------------------------------------------
# BLAS
# -----------------------------------------------------------------------------

# TODO: Check for BLAS libraries
# if(VENDOR_BLAS)
# endif()

# -----------------------------------------------------------------------------
# VISIT
# -----------------------------------------------------------------------------

# TODO: Did visit libraries, but need to do visit headers
if(VISIT)
  target_include_directories(libnek5000 PUBLIC
      ${VISIT_INSTALL}/libsim/V2/include/
      ${CMAKE_CURRENT_LIST_DIR}/3rd_party/)
  target_include_directories(nek5000 PUBLIC
      ${VISIT_INSTALL}/libsim/V2/include/
      ${CMAKE_CURRENT_LIST_DIR}/3rd_party/)
endif()
