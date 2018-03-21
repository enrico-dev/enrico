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
# Build Nek5000 library and executable
# =============================================================================

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/core)

# =============================================================================
# Build tools
# =============================================================================

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/tools)

# =============================================================================
# Build 3rd-party libraries
# =============================================================================

# -----------------------------------------------------------------------------
# gslib
# -----------------------------------------------------------------------------

# Setting USREXIT is necessary to override gslib's defaults
if(NOT DEFINED USREXIT)
  set(USREXIT TRUE)
endif()
add_subdirectory(3rd_party/gslib/gslib-1.0.1/)

# -----------------------------------------------------------------------------
# CVODE
# -----------------------------------------------------------------------------

# TODO: Build CVODE from scratch
