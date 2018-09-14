cmake_minimum_required(VERSION 3.1)
project(gslib C Fortran)

enable_testing()

# =============================================================================
# Preprocessor Symbols
# =============================================================================

if(MPI OR NOT DEFINED MPI)
  set(MPI TRUE)
  add_definitions(-DMPI)
endif()

if(MPI AND NOT NOMPIIO)
  add_definitions(-DUSEMPIIO)
endif()

include(FortranCInterface)
FortranCInterface_VERIFY()
if(FortranCInterface_GLOBAL_SUFFIX STREQUAL "_")
  add_definitions(-DUNDERSCORE)
endif()
# set(CMAKE_Fortran_FORMAT "FIXED")

if(BLAS)
  add_definitions(-DUSE_CBLAS)
else()
  add_definitions(-DUSE_NAIVE_BLAS)
endif()

if(USREXIT)
  add_definitions(-DUSE_USR_EXIT)
endif()

if(NOT DEFINED PREFIX)
  set(PREFIX gslib_)
endif()
add_definitions(-DPREFIX=${PREFIX})

if(NOT DEFINED FPREFIX)
  set(FPREFIX fgslib_)
endif()
add_definitions(-DFPREFIX=${FPREFIX})

add_definitions(-DGLOBAL_LONG_LONG)

# =============================================================================
# Headers
# =============================================================================

include_directories(${CMAKE_CURRENT_LIST_DIR}/src)

# =============================================================================
# Sources and Targets
# =============================================================================

# Adds gslib sources/target
include(${CMAKE_CURRENT_LIST_DIR}/src/CMakeLists.txt)

# For tests
include(${CMAKE_CURRENT_LIST_DIR}/tests/CMakeLists.txt)
