cmake_minimum_required(VERSION 3.1)
project(Nek5000 Fortran C)

# =============================================================================
# Default compilers
# =============================================================================

if(NOT MAXNEL)
  add_definitions(-DMAXNEL=150000)
endif()

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
# ROR: 2018-03-14: May need to make this target-specific with Fortran_FORMAT property
# set(CMAKE_Fortran_FORMAT "FIXED")

# Use preprocessor and 8-byte reals
# ROR: 2018-03-14: Maybe possible to use CMAKE's own CMAKE_Fortran_PREPROCESS_SOURCE flags
# see https://github.com/Kitware/CMake/blob/master/Modules/Compiler/GNU-Fortran.cmake
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
  if(NOT APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcmodel=medium")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mpreprocess")
  if(NOT APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcmodel=medium")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpconstant -fpp")
  if(NOT APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium -shared-intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcmodel=medium -shared-intel")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL XL)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsuffix=cpp=f")
  if(NOT APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -q64")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcmodel=medium -q64")
  endif()
else()
  message(FATAL_ERROR "Specified compiler (${CMAKE_Fortran_COMPILER_ID}) is not supported by Nek5000")
endif()

# =============================================================================
# Targets
# =============================================================================

include(${CMAKE_CURRENT_LIST_DIR}/genbox/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/genmap/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/reatore2/CMakeLists.txt)
