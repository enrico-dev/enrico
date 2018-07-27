cmake_minimum_required(VERSION 3.1)
project(stream Fortran C CXX)

# =============================================================================
# Discover name-mangling for routines in libnek5000
# =============================================================================
include(FortranCInterface)
FortranCInterface_VERIFY()
FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/nek_mangling.h
    MACRO_NAMESPACE C2F_
    SYMBOL_NAMESPACE C2F_
    SYMBOLS
    nek_init
    nek_end
    nek_solve)

# =============================================================================
# Headers for all targets
# =============================================================================
include_directories(
    ${CMAKE_CURRENT_LIST_DIR}/include
    ${CMAKE_CURRENT_LIST_DIR}/src/openmc/include
    ${CMAKE_CURRENT_LIST_DIR}/src/nek5000/core
    ${CMAKE_CURRENT_LIST_DIR}/src/nek5000/3rd_party/gslib/gslib-1.0.1/src
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface
    ${CMAKE_BINARY_DIR})

# =============================================================================
# Recursively build libnek5000 and libopenmc
# =============================================================================
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/nek5000)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/openmc)

# =============================================================================
# Build libstream
# =============================================================================
# Discover and set necessary flags for building fixed-format Fortran components
# of libstream.
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(FIXEDFORMAT_FLAGS "-fdefault-real-8 -fdefault-double-8 -cpp")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  set(FIXEDFORMAT_FLAGS "-r8 -Mpreprocess")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set(FIXEDFORMAT_FLAGS "-r8 -fpconstant -fpp")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL XL)
  set(FIXEDFORMAT_FLAGS "-qrealsize=8 -qdpc=e -qsuffix=cpp=f -qsuppress=cmpmsg")
else()
  message(FATAL_ERROR "Specified compiler (${CMAKE_Fortran_COMPILER_ID}) is not supported by Nek5000")
endif()
set(FIXEDFORMAT_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface/nek_parallel.F
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface/nek_geom.F
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface/nek_interface_fixedform.F)
set_source_files_properties(${FIXEDFORMAT_SOURCES} PROPERTIES
    Fortran_FORMAT FIXED
    COMPILE_FLAGS ${FIXEDFORMAT_FLAGS})

# Build all compoments of libstream
add_library(libstream
    ${CMAKE_CURRENT_LIST_DIR}/src/drivers.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/message_passing.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface/nek_geom.F
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface/nek_interface.F90
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface/nek_interface_fixedform.F
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface/nek_mass.F
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface/nek_parallel.F
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface/nek_size.F
    ${CMAKE_CURRENT_LIST_DIR}/include/drivers.h
    ${CMAKE_CURRENT_LIST_DIR}/include/message_passing.h
    ${CMAKE_CURRENT_LIST_DIR}/include/nek_interface.h
    ${CMAKE_CURRENT_LIST_DIR}/include/comm.h
    ${CMAKE_CURRENT_LIST_DIR}/include/stream_geom.h)
set_target_properties(libstream PROPERTIES OUTPUT_NAME stream)
target_link_libraries(libstream PUBLIC libnek5000 libopenmc)

# =============================================================================
# Build STREAM driver
# =============================================================================
add_executable(stream ${CMAKE_CURRENT_LIST_DIR}/src/main.cpp)
target_link_libraries(stream PUBLIC libstream)

# =============================================================================
# Build STREAM tests and demos
# =============================================================================
add_executable(test_openmc_api ${CMAKE_CURRENT_LIST_DIR}/tests/api_tests/test_openmc_api.cpp)
target_link_libraries(test_openmc_api PUBLIC libstream)

add_executable(test_nek5000_api ${CMAKE_CURRENT_LIST_DIR}/tests/api_tests/test_nek5000_api.cpp)
target_link_libraries(test_nek5000_api PUBLIC libstream)

add_executable(test_coupled_api ${CMAKE_CURRENT_LIST_DIR}/tests/api_tests/test_coupled_api.cpp)
target_link_libraries(test_coupled_api PUBLIC libstream)

add_executable(comm_split_demo ${CMAKE_CURRENT_LIST_DIR}/tests/comm_split_demo/main.cpp)
