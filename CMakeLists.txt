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
    include/
    vendor/openmc/include
    vendor/nek5000/core
    vendor/nek5000/3rd_party/gslib/gslib-1.0.1/src
    src/
    ${CMAKE_BINARY_DIR})

# =============================================================================
# Recursively build libnek5000 and libopenmc
# =============================================================================
add_subdirectory(vendor/nek5000)
add_subdirectory(vendor/openmc)
add_subdirectory(vendor/gsl)

# =============================================================================
# Build libstream
# =============================================================================
# Build all compoments of libstream
add_library(libstream
    src/base_drivers.cpp
    src/message_passing.cpp
    src/openmc_driver.cpp
    src/openmc_interface.cpp
    src/openmc_nek_driver.cpp
    src/nek_driver.cpp
    include/base_drivers.h
    include/message_passing.h
    include/nek_interface.h
    include/nek_driver.h
    include/openmc_driver.h
    include/openmc_interface.h
    include/openmc_nek_driver.h
    include/comm.h
    include/stream_const.h
    include/stream_geom.h)
set_target_properties(libstream PROPERTIES OUTPUT_NAME stream)
target_link_libraries(libstream PUBLIC libnek5000 libopenmc gsl-lite)

# =============================================================================
# Build STREAM driver
# =============================================================================
add_executable(stream src/main.cpp)
target_link_libraries(stream PUBLIC libstream)

# =============================================================================
# Build STREAM tests and demos
# =============================================================================
add_executable(test_openmc_api tests/api_tests/test_openmc_api.cpp)
target_link_libraries(test_openmc_api PUBLIC libstream)

add_executable(test_nek5000_api tests/api_tests/test_nek5000_api.cpp)
target_link_libraries(test_nek5000_api PUBLIC libstream)

add_executable(test_coupled_api tests/api_tests/test_coupled_api.cpp)
target_link_libraries(test_coupled_api PUBLIC libstream)

add_executable(comm_split_demo tests/comm_split_demo/main.cpp)

add_executable(test_nek5000_singlerod tests/singlerod/short/test_nek5000.cpp)
target_link_libraries(test_nek5000_singlerod PUBLIC libstream)

add_executable(test_openmc_singlerod tests/singlerod/short/test_openmc.cpp)
target_link_libraries(test_openmc_singlerod PUBLIC libstream)

add_executable(test_coupled_singlerod tests/singlerod/short/test_coupled.cpp)
target_link_libraries(test_coupled_singlerod PUBLIC libstream)
