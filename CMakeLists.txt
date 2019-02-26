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
# Heat transfer surrogate
# =============================================================================
add_library(heat_xfer vendor/surrogates/heat_xfer_backend.cpp)
target_include_directories(heat_xfer PUBLIC vendor/surrogates/)

# =============================================================================
# Build libstream
# =============================================================================
# Build all compoments of libstream
add_library(libstream
    src/base_drivers.cpp
    src/heat_driver.cpp
    src/message_passing.cpp
    src/openmc_driver.cpp
    src/openmc_interface.cpp
    src/openmc_nek_driver.cpp
    src/nek_driver.cpp
    src/vtk_viz.cpp
    src/openmc_heat_driver.cpp)
set_target_properties(libstream PROPERTIES OUTPUT_NAME stream)
target_link_libraries(libstream PUBLIC libnek5000 libopenmc xtensor heat_xfer pugixml gsl-lite)

# =============================================================================
# Build STREAM driver
# =============================================================================
add_executable(stream src/main.cpp)
target_link_libraries(stream PUBLIC libstream)

# =============================================================================
# Build STREAM tests and demos
# =============================================================================
add_executable(comm_split_demo tests/comm_split_demo/main.cpp)

add_executable(test_nek5000_singlerod tests/singlerod/short/test_nek5000.cpp)
target_link_libraries(test_nek5000_singlerod PUBLIC libstream)

add_executable(test_openmc_singlerod tests/singlerod/short/test_openmc.cpp)
target_link_libraries(test_openmc_singlerod PUBLIC libstream)

# Ensure C++14 standard is used
set_target_properties(
    stream libstream
    comm_split_demo
    test_nek5000_singlerod
    test_openmc_singlerod
    PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF)
