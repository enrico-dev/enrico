cmake_minimum_required(VERSION 3.3)
project(enrico Fortran C CXX)

option(USE_NEK  "Determines if Nek5000 will be included" ON)

# On BG/Q, linking requires -dynamic, not -rdynamic
if (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "ppc64")
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-dynamic")
endif ()

# =============================================================================
# Check for SCALE (optional)
# =============================================================================
find_package(SCALE PATHS ${SCALE_DIR} COMPONENTS Shift OmnibusDriver Tpetra)

if (SCALE_FOUND)
    MESSAGE(STATUS "SCALE found!")
    MESSAGE(STATUS "SCALE_DIR = ${SCALE_DIR}")
    MESSAGE(STATUS "SCALE_VERSION = ${SCALE_VERSION}")
endif ()

# =============================================================================
# Discover name-mangling for routines in libnek5000
# =============================================================================
if (USE_NEK)
    include(FortranCInterface)
    FortranCInterface_VERIFY()
    FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/nek_mangling.h
        MACRO_NAMESPACE C2F_
        SYMBOL_NAMESPACE C2F_
        SYMBOLS
        nek_init
        nek_end
        nek_solve)
endif()

# =============================================================================
# Headers for all targets
# =============================================================================
include_directories(
    include/
    vendor/iapws/include
    vendor/openmc/include
    vendor
    src/
    ${CMAKE_BINARY_DIR})

if (SCALE_FOUND)
    include_directories(${SCALE_INCLUDE_DIRS})
    include_directories(${SCALE_TPL_INCLUDE_DIRS})
endif ()

# =============================================================================
# Recursively build libnek5000 and libopenmc
# =============================================================================
if (USE_NEK)
    add_subdirectory(vendor/nek5000)
endif()
add_subdirectory(vendor/openmc)

# =============================================================================
# IAPWS Correlations
# =============================================================================
add_library(iapws vendor/iapws/iapws.cpp)
target_link_libraries(iapws PRIVATE gsl-lite)
target_compile_definitions(iapws PRIVATE GSL_THROW_ON_CONTRACT_VIOLATION)

# =============================================================================
# Heat transfer surrogate
# =============================================================================
add_library(heat_xfer vendor/surrogates/heat_xfer_backend.cpp)
target_include_directories(heat_xfer PUBLIC vendor/surrogates/)
target_link_libraries(heat_xfer PUBLIC iapws)

# =============================================================================
# Build libenrico
# =============================================================================

set(SOURCES
    src/driver.cpp
    src/coupled_driver.cpp
    src/comm_split.cpp
    src/surrogate_heat_driver.cpp
    src/mpi_types.cpp
    src/openmc_driver.cpp
    src/cell_instance.cpp
    src/vtk_viz.cpp
    src/heat_fluids_driver.cpp)

if (USE_NEK)
    list(APPEND SOURCES src/nek_driver.cpp)
endif()

if (SCALE_FOUND)
    list(APPEND SOURCES src/shift_driver.cpp)
endif ()

# Build all compoments of libenrico
add_library(libenrico
            ${SOURCES})

set_target_properties(libenrico PROPERTIES OUTPUT_NAME enrico)
set(LIBRARIES
    iapws
    libopenmc
    xtensor
    heat_xfer
    pugixml
    gsl-lite)

target_compile_definitions(libenrico PRIVATE GSL_THROW_ON_CONTRACT_VIOLATION)

if (USE_NEK)
    target_compile_definitions(libenrico PUBLIC USE_NEK)
    list(APPEND LIBRARIES libnek5000)
endif()

if (SCALE_FOUND)
  target_compile_definitions(libenrico PRIVATE USE_SHIFT)
  list(APPEND LIBRARIES ${SCALE_LIBRARIES} ${SCALE_TPL_LIBRARIES})
endif ()

target_link_libraries(libenrico PUBLIC ${LIBRARIES})

# =============================================================================
# Build enrico driver
# =============================================================================
add_executable(enrico src/main.cpp)
target_link_libraries(enrico PUBLIC libenrico)

# =============================================================================
# Build enrico tests and demos
# =============================================================================
add_executable(comm_split_demo tests/comm_split_demo/main.cpp)

add_executable(test_openmc_singlerod tests/singlerod/short/test_openmc.cpp)
target_link_libraries(test_openmc_singlerod PUBLIC libenrico)

# Ensure C++14 standard is used
set_target_properties(
        enrico libenrico
        comm_split_demo
        heat_xfer
        iapws
        test_openmc_singlerod
        PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF)


if (USE_NEK)
    add_executable(test_nek5000_singlerod tests/singlerod/short/test_nek5000.cpp)
    target_link_libraries(test_nek5000_singlerod PUBLIC libenrico)
    set_target_properties(
        test_nek5000_singlerod
        PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF)
endif()

# =============================================================================
# Build unit tests
# =============================================================================
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE vendor/catch/single_include/catch2)

add_executable(unittests
  tests/unit/catch.cpp
  tests/unit/test_surrogate_th.cpp)
target_link_libraries(unittests PUBLIC Catch pugixml libenrico)
set_target_properties(unittests PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF)
