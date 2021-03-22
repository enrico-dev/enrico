cmake_minimum_required(VERSION 3.3)
project(enrico Fortran C CXX)

#===============================================================================
# RPATH and install information
#===============================================================================

include(GNUInstallDirs)

# NekRS must be installed (via `make install`) to have all headers/libs available.
# For ENRICO, we'll install in the build directory by default
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "" FORCE)
endif ()

# This block of code ensures that dynamic libraries can be found via the RPATH
# whether the executable is the original one from the build directory or the
# installed one in CMAKE_INSTALL_PREFIX. Ref:
# https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)


# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_INSTALL_RPATH
  "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:${CMAKE_INSTALL_PREFIX}/lib")


# =============================================================================
# Nek Options
# =============================================================================

set(NEK_DIST
        "nek5000" CACHE STRING
        "Distribution of Nek5000 or NekRS [nek5000|nekrs|none]")
string(TOLOWER "${NEK_DIST}" NEK_DIST)
if (NOT (NEK_DIST STREQUAL "nek5000" OR
        NEK_DIST STREQUAL "nekrs" OR
        NEK_DIST STREQUAL "none"))
    message(FATAL_ERROR "Invalid value for NEK_DIST.  Must be 'nek5000', 'nekrs', or 'none' (case-insensitive)")
else ()
    set(USE_NEK5000 FALSE)
    set(USE_NEKRS FALSE)
    if (NEK_DIST STREQUAL "none")
        message(STATUS "Compiling without NekRS or Nek5000")
    elseif (NEK_DIST STREQUAL "nek5000")
        set(USE_NEK5000 TRUE)
        message(STATUS "Compiling with Nek5000")
    elseif (NEK_DIST STREQUAL "nekrs")
        set(USE_NEKRS TRUE)
        message(STATUS "Compiling with NekRS")
    endif ()
endif ()

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
if (USE_NEK5000)
    include(FortranCInterface)
    FortranCInterface_VERIFY()
    FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/nek_mangling.h
            MACRO_NAMESPACE C2F_
            SYMBOL_NAMESPACE C2F_
            SYMBOLS
            nek_init
            nek_end
            nek_solve)
endif ()

# =============================================================================
# Verify and create headers for NekRS integration
# =============================================================================

if (USE_NEKRS)
    include(CheckIncludeFileCXX)
    check_include_file_cxx(dlfcn.h DLOPEN_H_FOUND)
    if (NOT ${DLOPEN_H_FOUND})
        message(FATAL_ERROR "CMake could not find dlfcn.h, the header for libdl.")
    endif ()

    set(NEKRS_HOME ${CMAKE_INSTALL_PREFIX})
    configure_file(config/nekrs_home.h.in ${CMAKE_BINARY_DIR}/include/nekrs_home.h @ONLY)
endif ()

# =============================================================================
# Headers for all targets
# =============================================================================
include_directories(
        include/
        vendor/iapws/include
        vendor/openmc/include
        vendor
        src/
        ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}/include)

if (SCALE_FOUND)
    include_directories(${SCALE_INCLUDE_DIRS})
    include_directories(${SCALE_TPL_INCLUDE_DIRS})
endif ()

# =============================================================================
# Recursively build libnek5000/libnekrs and libopenmc
# =============================================================================
if (USE_NEK5000)
    add_subdirectory(vendor/nek5000)
elseif (USE_NEKRS)
    add_subdirectory(vendor/nekRS)
    # TODO: An upstream fix to the NekRS CMakeLists should let us to omit this
    include_directories(vendor/nekRS/src/lib)
endif ()
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

if (USE_NEK5000)
    list(APPEND SOURCES src/nek5000_driver.cpp)
elseif (USE_NEKRS)
    list(APPEND SOURCES src/nekrs_driver.cpp)
endif ()

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

if (USE_NEK5000)
    target_compile_definitions(libenrico PUBLIC USE_NEK5000)
    list(APPEND LIBRARIES libnek5000)
elseif (USE_NEKRS)
    target_compile_definitions(libenrico PUBLIC USE_NEKRS)
    list(APPEND LIBRARIES nekrs-lib)
endif ()

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

add_executable(test_openmc_singlerod tests/singlerod/short/openmc/test_openmc.cpp)
target_link_libraries(test_openmc_singlerod PUBLIC libenrico)

# Ensure C++14 standard is used
set_target_properties(
        enrico libenrico
        comm_split_demo
        heat_xfer
        iapws
        test_openmc_singlerod
        PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF)


if (USE_NEK5000)
    add_executable(test_nek5000_singlerod tests/singlerod/short/nek5000/test_nek5000.cpp)
    target_link_libraries(test_nek5000_singlerod PUBLIC libenrico)
    set_target_properties(
            test_nek5000_singlerod
            PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF)
endif ()

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

#################################################################################
# Install targets
#################################################################################

set(INSTALL_TARGETS 
  enrico 
  libenrico 
  iapws 
  heat_xfer 
  Catch 
  unittests 
  comm_split_demo
  test_openmc_singlerod)

if (USE_NEK5000)
  list(APPEND INSTALL_TARGETS test_nek5000_singlerod)
endif()

install(TARGETS ${INSTALL_TARGETS}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
)
