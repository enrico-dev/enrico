cmake_minimum_required(VERSION 3.1)
project(stream Fortran C CXX)

# Produces a C header file (name_mangling.h) to declare Nek5000's Fortran routines
include(FortranCInterface)
FortranCInterface_VERIFY()
FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/nek_mangling.h
    MACRO_NAMESPACE C2F_
    SYMBOL_NAMESPACE C2F_
    SYMBOLS
    nek_init
    nek_init_step
    nek_step
    nek_finalize_step
    nek_end)

# Headers for all targets
include_directories(
    ${CMAKE_CURRENT_LIST_DIR}/include
    ${CMAKE_CURRENT_LIST_DIR}/src/openmc/include
    ${CMAKE_CURRENT_LIST_DIR}/src/nek5000/core
    ${CMAKE_CURRENT_LIST_DIR}/src/nek5000/3rd_party/gslib/gslib-1.0.1/src
    ${CMAKE_BINARY_DIR})

# For building Nek5000. If there is no input deck, put dummies there.
option(CASENAME    "Name of problem case (.usr and .rea files)")
if(NOT CASENAME)
  message(WARNING "A Nek5000 case was *NOT* provided with '-DCASENAME=<casename>'.  Nek5000 will be configured with a generic input deck.'")
  set(CASENAME zero)
  file(COPY ${CMAKE_CURRENT_LIST_DIR}/config/zero.usr DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  find_file(FOUND_SIZE ${CMAKE_CURRENT_BINARY_DIR}/SIZE)
  if(NOT FOUND_SIZE)
    file(COPY ${CMAKE_CURRENT_LIST_DIR}/config/SIZE DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endif()
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/nek5000)

# For building OpenMC
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/openmc)

# Target for STREAM library
add_library(libstream ${CMAKE_CURRENT_LIST_DIR}/src/drivers.cpp )
set_target_properties(libstream PROPERTIES OUTPUT_NAME stream)
target_link_libraries(libstream PUBLIC libnek5000 libopenmc)

# Target for STREAM executable
add_executable(stream ${CMAKE_CURRENT_LIST_DIR}/src/main.cpp)
target_link_libraries(stream PUBLIC libstream)

# For building tests
add_executable(test_openmc_api ${CMAKE_CURRENT_LIST_DIR}/tests/openmc_api/main.cpp)
target_link_libraries(test_openmc_api PUBLIC libstream)
