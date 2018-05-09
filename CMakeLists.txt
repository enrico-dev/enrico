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
    nek_end
    nek_solve)

# Headers for all targets
include_directories(
    ${CMAKE_CURRENT_LIST_DIR}/include
    ${CMAKE_CURRENT_LIST_DIR}/src/openmc/include
    ${CMAKE_CURRENT_LIST_DIR}/src/nek5000/core
    ${CMAKE_CURRENT_LIST_DIR}/src/nek5000/3rd_party/gslib/gslib-1.0.1/src
    ${CMAKE_BINARY_DIR})

# For building libnek5000 and libopenmc
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/nek5000)
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

add_executable(test_nek5000_api ${CMAKE_CURRENT_LIST_DIR}/tests/nek5000_api/main.cpp)
target_link_libraries(test_nek5000_api PUBLIC libstream)

add_executable(comm_split_demo ${CMAKE_CURRENT_LIST_DIR}/tests/comm_split_demo/main.cpp)

add_executable(test_coupled_driver ${CMAKE_CURRENT_LIST_DIR}/tests/coupled_driver/main.cpp)
target_link_libraries(test_coupled_driver PUBLIC libstream)
