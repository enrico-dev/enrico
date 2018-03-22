cmake_minimum_required(VERSION 3.1)
project(stream Fortran C CXX)

# ROR: 2018-03-22: The header is placed in-source.
# Other alternative is to place in build directory.
include(FortranCInterface)
FortranCInterface_VERIFY()
FortranCInterface_HEADER(${CMAKE_CURRENT_LIST_DIR}/include/nek_interface.h
    MACRO_NAMESPACE C2F_
    SYMBOL_NAMESPACE C2F_
    SYMBOLS
    nek_init
    nek_init_step
    nek_step
    nek_finalize_step
    nek_end)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/nek5000)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/openmc)

add_executable(stream
    ${CMAKE_CURRENT_LIST_DIR}/src/drivers.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/nek_interface.f)
target_include_directories(stream PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include
    ${CMAKE_CURRENT_LIST_DIR}/src/openmc/include
    ${CMAKE_BINARY_DIR})
target_link_libraries(stream nek)
