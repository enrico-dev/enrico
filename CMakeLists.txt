cmake_minimum_required(VERSION 3.0)
project(stream Fortran C CXX)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/nek5000)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/openmc)
