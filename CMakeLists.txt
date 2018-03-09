cmake_minimum_required(VERSION 3.0)
project(stream CXX)

include_directories(stream ${PROJECT_SOURCE_DIR}/include/ /usr/include/mpich)
add_executable(stream
        ${PROJECT_SOURCE_DIR}/src/main.cpp
        ${PROJECT_SOURCE_DIR}/src/drivers.cpp)
