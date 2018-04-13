#===============================================================================
# CMakeLists for STREAM testing executables
#
# This builds the CPP executables for running the STREAM tests.  It is not
# a testing harness itself; it simply builds the executables for use by
# an external testing harness.
#
# This CMakeLists does not build libopenmc nor libnek5000.  It assumes they
# have were beforehand.  Hence, the external testing harness must ensure that the
# libs are built before the test exectutables.
#===============================================================================

#find_library(LIBOPENMC_FOUND libopenmc PATHS ${CMAKE_CURRENT_BINARY_DIR} NO_DEFAULT_PATH)
#if(NOT LIBOPENMC_FOUND)
#  message(FATAL_ERROR "Could not find libopenmc")
#endif()

#find_library(LIBNEK5000_FOUND libnek5000 PATHS ${CMAKE_CURRENT_BINARY_DIR} NO_DEFAULT_PATH)
#if(NOT LIBNEK5000_FOUND)
#  message(FATAL_ERROR "Could not find libnek5000")
#endif()

#add_subdirectory(src/nek5000_api)
add_executable(test_openmc_api src/openmc_api/main.cpp)
