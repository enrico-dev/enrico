cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(gslib C Fortran)

enable_testing()

# =============================================================================
# Sources and Targets
# =============================================================================

# Adds gslib sources/target
set(gs_sources
    src/comm.c
    src/crystal.c
    src/fail.c
    src/gs.c
    src/gs_local.c
    src/sarray_sort.c
    src/sarray_transfer.c
    src/sort.c
    src/tensor.c)
set(xxt_sources
    src/sparse_cholesky.c
    src/xxt.c)
set(amg_sources
    src/amg.c)
set(fwrapper_sources
    src/fcrs.c
    src/fcrystal.c
    src/findpts.c)
set(intp_sources
    src/findpts_el_2.c
    src/findpts_el_3.c
    src/findpts_local.c
    src/lob_bnd.c
    src/obbox.c
    src/poly.c)

add_library(libgs
    ${gs_sources}
    ${xxt_sources}
    ${amg_sources}
    ${fwrapper_sources}
    ${intp_sources}
    src/rand_elt_test.c)
set_target_properties(libgs PROPERTIES OUTPUT_NAME gs)

# =============================================================================
# Headers
# =============================================================================

target_include_directories(libgs PUBLIC src)

# =============================================================================
# Preprocessor Symbols
# =============================================================================

if(MPI OR NOT DEFINED MPI)
  set(MPI TRUE)
  target_compile_definitions(libgs PRIVATE -DMPI)
endif()

if(MPI AND NOT NOMPIIO)
  target_compile_definitions(libgs PRIVATE -DUSEMPIIO)
endif()

include(FortranCInterface)
FortranCInterface_VERIFY()
if(FortranCInterface_GLOBAL_SUFFIX STREQUAL "_")
  target_compile_definitions(libgs PRIVATE -DUNDERSCORE)
endif()
# set(CMAKE_Fortran_FORMAT "FIXED")

if(BLAS)
  target_compile_definitions(libgs PRIVATE -DUSE_CBLAS)
else()
  target_compile_definitions(libgs PRIVATE -DUSE_NAIVE_BLAS)
endif()

if(USREXIT)
  target_compile_definitions(libgs PRIVATE -DUSE_USR_EXIT)
endif()

if(NOT DEFINED PREFIX)
  set(PREFIX gslib_)
endif()
target_compile_definitions(libgs PRIVATE -DPREFIX=${PREFIX})

if(NOT DEFINED FPREFIX)
  set(FPREFIX fgslib_)
endif()
target_compile_definitions(libgs PRIVATE -DFPREFIX=${FPREFIX})

target_compile_definitions(libgs PRIVATE -DGLOBAL_LONG_LONG)

# =============================================================================

# For tests
# include(tests/CMakeLists.txt)
