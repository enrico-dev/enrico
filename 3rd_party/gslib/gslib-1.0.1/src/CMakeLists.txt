# =============================================================================
# Targets and Sources
# =============================================================================

add_library(gs "")

set(gs_sources
    ${CMAKE_CURRENT_LIST_DIR}/comm.c
    ${CMAKE_CURRENT_LIST_DIR}/crystal.c
    ${CMAKE_CURRENT_LIST_DIR}/fail.c
    ${CMAKE_CURRENT_LIST_DIR}/gs.c
    ${CMAKE_CURRENT_LIST_DIR}/gs_local.c
    ${CMAKE_CURRENT_LIST_DIR}/sarray_sort.c
    ${CMAKE_CURRENT_LIST_DIR}/sarray_transfer.c
    ${CMAKE_CURRENT_LIST_DIR}/sort.c
    ${CMAKE_CURRENT_LIST_DIR}/tensor.c)
set(xxt_sources
    ${CMAKE_CURRENT_LIST_DIR}/sparse_cholesky.c
    ${CMAKE_CURRENT_LIST_DIR}/xxt.c)
set(amg_sources
    ${CMAKE_CURRENT_LIST_DIR}/amg.c)
set(fwrapper_sources
    ${CMAKE_CURRENT_LIST_DIR}/fcrs.c
    ${CMAKE_CURRENT_LIST_DIR}/fcrystal.c
    ${CMAKE_CURRENT_LIST_DIR}/findpts.c)
set(intp_sources
    ${CMAKE_CURRENT_LIST_DIR}/findpts_el_2.c
    ${CMAKE_CURRENT_LIST_DIR}/findpts_el_3.c
    ${CMAKE_CURRENT_LIST_DIR}/findpts_local.c
    ${CMAKE_CURRENT_LIST_DIR}/lob_bnd.c
    ${CMAKE_CURRENT_LIST_DIR}/obbox.c
    ${CMAKE_CURRENT_LIST_DIR}/poly.c)

target_sources(gs PUBLIC
    ${gs_sources}
    ${xxt_sources}
    ${amg_sources}
    ${fwrapper_sources}
    ${intp_sources}
    ${CMAKE_CURRENT_LIST_DIR}/rand_elt_test.c)
