enable_testing()

# =============================================================================
# Sources, Targets
# =============================================================================

function(add_gslib_test target)
  # TODO: This appears to rebuild the object files. Doesn't seem to handle linking to gslib properly.
  add_executable(${target}  EXCLUDE_FROM_ALL ${CMAKE_CURRENT_LIST_DIR}/${target}.c)
  target_link_libraries(${target} PUBLIC gs)
  add_test(NAME run_${target}
      COMMAND mpirun -np 2 ./${target}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

add_gslib_test(comm_test)
add_gslib_test(crs_test)
add_gslib_test(crystal_test)
add_gslib_test(findpts_el_2_test)
add_gslib_test(findpts_el_2_test2)
add_gslib_test(findpts_el_3_test)
add_gslib_test(findpts_el_3_test2)
add_gslib_test(findpts_local_test)
add_gslib_test(findpts_test)
add_gslib_test(gs_test)
add_gslib_test(gs_unique_test)
add_gslib_test(lob_bnd_test)
add_gslib_test(obbox_test)
add_gslib_test(poly_test)
add_gslib_test(poly_test2)
add_gslib_test(sarray_sort_test)
add_gslib_test(sarray_transfer_test)
add_gslib_test(sort_test)
add_gslib_test(sort_test2)
add_gslib_test(spchol_test)
add_gslib_test(xxt_test)
add_gslib_test(xxt_test2)
