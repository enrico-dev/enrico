if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/genbox.f PROPERTIES COMPILE_FLAGS "-fdefault-real-8")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/genbox.f PROPERTIES COMPILE_FLAGS "-r8")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/genbox.f PROPERTIES COMPILE_FLAGS "-r8")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL XL)
  set_source_files_properties(${CMAKE_CURRENT_LIST_DIR}/genbox.f PROPERTIES COMPILE_FLAGS "-qrealsize=8")
else()
  message(FATAL_ERROR "Specified compiler (${CMAKE_Fortran_COMPILER_ID}) is not supported by Nek5000")
endif()

add_executable(reatore2 EXCLUDE_FROM_ALL
    ${CMAKE_CURRENT_LIST_DIR}/reatore2.f
    ${CMAKE_CURRENT_LIST_DIR}/strings.f
    ${CMAKE_CURRENT_LIST_DIR}/byte.c)

add_executable(rea2torea EXCLUDE_FROM_ALL
    ${CMAKE_CURRENT_LIST_DIR}/re2torea.f
    ${CMAKE_CURRENT_LIST_DIR}/strings.f
    ${CMAKE_CURRENT_LIST_DIR}/byte.c)
