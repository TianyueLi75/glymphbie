macro(glymphbie_add_test FILENAME)
  # Strip the path and extension from the filename to get the test name
  set(TESTNAME ${FILENAME})
  get_filename_component(TESTNAME ${TESTNAME} NAME_WE)
  get_filename_component(TESTNAME ${TESTNAME} NAME_WLE)

  add_executable(${TESTNAME} ${FILENAME})

  target_link_libraries(${TESTNAME} PRIVATE 
    GLYMPHBIE
    csbq
    pvfmm
    ${EXTRA_LIBS}
  ) 

  set_target_properties(${TESTNAME} PROPERTIES 
    CXX_STANDARD ${GLYMPHBIE_CXX_STANDARD}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  )

  add_test(
    NAME ${TESTNAME} 
    COMMAND $<TARGET_FILE:${TESTNAME}>
  )
endmacro()

macro(glymphbie_add_test_mpi FILENAME NPROC)
  # Strip the path and extension from the filename to get the test name
  set(TESTNAME ${FILENAME})
  get_filename_component(TESTNAME ${TESTNAME} NAME_WE)
  get_filename_component(TESTNAME ${TESTNAME} NAME_WLE)
  set(TEST_ARGS ${ARGN})

  add_executable(${TESTNAME} ${FILENAME})

  target_link_libraries(${TESTNAME} PRIVATE 
    GLYMPHBIE
    csbq
    pvfmm
    ${EXTRA_LIBS}
  ) 

  set_target_properties(${TESTNAME} PROPERTIES 
    CXX_STANDARD ${GLYMPHBIE_CXX_STANDARD}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  )

  add_test(
    NAME ${TESTNAME}
    COMMAND mpirun -np ${NPROC} $<TARGET_FILE:${TESTNAME}> ${TEST_ARGS}
  )
endmacro()
