include(${PROJECT_SOURCE_DIR}/cmake/add_test.cmake)

# glymphbie_add_test(./Annular_tests.cpp)
glymphbie_add_test_mpi(./Annular_tests.cpp 2)
glymphbie_add_test(./FuncWall_tests.cpp)
glymphbie_add_test(./IntegratorWall_tests.cpp)
# glymphbie_add_test_mpi(./Annular_mpi_tests.cpp 2)


