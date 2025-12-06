include(${PROJECT_SOURCE_DIR}/cmake/add_test.cmake)

glymphbie_add_test(./Annular_tests.cpp)
glymphbie_add_test_mpi(./Annular_mpi_tests.cpp 2)
# glymphbie_add_test_mpi(./SpatialSolver_tests.cpp 1 1 10 16)
glymphbie_add_test_mpi(./SpatialSolver_tests.cpp 1 1 10 16 0.3 0.48) # Current SpatialSolver_tests in Annular branch still takes R_in and R_out as arguments


