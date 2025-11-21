
cmake_minimum_required(VERSION 3.16)

project(GlymphBIE VERSION 0.1.0 LANGUAGES CXX)


#===============================================================================
# Basic CMake configuration
#===============================================================================

# Disable in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR
    "In-source builds are not allowed."
    " Create a separate directory for build files and delete CMakeCache.txt.")
endif()

# Set the C++ standard
set(GLYMPHBIE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Use mpic++ for compiling
set(CMAKE_CXX_COMPILER mpic++)

# -------------------------------
#  User options
# -------------------------------
option(DEBUG "Enable debug build with sanitizers" OFF)
option(ENABLE_BLAS "Enable BLAS support" ON)
option(ENABLE_LAPACK "Enable LAPACK support" ON)
option(ENABLE_FFTW "Enable FFTW (double) support" ON)
option(ENABLE_FFTWF "Enable FFTW (float) support" ON)
option(ENABLE_FFTWL "Enable FFTW (long double) support" OFF)

# Paths
set(PVFMM_DIR        "${CMAKE_SOURCE_DIR}/extern/pvfmm")
set(SCTL_INCLUDE_DIR "${PVFMM_DIR}/SCTL/include")
set(CSBQ_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/extern/csbq/include")
set(SCTL_DATA_PATH   ${CMAKE_SOURCE_DIR}/extern/csbq/data)

set(PVFMM_INC_DIR "${PVFMM_DIR}/include")
set(PVFMM_LIB_DIR "${PVFMM_DIR}/lib/.libs")
set(PVFMM_LIB "${PVFMM_LIB_DIR}/libpvfmm.so")

# Check that PVFMM library exists
if(NOT EXISTS "${PVFMM_LIB}")
    message(FATAL_ERROR "PVFMM library not found at ${PVFMM_LIB}. Build PVFMM first.")
endif()

add_library(pvfmm SHARED IMPORTED)
set_target_properties(pvfmm PROPERTIES
    IMPORTED_LOCATION "${PVFMM_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${PVFMM_INC_DIR};${SCTL_INCLUDE_DIR}"
)

add_library(csbq INTERFACE)


# -------------------------------
# Extra libraries from MakeVariables
# -------------------------------
# FFTW + MKL
set(PVFMM_FFTW_LIBS
    mkl_intel_thread
    mkl_core
    iomp5
    pthread
)
# BLAS + LAPACK
set(PVFMM_BLAS_LAPACK_LIBS
    mkl_intel_lp64
    mkl_sequential
    mkl_core
    pthread
    dl
    stdc++
    gfortran
    m
    quadmath
)
# X11 (if needed)
set(PVFMM_X_LIBS X11)

# Combine all libraries in the order of MakeVariables
set(PVFMM_DEP_LIBS
    ${PVFMM_LIB}            # PVFMM itself
    ${PVFMM_FFTW_LIBS}      # FFTW / MKL
    ${PVFMM_BLAS_LAPACK_LIBS} # BLAS/LAPACK
    ${PVFMM_X_LIBS}         # X11
)

target_link_libraries(pvfmm INTERFACE 
    ${PVFMM_DEP_LIBS}
)

target_compile_definitions(pvfmm INTERFACE
    SCTL_HAVE_PVFMM
    SCTL_HAVE_BLAS
    SCTL_HAVE_LAPACK
    SCTL_HAVE_FFTW
    SCTL_PROFILE=5
    SCTL_VERBOSE
)

target_include_directories(pvfmm INTERFACE
    ${PVFMM_INC_DIR}
    ${SCTL_INCLUDE_DIR}
)
target_include_directories(csbq INTERFACE
    ${CSBQ_INCLUDE_DIR}
)

# -------------------------------
# SCTL / PVFMM Base Defines
# -------------------------------
add_compile_definitions(
    SCTL_HAVE_MPI
    SCTL_DATA_PATH=${SCTL_DATA_PATH}
    SCTL_QUAD_T=__float128
)

#===============================================================================
# Compiler flags
#===============================================================================

# A macro to check if a compiler flag is supported (provided by CMake)
# https://cmake.org/cmake/help/latest/module/CheckCXXCompilerFlag.html
include(CheckCXXCompilerFlag)

# A macro to check if a compiler flag is supported, and if so, add it
# to the C++ compiler flags
macro(GLYMPHBIE_add_cxx_compiler_flag FLAG)
  # SFLAG is the flag without the hyphen, e.g., "-Wall" -> "Wall"
  string(REGEX REPLACE "-" "" SFLAG ${FLAG})
  # Check if the compiler supports the flag
  check_cxx_compiler_flag(${FLAG} COMPILER_SUPPORT_${SFLAG})
  # If the compiler supports the flag, add it to the C++ compiler flags
  if(COMPILER_SUPPORT_${SFLAG})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
  endif()
endmacro()

# Set warning flags
# GLYMPHBIE_add_cxx_compiler_flag("-Wall")
# GLYMPHBIE_add_cxx_compiler_flag("-Wfloat-conversion")

# -------------------------------
# OpenMP
# -------------------------------
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Using OpenMP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # target_link_libraries(GLYMPHBIE PUBLIC OpenMP::OpenMP_CXX)
endif()

if(DEBUG)
    message(STATUS "Debug mode enabled")

    # Optimization & debugging flags
    GLYMPHBIE_add_cxx_compiler_flag("-O2")
    GLYMPHBIE_add_cxx_compiler_flag("-g")
    GLYMPHBIE_add_cxx_compiler_flag("-fno-omit-frame-pointer")

    # Sanitizers (check availability)
    foreach(flag 
        "-fsanitize=address"
        "-fsanitize=leak"
        "-fsanitize=undefined"
        "-fsanitize=pointer-compare"
        "-fsanitize=pointer-subtract"
        "-fsanitize=float-divide-by-zero"
        "-fsanitize=float-cast-overflow"
        "-fno-sanitize-recover=all"
        "-fstack-protector"
    )
        GLYMPHBIE_add_cxx_compiler_flag(${flag})
    endforeach()

    add_compile_definitions(SCTL_MEMDEBUG)
else()
    GLYMPHBIE_add_cxx_compiler_flag("-O3")
    GLYMPHBIE_add_cxx_compiler_flag("-march=native")
    add_compile_definitions(NDEBUG)
endif()

# platform stack-trace flags
GLYMPHBIE_add_cxx_compiler_flag("-rdynamic")
if(APPLE)
    GLYMPHBIE_add_cxx_compiler_flag("-Wl,-no_pie")
else()
    GLYMPHBIE_add_cxx_compiler_flag("-gdwarf-4")
endif()

#===============================================================================
# Create the library
#===============================================================================

# Source files
# set(GLYMPHBIE_SOURCES
#   "src/Annular.cpp"
#   "src/Particle.cpp"
#   "src/StokesBIO.cpp"
#   "src/Wall.cpp"
# )
set(GLYMPHBIE_SOURCES
    "src/Annular.cpp"
    "src/StokesBIO.cpp"
)

add_library(GLYMPHBIE SHARED ${GLYMPHBIE_SOURCES})

target_include_directories(GLYMPHBIE
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(GLYMPHBIE PRIVATE
    csbq
    pvfmm
    ${EXTRA_LIBS}
)

# -------------------------------
#  Source discovery, build tests
# -------------------------------
# Set target test files manually
set(TEST_SOURCES
    "${CMAKE_SOURCE_DIR}/tests/Annular_tests.cpp"
    "${CMAKE_SOURCE_DIR}/tests/StokesBIO_tests.cpp"
    "${CMAKE_SOURCE_DIR}/tests/SpatialSolver_tests.cpp"
)

foreach(src ${TEST_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})

    target_link_libraries(${name} PRIVATE
        GLYMPHBIE
        csbq
        pvfmm
        ${EXTRA_LIBS}
    )

    set_target_properties(${name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endforeach()