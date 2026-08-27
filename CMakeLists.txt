cmake_minimum_required(VERSION 3.16)

# -------------------------------
#  MPI / PVFMM optional toolchain
# -------------------------------
# ENABLE_PVFMM builds the periodic BIE library + its examples/tests against PVFMM
# (the far field of every BoundaryIntegralOp). PVFMM has no serial build, so it
# forces ENABLE_MPI on. With both OFF the project builds only the standalone,
# self-contained quad-junctions examples (e.g. the capped-annulus watertightness
# test), which need neither MPI nor PVFMM -- just a serial g++ + MKL + FFTW.
# These are declared before project() so the compiler choice takes effect.
option(ENABLE_MPI   "Build with MPI (mpic++ + SCTL_HAVE_MPI)" ON)
option(ENABLE_PVFMM "Build the PVFMM-backed library, examples and tests (implies ENABLE_MPI)" ON)

if(ENABLE_PVFMM AND NOT ENABLE_MPI)
  message(STATUS "ENABLE_PVFMM=ON requires MPI; forcing ENABLE_MPI=ON")
  set(ENABLE_MPI ON CACHE BOOL "Build with MPI (mpic++ + SCTL_HAVE_MPI)" FORCE)
endif()

# Compiler selection must precede project(): mpic++ for MPI builds, g++ otherwise.
if(ENABLE_MPI)
  set(CMAKE_CXX_COMPILER mpic++)
else()
  if(NOT DEFINED CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER g++)
  endif()
endif()

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

# -------------------------------
#  User options
# -------------------------------
option(DEBUG "Enable debug build with sanitizers" OFF)
option(ENABLE_BLAS "Enable BLAS support" ON)
option(ENABLE_LAPACK "Enable LAPACK support" ON)
option(ENABLE_FFTW "Enable FFTW (double) support" ON)
option(ENABLE_FFTWF "Enable FFTW (float) support" ON)
option(ENABLE_FFTWL "Enable FFTW (long double) support" OFF)
option(ENABLE_GCOV "Enable code coverage report" OFF)

# Paths
message("Setting path to pvfmm and sctl")
# NB: SCTL headers live under the pvfmm submodule (extern/pvfmm/SCTL/include) and
# are needed even for the serial (no-PVFMM) build, so the pvfmm submodule must be
# checked out either way -- only the *built* libpvfmm.so is PVFMM-build specific.
set(PVFMM_DIR        "${CMAKE_SOURCE_DIR}/extern/pvfmm")
set(SCTL_INCLUDE_DIR "${PVFMM_DIR}/SCTL/include")
set(CSBQ_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/extern/csbq/include")
set(QJ_INCLUDE_DIR   "${CMAKE_SOURCE_DIR}/extern/quad-junctions/include")
set(QJ_COMPAT_HEADER "${QJ_INCLUDE_DIR}/quad_junctions/csbq_sctl_compat.hpp")
set(SCTL_DATA_PATH   ${CMAKE_SOURCE_DIR}/extern/csbq/data)

set(PVFMM_INC_DIR "${PVFMM_DIR}/include")
set(PVFMM_LIB_DIR "${PVFMM_DIR}/lib/.libs")
set(PVFMM_LIB "${PVFMM_LIB_DIR}/libpvfmm.so")

add_library(csbq INTERFACE)
target_include_directories(csbq INTERFACE
    ${CSBQ_INCLUDE_DIR}
)

if(ENABLE_PVFMM)
    # Check that PVFMM library exists
    if(NOT EXISTS "${PVFMM_LIB}")
        message(FATAL_ERROR "PVFMM library not found at ${PVFMM_LIB}. Build PVFMM first, or configure with -DENABLE_PVFMM=OFF -DENABLE_MPI=OFF for the serial (quad-junctions) examples only.")
    endif()

    add_library(pvfmm SHARED IMPORTED)
    set_target_properties(pvfmm PROPERTIES
        IMPORTED_LOCATION "${PVFMM_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${PVFMM_INC_DIR};${SCTL_INCLUDE_DIR}"
    )

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
endif()

# -------------------------------
# SCTL Base Defines
# -------------------------------
if(ENABLE_MPI)
    add_compile_definitions(SCTL_HAVE_MPI)
endif()
add_compile_definitions(
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

if(ENABLE_GCOV)
    message(STATUS "Test coverage enabled")
    add_compile_options(--coverage)
    # GLYMPHBIE_add_cxx_compiler_flag("--coverage") # this will remove one hyphen and confuse cxx.
    add_link_options(--coverage)
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

set(GLYMPHBIE_HEADERS
    "include/Annular.hpp"
    "include/PeriodicBIO.hpp"
    "include/FuncWall.hpp"
    "include/detail/FuncWall.tpp"
    "include/Particle.hpp"
    "include/FlowVis.hpp"
    "include/Wall.hpp"
    "include/detail/Wall.tpp"
    "include/IntegratorWall.hpp"
    "include/detail/IntegratorWall.tpp"
)
set(GLYMPHBIE_SOURCES
    "src/Annular.cpp"
    "src/StokesBIO.cpp"
)

# The GLYMPHBIE library, its CTest suite, and the periodic BIE examples all use
# the PVFMM-backed BoundaryIntegralOp, so they are only built when ENABLE_PVFMM.
if(ENABLE_PVFMM)
    add_library(GLYMPHBIE SHARED ${GLYMPHBIE_SOURCES} ${GLYMPHBIE_HEADERS})

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
    # Enable testing
    include(CTest)
    enable_testing()
    add_subdirectory(tests)

    # Include (periodic) examples
    set(EXAMPLE_SOURCES
        "${CMAKE_SOURCE_DIR}/examples/SpatialSolver_tests.cpp" # for convergence results
        "${CMAKE_SOURCE_DIR}/examples/FuncWall_SpatialSolver_ex.cpp" # example with functional walls
        "${CMAKE_SOURCE_DIR}/examples/SpatialSolver_selfconv.cpp"  # Self convergence
    )

    foreach(src ${EXAMPLE_SOURCES})
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
endif()

# -------------------------------
#  Standalone quad-junctions examples (no PVFMM, no MPI required)
# -------------------------------
# These link neither GLYMPHBIE nor pvfmm: they use only the header-only CSBQ +
# quad-junctions surface meshers with MKL (BLAS/LAPACK) + FFTW, so they build in
# a fully serial configuration (-DENABLE_PVFMM=OFF -DENABLE_MPI=OFF).
#
# Include order is LOAD-BEARING: quad-junctions/include must precede
# SCTL_INCLUDE_DIR so the vendored include/sctl/experimental/quad_element.{hpp,cpp}
# shadows upstream SCTL's older ancestor (see extern/quad-junctions/README.md).
# The compat shim is force-included so CSBQ finds SCTL_QUOTEME.
set(STANDALONE_EXAMPLE_SOURCES
    "${CMAKE_SOURCE_DIR}/examples/CappedAnnulus_watertight.cpp"
    "${CMAKE_SOURCE_DIR}/examples/CappedAnnulus_SpatialSolver_tests.cpp"
    "${CMAKE_SOURCE_DIR}/examples/CappedAnnulus_FuncWall_ex.cpp"
    "${CMAKE_SOURCE_DIR}/examples/CappedAnnulus_TracerParticles_verify.cpp"
    "${CMAKE_SOURCE_DIR}/examples/CappedAnnulus_TracerParticles_ex.cpp"
)

# The tracer-particles example writes its whole time series into a single VTKHDF
# (.vtkhdf) container (include/VtkHdfWriter.hpp), which needs the HDF5 C library.
# On HPC, load a serial hdf5 module first (e.g. `module load hdf5`), which sets
# HDF5_ROOT. We resolve the header + library file ourselves rather than using the
# HDF5-exported imported target: the HPC module's CMake package config carries a
# broken link interface, and this path also works for a system libhdf5-dev
# install (CI). Only this one example links HDF5.
find_package(HDF5 QUIET COMPONENTS C)
find_path(GLYMPHBIE_HDF5_INCLUDE hdf5.h
    HINTS $ENV{HDF5_ROOT} ${HDF5_ROOT} ${HDF5_C_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS}
    PATH_SUFFIXES include)
find_library(GLYMPHBIE_HDF5_LIB NAMES hdf5 hdf5_serial
    HINTS $ENV{HDF5_ROOT} ${HDF5_ROOT}
    PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu)
if(NOT GLYMPHBIE_HDF5_INCLUDE OR NOT GLYMPHBIE_HDF5_LIB)
    message(FATAL_ERROR "HDF5 (C) not found for CappedAnnulus_TracerParticles_ex. "
        "Load a serial hdf5 module (`module load hdf5`) or install libhdf5-dev. "
        "Resolved include='${GLYMPHBIE_HDF5_INCLUDE}' lib='${GLYMPHBIE_HDF5_LIB}'.")
endif()
message(STATUS "VTKHDF: using HDF5 include=${GLYMPHBIE_HDF5_INCLUDE} lib=${GLYMPHBIE_HDF5_LIB}")

set(QJ_SERIAL_LIBS
    mkl_intel_lp64 mkl_sequential mkl_core pthread
    fftw3_omp fftw3 fftw3f fftw3l dl
)

foreach(src ${STANDALONE_EXAMPLE_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})

    target_include_directories(${name} BEFORE PRIVATE
        ${QJ_INCLUDE_DIR}
        ${CSBQ_INCLUDE_DIR}
        ${SCTL_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/src
    )
    target_compile_options(${name} PRIVATE
        "-include" "${QJ_COMPAT_HEADER}"
        -mno-avx512fp16
    )
    target_compile_definitions(${name} PRIVATE
        SCTL_HAVE_BLAS SCTL_HAVE_LAPACK
        SCTL_HAVE_FFTW SCTL_HAVE_FFTWF SCTL_HAVE_FFTWL SCTL_FFTW_THREADS
        SCTL_GLOBAL_MEM_BUFF=0
        # Read CSBQ's checked-in precomputed toroidal/special/log quadrature rules
        # (extern/csbq/data) instead of the cwd-relative ./data default, so these
        # binaries never recompute the (slow) singular quadrature tables no matter
        # where they are run from.
        CSBQ_DATA_PATH=${SCTL_DATA_PATH}
    )
    target_link_libraries(${name} PRIVATE ${QJ_SERIAL_LIBS})

    # Only the tracer-particles example uses the VTKHDF writer -> link HDF5 to it.
    if(name STREQUAL "CappedAnnulus_TracerParticles_ex")
        target_include_directories(${name} PRIVATE ${GLYMPHBIE_HDF5_INCLUDE})
        target_link_libraries(${name} PRIVATE ${GLYMPHBIE_HDF5_LIB})
    endif()

    set_target_properties(${name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endforeach()