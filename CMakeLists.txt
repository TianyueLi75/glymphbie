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

# Check that the compiler supports the C++ standard
set(GLYMPHBIE_MIN_GCC_VERSION 8.0)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${GLYMPHBIE_MIN_GCC_VERSION})
    message(FATAL_ERROR "GCC version must be at least ${GLYMPHBIE_MIN_GCC_VERSION}")
  endif()
else()
  message(WARNING 
          "You are using an unsupported compiler! "
          "Please use GCC >= ${GLYMPHBIE_MIN_GCC_VERSION}")
endif()

# -------------------------------
#  User options
# -------------------------------
option(DEBUG "Enable debug build with sanitizers" OFF)
option(ENABLE_BLAS "Enable BLAS support" OFF)
option(ENABLE_LAPACK "Enable LAPACK support" OFF)
option(ENABLE_FFTW "Enable FFTW (double) support" OFF)
option(ENABLE_FFTWF "Enable FFTW (float) support" OFF)
option(ENABLE_FFTWL "Enable FFTW (long double) support" OFF)

# Paths
set(PVFMM_DIR        "${CMAKE_SOURCE_DIR}/extern/pvfmm" CACHE PATH "")
set(SCTL_INCLUDE_DIR "${PVFMM_DIR}/SCTL/include"         CACHE PATH "")
set(CSBQ_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/extern/CSBQ/include" CACHE PATH "")
set(SCTL_DATA_PATH   "${CMAKE_SOURCE_DIR}/extern/CSBQ/data"    CACHE PATH "")

set(PVFMM_INC_DIR "${PVFMM_DIR}/include")
set(PVFMM_LIB_DIR "${PVFMM_DIR}/lib/.libs")

# -------------------------------
# SCTL / PVFMM Base Defines
# -------------------------------
add_compile_definitions(
    SCTL_PROFILE=5
    SCTL_VERBOSE
    SCTL_HAVE_MPI
    SCTL_HAVE_PVFMM
    SCTL_DATA_PATH="${SCTL_DATA_PATH}"
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
GLYMPHBIE_add_cxx_compiler_flag("-Wall")
GLYMPHBIE_add_cxx_compiler_flag("-Wextra")
GLYMPHBIE_add_cxx_compiler_flag("-Wcast-align")
GLYMPHBIE_add_cxx_compiler_flag("-Wcast-qual")
GLYMPHBIE_add_cxx_compiler_flag("-Wdouble-promotion")
GLYMPHBIE_add_cxx_compiler_flag("-Wfloat-equal")
GLYMPHBIE_add_cxx_compiler_flag("-Wimplicit-fallthrough")
GLYMPHBIE_add_cxx_compiler_flag("-Wnon-virtual-dtor")
GLYMPHBIE_add_cxx_compiler_flag("-Wshadow")
GLYMPHBIE_add_cxx_compiler_flag("-Wunused")
GLYMPHBIE_add_cxx_compiler_flag("-Wvla")
GLYMPHBIE_add_cxx_compiler_flag("-Wzero-as-null-pointer-constant")
GLYMPHBIE_add_cxx_compiler_flag("-Wunreachable-code")
GLYMPHBIE_add_cxx_compiler_flag("-Wmissing-noreturn")
GLYMPHBIE_add_cxx_compiler_flag("-Wtautological-compare")
GLYMPHBIE_add_cxx_compiler_flag("-Wundef")
GLYMPHBIE_add_cxx_compiler_flag("-Werror")
GLYMPHBIE_add_cxx_compiler_flag("-Wold-style-cast")
GLYMPHBIE_add_cxx_compiler_flag("-pedantic")

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

# ------------------------------------------------------
# BLAS
# ------------------------------------------------------
if(ENABLE_BLAS)
    find_package(BLAS REQUIRED)

    if(BLAS_FOUND)
        message(STATUS "BLAS found: ${BLAS_LIBRARIES}")
        add_compile_definitions(SCTL_HAVE_BLAS)
        list(APPEND EXTRA_LIBS ${BLAS_LIBRARIES})
    else()
        message(FATAL_ERROR "BLAS requested but not found")
    endif()
endif()

# ------------------------------------------------------
# LAPACK
# ------------------------------------------------------
if(ENABLE_LAPACK)
    find_package(LAPACK REQUIRED)

    if(LAPACK_FOUND)
        message(STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
        add_compile_definitions(SCTL_HAVE_LAPACK)
        list(APPEND EXTRA_LIBS ${LAPACK_LIBRARIES})
    else()
        message(FATAL_ERROR "LAPACK requested but not found")
    endif()
endif()

# ------------------------------------------------------
# FFTW (double)
# ------------------------------------------------------
if(ENABLE_FFTW)
    find_package(FFTW3 REQUIRED COMPONENTS DOUBLE)

    if(FFTW3_FOUND)
        message(STATUS "Found FFTW (double): ${FFTW3_LIBRARIES}")
        add_compile_definitions(SCTL_HAVE_FFTW)
        list(APPEND EXTRA_LIBS ${FFTW3_LIBRARIES})
    else()
        message(FATAL_ERROR "FFTW (double) requested but not found")
    endif()
endif()

# ------------------------------------------------------
# FFTWF (float)
# ------------------------------------------------------
if(ENABLE_FFTWF)
    find_package(FFTW3f REQUIRED)

    if(FFTW3F_FOUND)
        message(STATUS "Found FFTWF (float): ${FFTW3F_LIBRARIES}")
        add_compile_definitions(SCTL_HAVE_FFTWF)
        list(APPEND EXTRA_LIBS ${FFTW3F_LIBRARIES})
    else()
        message(FATAL_ERROR "FFTWf (float) requested but not found")
    endif()
endif()

# ------------------------------------------------------
# FFTWL (long double)
# ------------------------------------------------------
if(ENABLE_FFTWL)
    find_package(FFTW3l REQUIRED)

    if(FFTW3L_FOUND)
        message(STATUS "Found FFTWl (long double): ${FFTW3L_LIBRARIES}")
        add_compile_definitions(SCTL_HAVE_FFTWL)
        list(APPEND EXTRA_LIBS ${FFTW3L_LIBRARIES})
    else()
        message(FATAL_ERROR "FFTWl (long double) requested but not found")
    endif()
endif()

#===============================================================================
# Create the library
#===============================================================================

# Source files
# set(GLYMPHBIE_SOURCES
#   "src/Annular.cpp"
#   "src/Particle.cpp"
#   "src/SpatialSolver.cpp"
#   "src/Wall.cpp"
# )
set(GLYMPHBIE_SOURCES
    "src/Annular.cpp"
)

add_library(GLYMPHBIE SHARED ${GLYMPHBIE_SOURCES})

target_include_directories(GLYMPHBIE
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# -------------------------------
#  Include + Library paths
# -------------------------------
target_include_directories(GLYMPHBIE PUBLIC
    ${SCTL_INCLUDE_DIR}
    ${CSBQ_INCLUDE_DIR}
    ${PVFMM_INC_DIR}
)

target_link_directories(GLYMPHBIE PRIVATE ${PVFMM_LIB_DIR})
target_link_libraries(GLYMPHBIE PRIVATE pvfmm)

# -------------------------------
#  Source discovery, build tests
# -------------------------------
file(GLOB TEST_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/tests/*.cpp"
)

foreach(src ${TEST_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})

    target_link_libraries(${name} PRIVATE GLYMPHBIE
        ${EXTRA_LIBS}
    )

    set_target_properties(${name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endforeach()