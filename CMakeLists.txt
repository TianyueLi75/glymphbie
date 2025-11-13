cmake_minimum_required(VERSION 3.20)

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

# If no build type is specified, default to Debug 
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, defaulting to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# Set the C++ standard
set(GLYMPHBIE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD ${GLYMPHBIE_CXX_STANDARD} CACHE STRING "Default C++ standard")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Require C++ standard")

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
GLYMPHBIE_add_cxx_compiler_flag("-Wconversion")
GLYMPHBIE_add_cxx_compiler_flag("-Wdouble-promotion")
GLYMPHBIE_add_cxx_compiler_flag("-Wfloat-equal")
GLYMPHBIE_add_cxx_compiler_flag("-Wimplicit-fallthrough")
GLYMPHBIE_add_cxx_compiler_flag("-Wmissing-noreturn")
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

# Add other compiling flags
# TODO

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

#===============================================================================
# Other directories
#===============================================================================

# Enable testing
include(CTest)
enable_testing()
add_subdirectory(tests)

# # Add examples
# add_subdirectory(examples)