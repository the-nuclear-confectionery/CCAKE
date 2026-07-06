###############################################################
cmake_minimum_required(VERSION 3.5)
project(ccake VERSION 2.0 LANGUAGES CXX C)

configure_file(include/config.h.in config.h)
#include_directories(${PROJECT_BINARY_DIR} /usr/include)

###############################################################
# require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

###############################################################
# Default to an optimized Release build when the user does not specify one.
# (Must come before the Debug/optimization checks below so they see the value.)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

###############################################################
# Setup debug build
if(CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "Configuring debug build.")
  set(CMAKE_DEBUG_POSTFIX d)
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -Wall -Wextra -Wpedantic -Werror") # <- Beware of -Werror when using this line
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  add_compile_definitions(DEBUG)
endif()

option(USE_CUDA "Enable compilation for NVIDIA GPUs" OFF)
if(USE_CUDA)
  enable_language(CUDA)
    set(CMAKE_CUDA_SEPARABLE_COMPILATION "ON")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --Xcudafe --display_error_number=0")
endif()

###############################################################
# Optimization flags (CPU host build)
# CCAKE_ARCH selects the target CPU for -march/-mtune. Default "native" tunes for
# the machine you compile ON. On Delta both login and cpu nodes are AMD Milan, so
# "native" resolves to znver3; set -DCCAKE_ARCH=znver3 (or x86-64-v3) to build a
# portable binary elsewhere.
set(CCAKE_ARCH "native" CACHE STRING "CPU architecture passed to -march/-mtune")

option(CCAKE_LTO "Enable link-time optimization (slower build/link)" OFF)

if(NOT USE_CUDA AND NOT CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "Enabling optimized flags (-O3, arch '${CCAKE_ARCH}', unroll).")
  # -O3 + arch tuning + loop unrolling on top of the build-type defaults (Release
  # already adds -O3 -DNDEBUG). -ffast-math / -Ofast are DELIBERATELY avoided: the
  # SPH root-finder, EoS inversion and causality checks depend on strict IEEE
  # semantics (NaN/inf handling, no reassociation), which -ffast-math would break.
  # -fno-math-errno is safe (we never read errno after libm calls) and lets the
  # compiler inline sqrt/etc.
  add_compile_options(-O3 -march=${CCAKE_ARCH} -mtune=${CCAKE_ARCH} -funroll-loops -fno-math-errno)
  # LTO is opt-in: across Kokkos/Cabana templates it can greatly increase link
  # time. Enable with -DCCAKE_LTO=ON once a plain build is confirmed working.
  if(CCAKE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _ccake_ipo_ok OUTPUT _ccake_ipo_msg)
    if(_ccake_ipo_ok)
      message(STATUS "LTO enabled.")
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    else()
      message(WARNING "CCAKE_LTO=ON but IPO unsupported: ${_ccake_ipo_msg}")
    endif()
  endif()
endif()

set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lstdc++")


###############################################################
# Find dependencies
find_package(GSL COMPONENTS C CXX REQUIRED)
find_package(HDF5 COMPONENTS CXX REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(Kokkos REQUIRED)

find_package(Cabana REQUIRED)

###############################################################
# add include directories
include_directories(include)
include_directories(${GSL_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(build)
###############################################################
# Hack to supress the use of -isystem flags
set(CMAKE_NO_SYSTEM_FROM_IMPORTED 1)

###############################################################
# Create a script to clean cmake
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/clean-cmake.sh "rm -rf CMakeCache.txt  CMakeFiles Makefile  ccake *.cmake  config.h src *_test")
# Set the permissions to execute for the owner
execute_process(COMMAND chmod 755 ${CMAKE_CURRENT_BINARY_DIR}/clean-cmake.sh)
add_custom_target(clean-cmake
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/clean-cmake.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Cleaning cmake files"
)

###############################################################
#list(APPEND CMAKE_PREFIX_PATH "/projects/jnorhos/local")
add_subdirectory(src)
add_executable(${PROJECT_NAME} src/ccake.cpp)
add_executable(generate_eos_table utilities/generate_eos_table.cpp)
add_executable(validate_interpolator utilities/validate_interpolator.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE src)
target_link_libraries(generate_eos_table PRIVATE eos)
target_link_libraries(validate_interpolator PRIVATE src)

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_SOURCE_DIR})

###############################################################
#Turn on testing with 'cmake -DTESTING=ON ..'
option(TESTING "Enable testing" OFF)
if(ENABLE_TESTING)
  # fetch googletest for unit testing
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/b796f7d44681514f58a683a3a71ff17c94edb0c1.zip
   )
   # For Windows: Prevent overriding the parent project's compiler/linker settings
   set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
   FetchContent_MakeAvailable(googletest)
   enable_testing()
   add_executable(vector_test unit_tests/vector_test.cpp)
   add_executable(input_test unit_tests/input_test.cpp)
   add_executable(kernel_test unit_tests/kernel_test.cpp)
   add_executable(smoothing_test unit_tests/smoothing_test.cpp)
   add_executable(grad_smoothing_test unit_tests/grad_smoothing_test.cpp)
   add_executable(eom_default_test unit_tests/eom_default_test.cpp)
   add_executable(evolver_test unit_tests/evolver_test.cpp)
   add_executable(eos_interpolator_test unit_tests/eos_interpolator_test.cpp)
   target_link_libraries(vector_test GTest::gtest_main)
   target_link_libraries(input_test GTest::gtest_main src ${YAML_CPP_LIBRARIES})
   target_link_libraries(kernel_test GTest::gtest_main src ${YAML_CPP_LIBRARIES})
   target_link_libraries(smoothing_test GTest::gtest_main src ${YAML_CPP_LIBRARIES})
   target_link_libraries(grad_smoothing_test GTest::gtest_main src ${YAML_CPP_LIBRARIES})
   target_link_libraries(eom_default_test GTest::gtest_main src ${YAML_CPP_LIBRARIES})
   target_link_libraries(evolver_test GTest::gtest_main src ${YAML_CPP_LIBRARIES})
   target_link_libraries(eos_interpolator_test GTest::gtest_main src ${YAML_CPP_LIBRARIES})

   include(GoogleTest)
   gtest_discover_tests(vector_test)
   gtest_discover_tests(input_test)
   gtest_discover_tests(kernel_test)
   gtest_discover_tests(grad_smoothing_test)
   gtest_discover_tests(eom_default_test)
   gtest_discover_tests(evolver_test)
   gtest_discover_tests(eos_interpolator_test)
endif()
