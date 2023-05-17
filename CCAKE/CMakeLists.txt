###############################################################
cmake_minimum_required(VERSION 3.4)
project(ccake VERSION 1.0 LANGUAGES CXX C)

configure_file(include/config.h.in config.h)
#include_directories(${PROJECT_BINARY_DIR} /usr/include)

###############################################################
# require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

###############################################################
# Setup debug build
if(CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "Configuring debug build.")
  set(CMAKE_DEBUG_POSTFIX d)
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -Wall -Wextra -Wpedantic -Werror") # <- Beware of -Werror when using this line
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -Wall")
  add_compile_definitions(DEBUG)
endif()

###############################################################
# include useful compiler/linker flags
if (CMAKE_BUILD_TYPE MATCHES Release )
  set(GCC_COVERAGE_COMPILE_FLAGS "-Ofast -s -march=native -mtune=native")
endif()
#set(GCC_COVERAGE_LINK_FLAGS    "-lstdc++")
set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

###############################################################
# set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

###############################################################
# add include directories
include_directories(include)

###############################################################
# Hack to supress the use of -isystem flags
set(CMAKE_NO_SYSTEM_FROM_IMPORTED 1)

###############################################################
#list(APPEND CMAKE_PREFIX_PATH "/projects/jnorhos/local")
add_executable(${PROJECT_NAME} src/ccake.cpp)

add_subdirectory(src)

target_link_libraries(${PROJECT_NAME} PRIVATE src)

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_SOURCE_DIR})

###############################################################
#Turn on testing with 'cmake -DTESTING=ON ..'
option(TESTING "Enable testing" OFF)
if(TESTING)
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
   target_link_libraries(vector_test GTest::gtest_main)
   target_link_libraries(input_test GTest::gtest_main src ${YAML_CPP_LIBRARIES})

   include(GoogleTest)
   gtest_discover_tests(vector_test)
   gtest_discover_tests(input_test)
endif()

###############################################################
# Generate helper script to clean cmake
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/clean-cmake.sh
  COMMAND echo '#!/bin/bash
  rm -r CMakeCache.txt  CMakeFiles Makefile  ccake cmake_install.cmake  config.h src' > ${CMAKE_CURRENT_BINARY_DIR}/clean-cmake.sh
)