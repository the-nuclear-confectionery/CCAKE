cmake_minimum_required(VERSION 3.4)
project(persephone VERSION 1.0 LANGUAGES CXX)

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "/projects/jnorhos/local")

# set up GSL libraries
#find_package(GSL)
#target_include_directories(${PROJECT_NAME} PUBLIC ${GSL_INCLUDE_DIR})
#target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})


# default build type: Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# default install prefix: ~/.local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local"
      CACHE PATH "Install path prefix, prepended onto install directories."
      FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

option(MORE_WARNINGS "enable more compiler warnings" ON)
if(MORE_WARNINGS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # adapted from http://stackoverflow.com/a/9862800
  string(APPEND CMAKE_CXX_FLAGS " -Wpedantic -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-conversion -Wsign-promo -Wstrict-overflow=2 -Wundef")
endif()

string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Werror")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Og")
endif()

#set(LIBRARY_NAME "lib${PROJECT_NAME}")
#add_library(${LIBRARY_NAME} "")

add_executable(${PROJECT_NAME} "")

add_subdirectory(src)

target_link_libraries(${PROJECT_NAME} PRIVATE src)


