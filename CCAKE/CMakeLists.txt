###############################################################
cmake_minimum_required(VERSION 3.4)
project(ccake VERSION 1.0 LANGUAGES CXX C)

configure_file(include/config.h.in config.h)
include_directories(${PROJECT_BINARY_DIR})

###############################################################
# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

###############################################################
# include useful compiler/linker flags
set(GCC_COVERAGE_COMPILE_FLAGS "-Ofast -s -march=native -mtune=native")
set(GCC_COVERAGE_LINK_FLAGS    "-lstdc++")
set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")


###############################################################
#list(APPEND CMAKE_PREFIX_PATH "/projects/jnorhos/local")

add_executable(${PROJECT_NAME} "")

add_subdirectory(src)

target_link_libraries(${PROJECT_NAME} PRIVATE src)

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_SOURCE_DIR})
