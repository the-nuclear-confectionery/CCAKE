cmake_minimum_required(VERSION 3.4)
project(hello VERSION 1.0 LANGUAGES CXX)

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "/projects/jnorhos/local")
string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra")

add_subdirectory(src)
