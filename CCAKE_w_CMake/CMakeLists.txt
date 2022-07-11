###############################################################
cmake_minimum_required(VERSION 3.4)
project(ccake VERSION 1.0 LANGUAGES CXX C)

###############################################################
# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

###############################################################
# set variables to set up installation
# (generic example)
#set(tool_dest "bin")
#set(include_dest "include/my_library-${MY_LIBRARY_VERSION}")
#set(main_lib_dest "lib/my_library-${MY_LIBRARY_VERSION}")
#set(dest_dir "./")

###############################################################
list(APPEND CMAKE_PREFIX_PATH "/projects/jnorhos/local")

add_executable(${PROJECT_NAME} "")

add_subdirectory(src)

target_link_libraries(${PROJECT_NAME} PRIVATE src)

install(TARGETS ccake DESTINATION "./")

