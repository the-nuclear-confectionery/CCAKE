###############################################################################
cmake_minimum_required(VERSION 3.4)
project(ccake VERSION 1.0 LANGUAGES CXX C)

configure_file(include/config.h.in config.h)
include_directories(${PROJECT_BINARY_DIR})

###############################################################################
# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

###############################################################################
#list(APPEND CMAKE_PREFIX_PATH "/projects/jnorhos/local")

add_executable(${PROJECT_NAME} "")



###############################################################################
# find and link GSL libraries
find_package(GSL COMPONENTS C CXX REQUIRED)
target_link_libraries(src PUBLIC ${GSL_LIBRARIES})
target_include_directories(src PUBLIC ${GSL_INCLUDE_DIR})

###############################################################################
# find and link HDF5 libraries
find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND)
	link_directories( ${HDF5_LIBRARY_DIRS} )
	include_directories( ${HDF5_INCLUDE_DIR} )

	target_link_libraries(src PUBLIC ${HDF5_CXX_LIBRARIES})
endif()

###############################################################################
# include useful compiler/linker flags
if (HDF5_FOUND)
  message(STATUS "HDF5 found!")
  add_definitions( -DHDF5 )
else()
  message(STATUS "HDF5 not found!")
endif()

add_compile_options(
     -Wall -Werror -Wno-error=maybe-uninitialized
     -Ofast -s -march=native -mtune=native )
add_link_options( -lstdc++ )



add_subdirectory(src)


###############################################################################
# build and install
target_link_libraries(${PROJECT_NAME} PRIVATE src)

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_SOURCE_DIR})
