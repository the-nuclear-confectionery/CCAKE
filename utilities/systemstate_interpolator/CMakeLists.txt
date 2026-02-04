cmake_minimum_required(VERSION 3.16)
project(systemstate_interpolator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -fopenmp")
find_package(Kokkos REQUIRED)
find_package(Cabana REQUIRED)

# Root of the main project (two levels up from utilities/systemstate_interpolator)
set(PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")

add_executable(system_state_sph
  main.cpp
  system_state_int.cpp
  ${PROJECT_ROOT}/src/kernel.cpp
)

target_include_directories(system_state_sph PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}          # for system_state.h (and any local headers)
  ${PROJECT_ROOT}/include            # for kernel.h
)

# If kernel.h uses Kokkos macros/types, you likely need Kokkos too.
# Uncomment and adapt if your main project already provides Kokkos targets.
# find_package(Kokkos REQUIRED)
# target_link_libraries(system_state_sph PRIVATE Kokkos::kokkos)
