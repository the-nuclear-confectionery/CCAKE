add_library(src OBJECT
	BSQHydro.cpp
	eos.cpp
	eos_interpolator.cpp
	eos_extension.cpp
	eos_table.cpp
	evolver.cpp
	formatted_output.cpp
	kernel.cpp
	matrix.cpp
	particle.cpp
	rootfinder.cpp
	eom_default.cpp
	sph_workstation.cpp
	system_state.cpp
	transport_coefficients.cpp
	vector.cpp
	utilities.cpp
	input.cpp
	output.cpp
)

###############################################################################
# find and link GSL libraries
#set (ENV{CPATH} ${GSL_INCLUDE_DIR})
target_link_libraries(src PUBLIC ${GSL_LIBRARIES})
target_include_directories(src PUBLIC ${GSL_INCLUDE_DIR})

###############################################################################
# find and link HDF5 libraries
link_directories( ${HDF5_LIBRARY_DIRS} )
#include_directories( ${HDF5_INCLUDE_DIR} )
target_link_libraries(src PUBLIC ${HDF5_CXX_LIBRARIES} )

###############################################################################
# Include cabana and kokkos libraries
target_link_libraries(src PUBLIC Kokkos::kokkos)
target_link_libraries(src PUBLIC Cabana::Core)

###############################################################################
# find and link yaml-cpp libraries
find_package(yaml-cpp REQUIRED)
link_directories( ${YAML_CPP_LIBRARY_DIRS} )
#include_directories( ${YAML_CPP_INCLUDE_DIR} )
target_link_libraries(src PUBLIC ${YAML_CPP_LIBRARIES} )
