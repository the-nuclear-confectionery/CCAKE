###############################################################################
add_library(src OBJECT
	BSQHydro.cpp
	ccake.cpp
	eos.cpp
	eos_extension.cpp
	eos_table.cpp
	formatted_output.cpp
	input_output.cpp
	kernel.cpp
	linklist.cpp
	matrix.cpp
	particle.cpp
	rootfinder.cpp
	sph_workstation.cpp
	system_state.cpp
	transport_coefficients.cpp
	vector.cpp
)


###############################################################################
# link to libraries and include directories
target_link_libraries(src PUBLIC ${GSL_LIBRARIES})
target_include_directories(src PUBLIC ${GSL_INCLUDE_DIR})

if (HDF5_FOUND)
  link_directories( ${HDF5_LIBRARY_DIRS} )
  include_directories( ${HDF5_INCLUDE_DIR} )
  target_link_libraries(src PUBLIC ${HDF5_CXX_LIBRARIES})
endif()

target_include_directories(src PUBLIC ../include)
