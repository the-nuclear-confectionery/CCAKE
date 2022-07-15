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
# find and link GSL libraries
#find_package(GSL COMPONENTS C CXX REQUIRED)
#target_link_libraries(src PUBLIC ${GSL_LIBRARIES})
#target_include_directories(src PUBLIC ${GSL_INCLUDE_DIR})

###############################################################################
# find and link HDF5 libraries
#find_package(HDF5 COMPONENTS CXX)
#if (HDF5_FOUND)
#	set(HDF5_FOUND ${HDF5_FOUND} PARENT_SCOPE)
#	link_directories( ${HDF5_LIBRARY_DIRS} )
#	include_directories( ${HDF5_INCLUDE_DIR} )
#
#	target_link_libraries(src PUBLIC ${HDF5_CXX_LIBRARIES})
#endif()

target_include_directories(src PUBLIC ../include)
