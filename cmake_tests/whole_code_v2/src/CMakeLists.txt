add_library(src OBJECT
	BSQHydro.cpp
	eos.cpp
	eos_derivatives.cpp
	eos_extension.cpp
	eos_initialization.cpp
	eos_table.cpp
	equations_of_motion.cpp
	input_output.cpp
	kernel.cpp
	linklist.cpp
	matrix.cpp
	new_matrix.cpp
	new_vector.cpp
	particle.cpp
	persephone.cpp
	rootfinder.cpp
	runge_kutta.cpp
	sph_workstation.cpp
	system_state.cpp
	vector.cpp
)

find_package(GSL)
target_link_libraries(src PUBLIC ${GSL_LIBRARIES})
target_include_directories(src PUBLIC ${GSL_INCLUDE_DIR})

target_include_directories(src PUBLIC .)
