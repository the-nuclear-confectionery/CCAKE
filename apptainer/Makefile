all: images/ubuntu.sif images/ccake.sif	 ../CCAKE/EoS/Houston/thermo.h5 images/ccake-dev.sif build-release build-debug

.PHONY: build-release build-debug run-gubser

# ==================================
# Recipes with example runs
# ==================================

run-gubser: images/ccake.sif ../CCAKE/build/ccake
	singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c " \
	cd /CCAKE/build ; \
	mkdir -p initial_conditions ; \
	../utilities/gubser_ic_converter.py ../misc/Gubser_checks/ac/Initial_Profile_tau=1fm.dat initial_conditions/gubser_tau_1.0.dat ; \
	./ccake ../input/input_parameters_gubser_checks.yaml output_gubser ; \
	../plotting_scripts/plot_gubser.py ../misc/Gubser_checks/ac output_gubser gubser.png ;"

run-trento: images/ccake.sif ../CCAKE/build/ccake
	singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c " \
	cd /CCAKE/build ; \
	mkdir -p initial_conditions ; \
	trento Pb Pb -o trento_example --grid-max 14 --grid-step=.02 --b-max 6 --b-min 5 --random-seed 42 --norm 50; \
	../utilities/trento2ccake.py trento_example/0.dat 0.1 14 initial_conditions/ccake-trento.dat; \
	rm -r trento_example; \
	./ccake ../input/input_parameters_ccake.yaml output_trento ; \
	"

# ==================================
# Recipes for building images
# ==================================
images/ubuntu.sif: defs/ubuntu.def 
	singularity build --fakeroot $@ $<

images/ccake.sif: defs/ccake.def images/ubuntu.sif
	singularity build --fakeroot $@ defs/ccake.def

images/ccake-dev.sif: defs/ccake-dev.def images/ccake.sif
	singularity build --fakeroot $@ defs/ccake-dev.def

# ==================================
# Recipes for compiling
# ==================================
../CCAKE/build/ccake: images/ccake.sif
	mkdir -p ../CCAKE/build; \
	singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c "cd /CCAKE/build ; cmake .. ; make -j";
	cd ../CCAKE/build; \
	ln -sf ../EoS ./

../CCAKE/build-debug/ccake: images/ccake-dev.sif
	mkdir -p ../CCAKE/build-debug
	singularity exec --bind ../CCAKE:/CCAKE images/ccake-dev.sif bash -c "cd /CCAKE/build-debug ; cmake -DCMAKE_BUILD_TYPE='Debug' .. ; make -j"
	cd ../CCAKE/build-debug; \
	ln -sf ../EoS ./

build-release: ../CCAKE/build/ccake
build-debug: ../CCAKE/build-debug/ccake

../CCAKE/EoS/Houston/thermo.h5:
	mkdir -p ../CCAKE/EoS/Houston
	wget -O $@  https://zenodo.org/record/6829115/files/thermo.h5?download=1

