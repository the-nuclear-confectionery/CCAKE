all: images/ubuntu.sif images/ccake.sif	 ../CCAKE/EoS/Houston/thermo.h5 images/ccake-dev.sif build-release build-debug

.PHONY: build-release build-debug run-gubser build-gpu run-trento

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

run-trento: images/ccake.sif ../CCAKE/build/ccake ../CCAKE/EoS/Houston/thermo.h5
	singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c " \
	cd /CCAKE/build ; \
	mkdir -p initial_conditions ; \
	trento Pb Pb -o trento_example --grid-max 14 --grid-step=.02 --b-max 6 --b-min 5 --random-seed 42 --norm 50; \
	../utilities/trento2ccake.py trento_example/0.dat 0.1 14 initial_conditions/ccake-trento.dat; \
	rm -r trento_example; \
	./ccake ../input/input_parameters_ccake.yaml output_trento ; \
	"

run: images/ccake-dev.sif ../CCAKE/build/ccake
	singularity exec --bind ../CCAKE:/CCAKE images/ccake-dev.sif zsh

run-debug: images/ccake-dev.sif ../CCAKE/build-debug/ccake
	singularity exec --bind ../CCAKE:/CCAKE images/ccake-dev.sif zsh

run-1D: images/ccake.sif ../CCAKE/build/ccake
	singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c " \
	cd /CCAKE/build ; \
	./ccake ../input/input_parameters_1D_check.yaml  output_1D "

run-gubser-gpu:
	singularity exec --nv --nvccli --bind ../CCAKE:/CCAKE images/ccake-gpu.sif bash -c " \
	cd /CCAKE/build ; \
	mkdir -p initial_conditions ; \
	../utilities/gubser_ic_converter.py ../misc/Gubser_checks/ac/Initial_Profile_tau=1fm.dat initial_conditions/gubser_tau_1.0.dat ; \
	./ccake ../input/input_parameters_gubser_checks.yaml output_gubser ; \
	../plotting_scripts/plot_gubser.py ../misc/Gubser_checks/ac output_gubser gubser.png ;"

run-gubser-charges-cosh: images/ccake.sif ../CCAKE/build/ccake
	singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c " \
	cd /CCAKE/build ; \
	./ccake ../input/input_parameters_gubser_checks_finite_B_cosh.yaml output_gubser_tau0_1p0_mu0_50_T0_200-cosh ;"

run-gubser-charges-conformal: images/ccake.sif ../CCAKE/build/ccake
	singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c " \
	cd /CCAKE/build ; \
	./ccake ../input/input_parameters_gubser_checks_finite_B_conformal.yaml output_gubser_tau0_1p0_mu0_50_T0_200-conformal ;"

run-gubser-debug: images/ccake-dev.sif ../CCAKE/build-debug/ccake
	singularity exec --bind ../CCAKE:/CCAKE images/ccake-dev.sif bash -c " \
	cd /CCAKE/build-debug ; \
	gdb --args ./ccake ../input/input_parameters_gubser_checks_finite_B_conformal.yaml output_gubser_tau0_1p0_mu0_50_T0_200-conformal ;"


# ==================================
# Recipes for building images
# ==================================
images/ubuntu.sif: defs/ubuntu.def
	singularity build --fakeroot $@ $<

images/ccake.sif: defs/ccake.def images/ubuntu.sif
	singularity build --fakeroot $@ defs/ccake.def

images/ccake-dev.sif: defs/ccake-dev.def images/ccake.sif
	singularity build --fakeroot $@ defs/ccake-dev.def

images/nvhpc.sif: defs/nvhpc.def
	singularity build --nv --fakeroot $@ $<

images/ccake-gpu.sif: defs/ccake-gpu.def images/nvhpc.sif
	singularity build --nv --fakeroot $@ defs/ccake-gpu.def

# ==================================
# Recipes for compiling
# ==================================
../CCAKE/build/ccake: images/ccake.sif ../CCAKE/src/* ../CCAKE/include/*
	mkdir -p ../CCAKE/build; \
	singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c "cd /CCAKE/build ; cmake .. ; make -j";
	cd ../CCAKE/build; \
	ln -sf ../EoS ./

../CCAKE/build-debug/ccake: images/ccake-dev.sif ../CCAKE/src/* ../CCAKE/include/*
	mkdir -p ../CCAKE/build-debug
	singularity exec --bind ../CCAKE:/CCAKE images/ccake-dev.sif bash -c "cd /CCAKE/build-debug ; cmake -DCMAKE_BUILD_TYPE='Debug' .. ; make -j"
	cd ../CCAKE/build-debug; \
	ln -sf ../EoS ./

build-release: ../CCAKE/build/ccake
build-debug: ../CCAKE/build-debug/ccake

../CCAKE/EoS/Houston/thermo.h5:
	mkdir -p ../CCAKE/EoS/Houston
	wget -O $@  https://zenodo.org/record/6829115/files/thermo.h5?download=1

build-gpu: images/ccake-gpu.sif
	mkdir -p ../CCAKE/build-gpu
	singularity exec --nv --bind ../CCAKE:/CCAKE images/ccake-gpu.sif bash -c "cd /CCAKE/build-gpu ; cmake .. ; make -j"
	cd ../CCAKE/build-gpu; \
	ln -sf ../EoS ./