all: images/ubuntu.sif images/ccake.sif  ../EoS/Houston/thermo.h5 images/ccake-dev.sif build-release build-debug
  
.PHONY: build-release build-debug run-gubser build-gpu run-trento run-bjorken

# ==================================
# Recipes with example runs
# ==================================

run-long: images/ccake.sif ../build/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c " \
        cd /CCAKE/build ; \
        ./ccake ../config/config_long.yaml output_long ;"
        
run-iccing: images/ccake.sif ../build/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c " \
        cd /CCAKE/build ; \
	./ccake ../input/input_parameters_iccing.yaml output_iccing ;"

run-bjorken: images/ccake.sif ../build/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c " \
        cd /CCAKE/build ; \
	./ccake ../input/config_bjorken.yaml output_bjorken ;"

run-landau: images/ccake.sif ../build/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c " \
        cd /CCAKE/build ; \
	./ccake ../input/config_landau.yaml output_landau ;"

run-bjorken-db: images/ccake-dev.sif ../build-debug/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake-dev.sif bash -c " \
        cd /CCAKE/build-debug ; \
        ./ccake ../config/config_bjorken.yaml output_bjorken ;"

run-long-debug: images/ccake-dev.sif ../build-debug/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake-dev.sif bash -c " \
        cd /CCAKE/build-debug ; \
        ./ccake ../config/config_long.yaml output_long ;"

run-gubser: images/ccake.sif ../build/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c " \
        cd /CCAKE/build ; \
        ./ccake ../input/input_parameters_gubser_checks.yaml output_gubser ;"

run-gubser-db: images/ccake.sif ../build-debug/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake-dev.sif bash -c " \
        cd /CCAKE/build-debug ; \
        mkdir -p initial_conditions ; \
        ../utilities/gubser_ic_converter.py ../misc/Gubser_checks/ac/Initial_Profile_tau=1fm.dat initial_conditions/gubser_tau_1.0.dat ; \
        ./ccake ../input/input_parameters_gubser_checks.yaml output_gubser ; "

run-trento: images/ccake.sif ../build/ccake ../EoS/Houston/thermo.h5
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c " \
        cd /CCAKE/build ; \
        mkdir -p initial_conditions ; \
        trento Pb Pb -o trento_example --grid-max 14 --grid-step=.02 --b-max 6 --b-min 5 --random-seed 42 --norm 50; \
        ../utilities/trento2ccake.py trento_example/0.dat 0.1 14 initial_conditions/ccake-trento.dat; \
        rm -r trento_example; \
        ./ccake ../input/input_parameters_ccake.yaml output_trento ; \
        "

run: images/ccake-dev.sif ../build/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake-dev.sif zsh

run-debug: images/ccake-dev.sif ../build-debug/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake-dev.sif zsh

run-1D: images/ccake.sif ../build/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c " \
        cd /CCAKE/build ; \
        ./ccake ../input/input_parameters_1D_check.yaml  output_1D "

run-gubser-gpu:
	singularity exec --cleanenv --nv --nvccli --bind ../:/CCAKE images/ccake-gpu.sif bash -c " \
	cd /CCAKE/build ; \
	mkdir -p initial_conditions ; \
	../utilities/gubser_ic_converter.py ../misc/Gubser_checks/ac/Initial_Profile_tau=1fm.dat initial_conditions/gubser_tau_1.0.dat ; \
	./ccake ../input/input_parameters_gubser_checks.yaml output_gubser ; \
        ../plotting_scripts/plot_gubser.py ../misc/Gubser_checks/ac output_gubser gubser.png ;"

run-gubser-charges-cosh: images/ccake.sif ../build/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c " \
        cd /CCAKE/build ; \
	cd /CCAKE/build ; \
        ./ccake ../input/input_parameters_gubser_checks_finite_B_cosh.yaml output_gubser_tau0_1p0_mu0_50_T0_200-cosh ;"

run-gubser-charges-conformal: images/ccake.sif ../build/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c " \
        cd /CCAKE/build ; \
        ./ccake ../input/input_parameters_gubser_checks_finite_B_conformal.yaml output_gubser_tau0_1p0_mu0_50_T0_250-conformal ;"

dbg-gubser-charges-conformal: images/ccake-dev.sif ../build-debug/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake-dev.sif bash -c " \
        cd /CCAKE/build-debug ; \
        gdb --args ./ccake ../input/input_parameters_gubser_checks_finite_B_conformal.yaml output_gubser_tau0_1p0_mu0_50_T0_250-conformal ;"


run-gubser-debug: images/ccake-dev.sif ../build-debug/ccake
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake-dev.sif bash -c " \
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
../build/ccake: images/ccake.sif ../src/* ../include/*
	mkdir -p ../build; \
        singularity exec --cleanenv --bind ../:/CCAKE images/ccake.sif bash -c "cd /CCAKE/build ; cmake .. ; make -j";
	cd ../build; \
        ln -sf ../EoS ./

../build-debug/ccake: images/ccake-dev.sif ../src/* ../include/*
	mkdir -p ../build-debug
	singularity exec --cleanenv --bind ../:/CCAKE images/ccake-dev.sif bash -c "cd /CCAKE/build-debug ; cmake -DCMAKE_BUILD_TYPE='Debug' .. ; make -j"
	cd ../build-debug; \
        ln -sf ../EoS ./

build-release: ../build/ccake ../EoS/Houston/thermo.h5
build-debug: ../build-debug/ccake

../EoS/Houston/thermo.h5:
	mkdir -p ../EoS/Houston
	wget -O $@  https://zenodo.org/record/6829115/files/thermo.h5?download=1

build-gpu: images/ccake-gpu.sif
	mkdir -p ../build-gpu
	singularity exec --cleanenv --nv --bind ../:/CCAKE images/ccake-gpu.sif bash -c "cd /CCAKE/build-gpu ; cmake .. ; make -j"
	cd ../build-gpu; \
        ln -sf ../EoS ./


