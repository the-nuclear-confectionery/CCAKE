#!/usr/bin/env bash
singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c "
cd /CCAKE/build
mkdir initial_conditions
../utilities/gubser_ic_converter.py ../misc/Gubser_checks/ac/Initial_Profile_tau=1fm.dat initial_conditions/gubser_tau_1.0.dat
./ccake ../input/input_parameters_gubser_checks.yaml output_gubser
../plotting_scripts/plot_gubser.py ../misc/Gubser_checks/ac output_gubser gubser.png"
