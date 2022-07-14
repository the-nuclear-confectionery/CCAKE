#!/bin/bash
#----------

ICdir=initial_conditions
EoSdir=EoS/Houston

#-------------------------------------------------------------------------------
# set up example ICCING initial conditions
mkdir -p ${ICdir}

ICCINGexample="https://zenodo.org/record/6829000/files/ICCING_example.dat?download=1"
wget -O ${ICdir}/ICCING_example.dat ${ICCINGexample}

#-------------------------------------------------------------------------------
# set up lattice-based EoS table (text and HDF format)
mkdir -p ${EoSdir}

EoSTXTTable="https://zenodo.org/record/6829115/files/thermo.dat?download=1"
EoSHDFTable="https://zenodo.org/record/6829115/files/thermo.h5?download=1"
wget -O ${EoSdir}/thermo.dat ${EoSTXTTable}
wget -O ${EoSdir}/thermo.h5 ${EoSHDFTable}


#-------------------------------------------------------------------------------
# make all shell scripts executable
chmod +x *.sh


#-------------------------------------------------------------------------------
# compile with CMake
./build_all.sh