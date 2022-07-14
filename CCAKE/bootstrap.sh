#!/bin/bash
#----------

ICdir=initial_conditions
EoSdir=EoS/Houston

#-------------------------------------------------------------------------------
# set up example ICCING initial conditions
# if IC file does not exist
if [ ! -f ${ICdir}/ICCING_example.dat ]
then
  mkdir -p ${ICdir}

  ICCINGexample="https://zenodo.org/record/6829000/files/ICCING_example.dat?download=1"
  wget -O ${ICdir}/ICCING_example.dat ${ICCINGexample}
fi

#-------------------------------------------------------------------------------
# set up lattice-based EoS tables
# if TXT file does not exist
if [ ! -f ${EoSdir}/thermo.dat ]
then
  mkdir -p ${EoSdir}
  EoSTXTTable="https://zenodo.org/record/6829115/files/thermo.dat?download=1"
  wget -O ${EoSdir}/thermo.dat ${EoSTXTTable}
fi
# if HDF file does not exist
if [ ! -f ${EoSdir}/thermo.h5 ]
then
  mkdir -p ${EoSdir}
  EoSHDFTable="https://zenodo.org/record/6829115/files/thermo.h5?download=1"
  wget -O ${EoSdir}/thermo.h5 ${EoSHDFTable}
fi

#-------------------------------------------------------------------------------
# make all shell scripts executable
chmod +x *.sh


#-------------------------------------------------------------------------------
# compile with CMake
./build_all.sh "$@"