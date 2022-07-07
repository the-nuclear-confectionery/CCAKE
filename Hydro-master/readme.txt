How to compile:

chmod +x cleancompile.sh
./cleancompile.sh 

How to run an event (one test event is included):

./charm.sh PbPb2_EOS21_0 0 0 shear EOS21/0 trento PbPb5020TeV

where the order is:  
name of the input file with parameters (i.e. the PbPb2_EOS21_0 goes into inputPbPb2_EOS21_0.dat where the format is ions, EOS, and bin)
first event
last event
viscosity type (shear includes both bulk and shear)
output folder/bin 
initial condition type
detector
If everything works you should have the spectra of stable particles within NHhydro/df/out/trento/PbPb5020TeV/shear/EOS21/0/ev0dsbvc_dNdphidpp.dat

within Nhydro/df you can also calculate the flow harmonics using:

./vns/fo PbPb2_EOS21_0.dat 0 0 chgv1 0.3 3 1

where the order is:
full name of the input folder (within df/input - input folders for decays are slightly different than for hydro)
first event
last event
particles to observe (all=all charge particles)
lower pt cut
upper pt cut
decays on or off (1=decays, 0=no decays)

The papers to cite for v-USPhydro:

Phys.Rev. C90 (2014) no.3, 034907 
Phys.Rev. C88 (2013) 044916 

The decays section is an adapted version of AZhydro:
P. F. Kolb, J. Sollfrank, and U. Heinz, Phys. Rev. C 62
(2000) 054909; P. F. Kolb and R. Rapp, Phys. Rev. C 67
(2003) 044903; P. F. Kolb and U. Heinz, nucl-th/0305084.

