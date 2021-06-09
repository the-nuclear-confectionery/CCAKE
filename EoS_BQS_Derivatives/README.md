# EoS_BQS
QCD Equation of State with baryon number, electric charge and strangeness chemical potential

This code was created by Paolo Parotto and Jamie Stafford, Department of Physics, University of Houston, Houston, TX 77204, US.

This code produces an EoS with chemical potentials associated to baryon number, electric charge and strangeness. 

INPUT:
The required input are the parameters for the parametrization of the Taylor coefficients up to order mu^4 (22 coefficients), contained in "Coefficients_Paramters.dat".

THE OUTPUT: 
The code creates a folder "Coefficients_Checks", where the 22 coefficients are printed as functions of the temperature in the range T= 30 - 800 MeV, as well as their first and second derivatives wrt the temperature. This is done in order to check the correct behavior of the parametrized coefficients and their derivatives, and can be convenient for plotting.
Then, it creates a folder "Thermodynamics", where:
- It creates a file "EoS_Taylor_AllMu", where thermodynamic quantities are calculated on the whole 4D computational grid. 
	In order: T, muB, muQ, muS, Pressure, entropy density, baryon density, strangeness density, charge density, energy density, speed of sound;
	
- It creates a file "AllTherm_No_QS_Taylor.dat", where thermodynamic quantities are calculated on the 2D (T,muB) computational grid, with muQ=muS=0. 
	In order: T, muB, muQ, muS, Pressure, entropy density, baryon density, strangeness density, charge density, energy density, speed of sound;
	
- It creates a file "AllTherm_StrNeutr_Taylor", where thermodynamic quantities are calculated on the 2D (T,muB) computational grid, with strangeness neutrality. 
	In order: T, muB, muQ, muS, Pressure, entropy density, baryon density, strangeness density, charge density, energy density, speed of sound;
	
- It creates a file "AllTherm_StrNeutr_muBTConst_Taylor.dat", where thermodynamic quantities are calculated for T = 30 - 800 MeV, along lines of constant muB/T, for muB/T = 0.5,1,1.5,2,2.5,3, with strangeness neutrality.
	In order: muB/T, T, muB, muQ, muS, Pressure, entropy density, baryon density, strangeness density, charge density, energy density, speed of sound;
	
- It creates a file "AllTherm_NoQS_muBTConst_Taylor.dat", where thermodynamic quantities are calculated for T = 30 - 800 MeV, along lines of constant muB/T, for muB/T = 0.5,1,1.5,2,2.5,3, with muQ=muS=0.
	In order: muB/T, T, muB, muQ, muS, Pressure, entropy density, baryon density, strangeness density, charge density, energy density, speed of sound;


COMPILING This code uses the standard libraries <stdio.h>, <stdlib.h>, <math.h>, <string.h>, <time.h>. In addition, the libraries sys/types.h, sys/stat.h and unistd.h were included in order to use the functions mkdir() and chdir(). 
You can compile and then run in main directory with:

```bash
make
./EoS_BQS Coefficients_Parameters.dat
```


CONTACT For problems, debugging, contact Paolo Parotto at paolo.parotto@gmail.com

WHEN USING THIS CODE When using this code, the following work should be cited: http://inspirehep.net/record/1720588
