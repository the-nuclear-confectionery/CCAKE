The attached files contain semi-analytic solutions of Israel-Stewart theory in the Gubser flow regime. 

The file "Initial_Profile_tau=1fm.dat" contains the initial condition and is designed for a grid size of dx=dy=0.05fm. 

The files "y=0_tau=1.2_SemiAnalytic.dat", "y=0_tau=1.5_SemiAnalytic.dat", and "y=0_tau=2.0_SemiAnalytic.dat" contain the semi-analytic solution at times tau=1.2 fm, 1.5 fm, and 2.0 fm respectively, and y=0.

The files "y=x_tau=1.2_SemiAnalytic.dat", "y=x_tau=1.5_SemiAnalytic.dat", and "y=x_tau=2.0_SemiAnalytic.dat" contain the semi-analytic solution at times tau=1.2 fm, 1.5 fm, and 2.0 fm respectively, and y=x.
 
All files are organized as: 

"x (fm)", "y (fm)", "T (GeV)", "u^x", "u^y", "pi^xx (GeV/fm^3)", "pi^yy (GeV/fm^3)", "pi^xy (GeV/fm^3)", "\tau^2 pi^\eta\eta (GeV/fm^3)"

where T is the Temperature, u^x is the x component of the 4-velocity, and pi^xx is the xx component of the shear stress tensor

The initial time for the initial profile in "Initial_Profile_tau=1fm.dat" is tau_0 = 1 fm. The shear viscosity to entropy density ratio is eta/s=0.2, while the shear relaxation time is tau_R = 5*eta/(e+P). Here, e is the energy density and P is the thermodynamic pressure. 
The equation of state must be the one of a conformal fluid, i.e., e=3P. These transport coefficients and equation of state have to be used in order to describe the semi-analytic solutions contained here and displayed in the paper. The actual equations of motion that must be solved correspond to equations (2) and (3) of the paper. 
