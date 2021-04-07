How to use the BSQ Equation of State code for Hydro:


To use the class, create an eos object. You can use the eos(string, string, int) constructor or run eos() followed by init(string string int).


NOTE: The code takes in an Equation of State in dimensionless quantities and IMMEDIATELY converts it to units of MeV and fm. All of the values returned by the functions have full units. 
quantityFile must be in dimensionless quantities and must be formatted "T  muB  muQ  muS  p  s  B  S  Q  e  cs2" 
derivFile should be formatted "T  muB  muQ  muS  d2p/dB2  d2p/dQ2  d2p/dS2  d2p/dBdQ  d2p/dBdS d2p/dQdS  d2p/dTdB  d2p/dTdQ  d2p/dTdS  d2p/dT2"
The EoS_BSQ code from Claudia's group (with the derivative addition from Debora) is formatted correctly. Other Equations of State may not be

Use the tbqs function to change where the splines are evaluated. For example if you want to find the interpolated values at T = 500 MeV, muB = 50 MeV, muQ = 83.7 MeV and muS = 9 MeV, you would call tbqs(500, 50, 83.7, 9);. Each time tbqs is called, it evaluates and initializes all of the variables at that point during the function call. Once you have called tbqs on a point, you don't need to call it again until you switch points, or use any of the functions which change the evaluation point which are the entropy based functions: update_s, s_out, and s_terms_T.

NOTE: Be careful with BSQ vs BQS. Because the equation of state input files order mu as (muB, muQ, muS), and rho as (rhoB, rhoS, rhoQ), this is how they are ordered throughout the code. Function calls with mu will always be BQS, function calls with rho will always be BSQ.


The functions p(), s(), B(), S(), Q(), e(), cs2(), T(), muB(), muQ(), muS(), dwds(), dwdB(), dwdS() and dwdQ() are getter functions which return their respective variable at the point determined by the last tbqs() call.


cs2out(setT, setmuB, setmuS, setmuQ) is identical to calling tbqs at the given point followed by cs2().


wfz returns the enthalpy at the given (T,muB,muQ,muS).


update_s and s_out are functions that use the rootfinder:
update_s finds a desired point in T and mu from a given entropy and BSQ, and calls tbqs to evaluate the splines at that point.
s_out returns the entropy at the location which matches the given energy and BSQ. The splines are initialized at that point in T and mu, and the entropy at that point is returned.

Rootfinder note: These functions are not guaranteed to succeed. If the rootfinder does not find the given point, tbqs is not called and the T and mu location will not change, so be careful! If the rootfinder fails, the point is not in the EoS range, it is too far from the initial guess, or there is too much error in the interpolation. The latter issue can be solved by raising the TOLERANCE parameter in the eos.h file or by using an input EoS with a denser grid. TOLERANCE = 1e-11 was the value used when testing the ICCING initial conditions.

The parameters STEPS and TYPE may also be updated depending on the Equation of State. STEPS is the maximum number of iterations the rootfinder will take before failing automatically. The rootfinder historically will either converge or fail before reaching STEPS = 1000, but if a sufficiently large or complex equation of state was used it may need to be increased. TYPE indicates the rootfinding method used. Options for type are provided by the GSL multiroots library. Thy type options are dnewton, hybrids, hybrid, and broyden. The hybrids method was chosen because it does not crash the program when it fails to converge as the dnewton and broyden methods do, and it is faster than the hybrid method.


s_terms_T(setT) returns the entropy at T = setT and muB = muQ = muS = 0


The eosin, w, A, efreeze and sfreeze either do nothing or return 0. Do not use them until they have been updated.


NOTE for hydro-code in 4D. The following function signatures have changed to include BSQ, so they will cause compilation errors until they are changed on the hydro end:

update_s(double sets); -> update_s(double sets, double setB, double setS, double setQ);
cs2out(double setT); -> cs2out(double setT, double setMuB, double setMuQ, double setMuS);
wfz(double setT); -> wfz(double setT, double setMuB, double setMuQ, double setMuS);
s_out(double setE); -> s_out(double setE, double setB, double setS, double setQ);
