The organization of v-USPhydro is the following:

1. Generate Initial conditions (typically TRENTO)
2. Run the hydrodynamics+decays (outputs the spectra of all stable particles)
3. Run the flow analysis (outputs integrated vn's, vn(pT), factorization breaking, and sorts vn(pT) into centralities).  Note, this will need to be re-ran multiple times if you are interested in PID (paricle identification) so normally run chgv1 (all charged stable particles) and then the PID: pi+, p+, K+, lambda, cascade, omega
4. The integrated vn's can be analyzed offline to do cumulants, symmetric cumulants, event plane correlations etc
