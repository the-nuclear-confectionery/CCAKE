#ifndef DRIVER_H
#define DRIVER_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

//#ifdef __cpluscplus
//   extern "C" {
//#endif

__BEGIN_DECLS

// initialization functions
void initialize(const char parameters_filename[]);
void initialize_thermodynamics(const char parameters_filename[]);

// these functions accept input in MeV and output in physical units
// and in the (non-standard) ordering 
void get_eBSQ_densities(double point[], double densities[]);
void get_sBSQ_densities(double point[], double densities[]);
void get_full_thermo(double point[], double thermodynamics[]);

// these functions expect input/output in units of fm
// and in the (standard) ordering of T,muB,muQ,muS
void STANDARD_get_eBSQ_densities(double point[], double densities[]);
void STANDARD_get_sBSQ_densities(double point[], double densities[]);
void STANDARD_get_full_thermo(double point[], double thermodynamics[]);

__END_DECLS
//#ifdef __cplusplus
//   }
//#endif

#endif
