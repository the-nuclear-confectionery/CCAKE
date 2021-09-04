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

void initialize(const char parameters_filename[]);
void get_densities(double point[], double densities[]);

__END_DECLS
//#ifdef __cplusplus
//   }
//#endif

#endif
