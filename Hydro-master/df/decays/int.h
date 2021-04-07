#ifndef INT_H_ 
#define INT_H_

double	gauss(int  n, double f(double, void *),double  xlo, double xhi, void	*optvec );/* gauss(int n, double (*func)(), double xlo, xhi, void *optvec); */
void	gausspts( int n, double xlo,double  xhi, double *xvec, double *wvec ); /* gausspts( int n, double xlo, xhi; double *x, *w ); */
double	gaussn(int  n, int ndiv, double	(*f)(double, void *), double xlo, double xhi,void *optvec );/* gaussn(int n, ndiv, double (*func)(),
			double xlo, double xhi, double *para ); */
double	gala(int  n, double	(*f)(double, void *), double xlo, double invslope,void  *optvec );	/* gala(int n, double (*func)(), double xlo,invslope,void *optvec);*/
double	gahe(int  n, double	(*f)(double, void *), double center, double  width, void *optvec );	/* gahe(int n, double (*func)(), double center, width,void *optvec );*/
double	gauche( int n, double	(*f)(double, void *), double base, double pole, void * optvec );/* gauche(int n,double (*func)(),double base,pole,void *optvec ); */
double	gaussp( int n, double	(*f)(double, void *),  double xlo, double xhi, double	para[] );
double	galap(int  n, double	(*f)(double, void *), double xlo, double invslope, double	para[] );
double	gahep( int n, double	(*f)(double, void *), double center,double  width,double	*para);

double	gaussnbyn(int  n, double (*f)(double, void *), double xlo, double xhi, double ylo, double yhi );	/* not yet tested */
#endif
