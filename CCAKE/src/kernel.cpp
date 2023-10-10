#include "kernel.h"

namespace ccake
{

template class SPHkernel<1>;
template class SPHkernel<2>;
template class SPHkernel<3>;
  
template<> const double SPHkernel<2>::knorm  = 10./(7.*M_PI);
template<> const double SPHkernel<2>::knorm2 = 5./(14.*M_PI);
template<> const double SPHkernel<2>::kgrad  = -15./(14.*M_PI);
template<> const double SPHkernel<2>::kgrad2 = 10./(7.*M_PI);

template<> const double SPHkernel<3>::knorm  = 1./(4.*M_PI); 
template<> const double SPHkernel<3>::knorm2 = 1./(16.*M_PI); 
template<> const double SPHkernel<3>::kgrad  = -3./(14.*M_PI); 
template<> const double SPHkernel<3>::kgrad2 = 1./(4.*M_PI); 

template<> const double SPHkernel<1>::knorm  = 1./6.;
template<> const double SPHkernel<1>::knorm2 = 1./24.;
template<> const double SPHkernel<1>::kgrad  = -3./24.; 
template<> const double SPHkernel<1>::kgrad2 = 1./6.; 



template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::kernel(double distance, double hT )
{ 
  double norm       = knorm/pow(hT,D);
  double norm2      = knorm2/pow(hT,D);
  double q = distance/hT;
  if ( q >= 2.0 )
    return 0.0;
  if ( q >= 1.0 )
    return norm2*(2.0-q)*(2.0-q)*(2.0-q);

  double qq=q*q;
  return norm*(1 - 1.5*qq + 0.75*q*qq);
}

template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::distance(const double* r1, const double* r2)
{
  double d=0; 
  for (int idir=0; idir<D; ++idir){
    double diff = r1[idir]-r2[idir];
    d += diff*diff;
  }
  return sqrt(d);
}

/// @brief Calculates the gradient of the kernel function
/// @tparam D The dimensionality of the simulation
/// @param rel_dist The relative distance between two particles
/// @param r Norm of the relative distance
/// @param h SPH smoothing length
/// @param grad The array to store the gradient of the kernel
template<unsigned int D>
KOKKOS_FUNCTION
void SPHkernel<D>::gradKernel(double const* rel_dist, double r, double h, double* grad )
{
  double norm_grad  = kgrad/pow(h,D+1);
  double norm_grad2 = kgrad2/pow(h,D+1);
  double q = r/h;

  if ( q >= 2.0 )
  {
    for (int idir=0; idir<2; idir++)
      grad[idir] = 0.0;
    return;
  }
  if ( q >= 1.0 )
  {
    for (int idir=0; idir<2; idir++)
      grad[idir] = (norm_grad/r)*(2.0-q)*(2.0-q)*rel_dist[idir];
    return;
  }
  for (int idir=0; idir<2; idir++)
      grad[idir] = norm_grad2*( -3.0+2.25*q )*rel_dist[idir];
  return;
}  


  




}