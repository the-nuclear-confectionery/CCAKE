#include "eom_default.h"

using namespace ccake;

/// @brief Calculates the gamma factor = u^0.
/// @details This is the default implementation of the gamma factor calculation.
/// It assumes that the last component of the velocity vector is the longitudinal
/// velocity, and that the metric is Milne.
/// @param u The velocity vector (space components only).
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return 
template<unsigned int D>
double EoM_default<D>::gamma_calc(double *u, const double &time_squared) {
    return sqrt(1.0+dot(u,u,time_squared));
}

/// @brief Computes the inner product of two vectors.
/// @details This considers only the space components of the vectors. This is the 
/// implementation for the special case D=2.
/// @param u The first vector
/// @param v The second vector
/// @param time_squared The square of the time where the gamma factor will be computed
/// @return u^i v^i
template<>
double EoM_default<2>::dot(double* u, double* v, const double &time_squared) {
  double s = 0;
  for (unsigned int i=0; i<2; i++)
    s+= u[i]*v[i];
  return s;
}

/// @brief Transforms a scalar from the lab frame to the LRF.
/// @tparam D 
/// @param lab The quantity in the lab (computational) frame.
/// @param gamma Lorentz contraction factor. 
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return The quantity in the fluid LRF (local rest frame).
template<unsigned int D>
double EoM_default<D>::get_LRF(const double &lab, const double &gamma, 
                               const double &time_squared) {
                                return lab/gamma/time_squared;
}

/// @brief Computes the inner product of two vectors.
/// @details This considers only the space components of the vectors. This is the 
/// general case. It assumes that the last component of the velocity vector is the
/// longitudinal velocity.
/// @param u The first vector
/// @param v The second vector
/// @param time_squared The square of the time where the gamma factor will be computed
/// @return u^i v^i
template<unsigned int D>
double EoM_default<D>::dot(double* u,double* v, const double &time_squared) {
  double s = 0;
  for (unsigned int i=0; i<D-1; i++)
    s+= u[i]*v[i];
  s += u[D-1]*v[D-1]*time_squared;
  return s;
}

/// @brief Ensures that the shear tensor will be traceless.
/// @details This is the default implementation of the traceless condition. The last 
/// component is assumed to be the longitudinal one and is modified as to ensure the tensor
/// is traceless.
/// @param pi_diag A D+1 dimensional array containing the diagonal elements of the shear tensor.
/// @return pi^{DD} = pi^{00} - \sum_{i=1}^{D-1} pi^{ii}/tau^2
template<unsigned int D>
double EoM_default<D>::get_shvDD(double* pi_diag, const double &time_squared){
    double s = pi_diag[0];
    for (unsigned int i=1; i<D; i++)
        s -= pi_diag[i];
    s /= time_squared;
    return s;
}

/// @brief Ensures that the shear tensor will be traceless.
/// @details This is the implementation for the special case (2+1)D. It assumes the existence of
/// a longitudinal component. We are using Milne coordinates.
/// @param pi_diag A D+1 dimensional array containing the diagonal elements of the shear tensor.
/// @return pi^{33} = pi^{00} - \sum_{i=1}^{D} pi^{ii}/\tau^2
/// \todo I may be wrong about this implementation. It is worth to double check.
template<>
double EoM_default<2>::get_shvDD(double* pi_diag, const double &time_squared){
    double s = pi_diag[0];
    for (unsigned int i=1; i<3; i++)
        s -= pi_diag[i];
    s /= time_squared;
    return s;
}

/// @brief Get the distance between two particles.
/// @details This is the default implementation of the distance calculation between two
/// particles. This is done in Milne coordinates.
/// @param x1 The position of the first particle.
/// @param x2 The position of the second particle.
/// @param time_squared The square of the time where the distance will be computed.
/// @return The distance between the two particles.
template<unsigned int D>
double EoM_default<D>::get_distance(double* x1, double* x2, const double &time_squared){
    double s = 0;
    for (unsigned int i=0; i<D-1; i++)
        s += (x1[i]-x2[i])*(x1[i]-x2[i]);
    s += (x1[D-1]-x2[D-1])*(x1[D-1]-x2[D-1])*time_squared;
    return sqrt(s);
}