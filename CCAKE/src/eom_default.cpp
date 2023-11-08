#include "eom_default.h"

namespace ccake{
  template class EoM_default<1>;
  template class EoM_default<2>;
  template class EoM_default<3>;

template<unsigned int D>
void EoM_default<D>::reset_pi_tensor(std::shared_ptr<SystemState<D>> sysPtr){
  double t2 = (sysPtr->t)*(sysPtr->t);
  CREATE_VIEW(device_,sysPtr->cabana_particles)

  #ifdef DEBUG
   auto particles_host =
        Cabana::create_mirror_view_and_copy( Kokkos::HostSpace(), sysPtr->cabana_particles);
  auto h_matrix_view =  Cabana::slice<particle_info::hydro_space_matrix_info>(particles_host);
  auto h_vec_view =  Cabana::slice<particle_info::hydro_vector_info>(particles_host);
  auto h_scalar_view =  Cabana::slice<particle_info::hydro_scalar_info>(particles_host);
  auto h_spacetime_matrix_view =  Cabana::slice<particle_info::hydro_spacetime_matrix_info>(particles_host);
  Kokkos::fence();
  
  //Print hydro scalars for a particle in the midfle of the event
  int c=11000;
  //Scalars
  cout << "---- Before regularization of shear ----" << std::endl;
  for(int i=0; i<D+1; ++i){
    for(int j=0; j<D+1; ++j) std::cout << h_spacetime_matrix_view(c,hydro_info::shv, i, j) << " ";
    cout << endl;
  }
  #endif
  ///TODO: Needs checking for 1+1 and 3+1 cases
  auto kokkos_ensure_consistency = KOKKOS_LAMBDA(const int is, int ia) 
  {
    //Declare caches
    milne::Vector<double,D> u;
    milne::Matrix<double,D+1,D+1> shv;

    for(int idir=0; idir<D; ++idir) u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
    for(int idir=0; idir<D+1; ++idir)
    for(int jdir=0; jdir<D+1; ++jdir)
      shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);

    //computations
    double gamma = Kokkos::sqrt(1+milne::inner(u,u));
    //pi^{0i} = \pi^{ij}u_j/gamma
    milne::Vector<double,D> u_cov = u;
    u_cov.make_covariant(t2);
    for(int idir=1; idir<D+1; ++idir){
      milne::Vector<double,D> colp1_shv;
      for(int jdir=1; jdir<D+1; ++jdir) colp1_shv(jdir-1) = shv(idir,jdir);
      shv(0,idir) = 1./gamma*milne::inner(u,colp1_shv);
    } 
    
    //Simmetrizes
    for( int i=0; i<D+1; i++ )
    for( int j=i+1; j<D+1; j++ )
      shv(j,i) = shv(i,j);
    
    milne::Matrix<double,D,D> pimin, uu;
    for( int idir=0; idir<D; idir++ )
    for( int jdir=0; jdir<D; jdir++ ){
      pimin(idir,jdir) = shv(idir+1,jdir+1);
      uu(idir,jdir) = u(idir)*u(jdir);
    }

    //Fill data into cabana data structure before proceeding
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir)
    {
      device_hydro_space_matrix.access(is, ia, hydro_info::uu, idir, jdir) = uu(idir,jdir);
      device_hydro_space_matrix.access(is, ia, hydro_info::pimin, idir, jdir) = pimin(idir,jdir);
    }

    uu.make_covariant(0,t2);
    uu.make_covariant(1,t2);
    
    //pi^00 = u_i u_j pi^{ij}/gamma^2
    shv(0,0) = 1./gamma/gamma*milne::con(uu,pimin);
  
    //pi^33 = (pi^00 - pi^11 - pi^22)/t2
    double shv33 = shv(0,0);
    switch (D)
    {
    case 2:
      shv33 = (shv(0,0) - shv(1,1) - shv(2,2))/t2;
      break;
    case 3:
      shv33 = (shv(0,0) - shv(1,1) - shv(2,2))/t2;
      shv(3,3) = shv33;
      break;
    default:
      shv33 = shv(0,0)/t2;
      shv(1,1) = shv33;
      break;
    }
    

    //Return shv
    for(int idir=0; idir<D+1; ++idir)
    for(int jdir=0; jdir<D+1; ++jdir)
      device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir) = shv(idir,jdir) ;
    
  };
    //Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace> simd_policy(0, particles.size());
  Cabana::simd_parallel_for(*(sysPtr->simd_policy),kokkos_ensure_consistency,"kokkos_ensure_consistency");
  Kokkos::fence();
  #ifdef DEBUG
  Cabana::deep_copy( particles_host, sysPtr->cabana_particles);
  Kokkos::fence();
  
  //Print hydro scalars for a particle in the midfle of the event
  //Scalars
  cout << "---- After regularization of shear ----" << std::endl;
  cout << "tt\ttx\tty\txt\txx\txy\tyt\tyx\tyy" << endl;
  for(int i=0; i<D+1; ++i){
    for(int j=0; j<D+1; ++j) std::cout << h_spacetime_matrix_view(c,hydro_info::shv, j, i) << " ";
    cout << endl;
  }
  #endif
}

/// @brief Calculates the gamma factor = u^0.
/// @details This is the default implementation of the gamma factor calculation.
/// It assumes that the last component of the velocity vector is the longitudinal
/// velocity, and that the metric is Milne.
/// @param u The velocity vector (space components only).
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return the value of gamma
template<unsigned int D> KOKKOS_FUNCTION
double EoM_default<D>::gamma_calc(double u[D], const double &time_squared) {
    return sqrt(1.0+EoM_default<D>::dot(u,u,time_squared));
}

/// @brief Computes the inner product of two vectors.
/// @details This considers only the space components of the vectors. This is the
/// implementation for the special case D=2. IMPORTANT: We are not carrying the
/// metric signal here, that is, we are evaluating -v_j u^j.
/// @param u The first vector
/// @param v The second vector
/// @param time_squared The square of the time where the gamma factor will be computed
/// @return u^i v^i
template<> KOKKOS_FUNCTION
double EoM_default<2>::dot(double v[2], double u[2], const double &time_squared) {
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
template<unsigned int D> KOKKOS_FUNCTION
double EoM_default<D>::get_LRF(const double &lab, const double &gamma,
                               const double &t) {
                                return lab/gamma/t;
}

/// @brief Computes the inner product of two vectors.
/// @details This considers only the space components of the vectors. This is the
/// general case. It assumes that the last component of the velocity vector is the
/// longitudinal velocity. IMPORTANT: We are not carrying the metric signal here, that
/// is, we are evaluating -v_j u^j.
/// @param u The first vector
/// @param v The second vector
/// @param time_squared The square of the time where the gamma factor will be computed
/// @return u^i v^i
template<unsigned int D>
KOKKOS_FUNCTION
double EoM_default<D>::dot(double v[D], double u[D], const double &time_squared) {
  double s = 0;
  for (unsigned int i=0; i<D-1; i++)
    s+= u[i]*v[i];
  s += u[D-1]*v[D-1]*time_squared;
  return s;
}

/// @brief Contracts a vector with a vector (LHS contraction).
/// @details We are assuming a contraction in the form x^j = v_i T^{i j}. Bear in
/// mind that we are considering only the space components.
/// @param v The vector
/// @param T The tensor
/// @param time_squared The square of the time where the gamma factor will be computed
/// @param x Vector where the result will be stored
template<unsigned int D>
KOKKOS_FUNCTION
void EoM_default<D>::dot(double v[D],double T[D][D], const double &time_squared, double *x) {
  for (unsigned int j=0; j<D; j++){
    x[j] = 0;
    for (unsigned int i=0; i<D-1; i++){
      x[j]-= v[i]*T[i][j];
    }
    x[j] -= v[D-1]*T[D-1][j]*time_squared;
  }
}

/// @brief Contracts a vector with a vector (LHS contraction).
/// @details We are assuming a contraction in the form x^j = v_i T^{i j}. Bear in
/// mind that we are considering only the space components.
/// @param v The vector
/// @param T The tensor
/// @param time_squared The square of the time where the gamma factor will be computed
/// @param x Vector where the result will be stored
template<>
KOKKOS_FUNCTION
void EoM_default<2>::dot(double v[2],double T[2][2], const double &time_squared, double *x) {
  for (unsigned int j=0; j<2; j++){
    x[j] = 0;
    for (unsigned int i=0; i<2; i++){
      x[j]-= v[i]*T[i][j];
    }
  }
}

/// @brief Contracts all indices of 2 rank 2 contravariant tensors
/// @details Computes $ a = M^{ij}g_{ij}T^{ij}
/// @param M The first vector
/// @param T The second vector
/// @param time_squared The square of the time where the gamma factor will be computed
/// @return A double with the result of the contraction
template<>
KOKKOS_FUNCTION
double EoM_default<2>::full_contraction(double M[2][2], double T[2][2], const double &time_squared) {
  double r = 0;
  for (unsigned int i=0; i<2; i++)
  for (unsigned int j=0; j<2; j++)
    r = -M[i][j]*T[i][j];
  return r;
}

/// @brief Contracts all indices of 2 rank 2 contravariant tensors
/// @details Computes $ a = M^{ij}g_{ij}T^{ij}
/// @param M The first vector
/// @param T The second vector
/// @param time_squared The square of the time where the gamma factor will be computed
/// @return A double with the result of the contraction
template<unsigned int D>
KOKKOS_FUNCTION
double EoM_default<D>::full_contraction(double M[D][D], double T[D][D], const double &time_squared) {
  double r = 0;
  for (unsigned int i=0; i<D-1; i++){
    for (unsigned int j=0; j<D-1; j++)
      r = - M[i][j]*T[i][j];
    r -= M[D-1][i]*T[D-1][i];
    r -= M[i][D-1]*T[i][D-1];
  }
  r -= M[D-1][D-1]*T[D-1][D-1];
  return r;
}

/// @brief Ensures that the shear tensor will be traceless.
/// @details This is the default implementation of the traceless condition. The last
/// component is assumed to be the longitudinal one and is modified as to ensure the tensor
/// is traceless.
/// @param pi_diag A D+1 dimensional array containing the diagonal elements of the shear tensor.
/// @return pi^{DD} = pi^{00} - \sum_{i=1}^{D-1} pi^{ii}/tau^2
template<unsigned int D>
KOKKOS_FUNCTION
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
KOKKOS_FUNCTION
double EoM_default<2>::get_shvDD(double* pi_diag, const double &time_squared){
    double s = pi_diag[0];
    for (unsigned int i=1; i<3; i++)
        s -= pi_diag[i];
    s /= time_squared;
    return s;
}

template<unsigned int D>
void EoM_default<D>::evaluate_time_derivatives( std::shared_ptr<SystemState<D>> sysPtr)
{
  double t = (sysPtr->t);
  double t2 = (sysPtr->t);
  CREATE_VIEW(device_,sysPtr->cabana_particles);
  bool using_shear = true; ///TODO: This should be retrieved from the input parameters
  
  #ifdef DEBUG
  
  auto d_M = Cabana::AoSoA< Cabana::MemberTypes<double[D][D]>, DeviceType, VECTOR_LENGTH>("d_M",sysPtr->cabana_particles.size()) ;
  auto dM = Cabana::slice<0>(d_M);
  auto h_M = Cabana::AoSoA< Cabana::MemberTypes<double[D][D]>, HostType, VECTOR_LENGTH>("d_M",sysPtr->cabana_particles.size()) ;
  auto hM = Cabana::slice<0>(h_M);
  #endif
  
  auto fill_auxiliary_variables = KOKKOS_LAMBDA(int const is, int const ia){

    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double dwds = device_thermo.access(is, ia, thermo_info::dwds);
    double dwdB = device_thermo.access(is, ia, thermo_info::dwdB);
    double dwdQ = device_thermo.access(is, ia, thermo_info::dwdQ);
    double dwdS = device_thermo.access(is, ia, thermo_info::dwdS);
    double T = device_thermo.access(is, ia, thermo_info::T);
    double w = device_thermo.access(is, ia, thermo_info::w);
    double s = device_thermo.access(is, ia, thermo_info::s);
    double rhoB = device_thermo.access(is, ia, thermo_info::rhoB);
    double rhoQ = device_thermo.access(is, ia, thermo_info::rhoQ);
    double rhoS = device_thermo.access(is, ia, thermo_info::rhoS);
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    double Bulk = device_hydro_scalar.access(is, ia, hydro_info::Bulk);
    double setas = device_hydro_scalar.access(is, ia, hydro_info::setas);
    double stauRelax = device_hydro_scalar.access(is, ia, hydro_info::stauRelax);
    double tauRelax = device_hydro_scalar.access(is, ia, hydro_info::tauRelax);
    double zeta = device_hydro_scalar.access(is, ia, hydro_info::zeta);
    double shv00 = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, 0,0);

    milne::Matrix<double,D,D> gradV;
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir)
      gradV(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, jdir);

    
    double dsigma_dt = -sigma * milne::tr(gradV,t2); //Ok
    double g2        = gamma*gamma; //Ok
    double gt        = gamma*t; //Ok
    double dwdsT1    = 1 - dwds/T; //Ok
    double bigPI     = Bulk*sigma/gt; //Ok
    double C         = w + bigPI;

    double eta_o_tau = using_shear ? setas/stauRelax : 0.0; //Ok

    double Agam      = w - dwds*( s + bigPI/T ) - zeta/tauRelax
                     - dwdB*rhoB - dwdS*rhoS - dwdQ*rhoQ;

    device_hydro_scalar.access(is, ia, hydro_info::dsigma_dt)     = dsigma_dt;
    device_hydro_scalar.access(is, ia, hydro_info::gamma_squared) = g2;
    device_hydro_scalar.access(is, ia, hydro_info::gamma_cube)    = gamma*g2;
    device_hydro_scalar.access(is, ia, hydro_info::gamma_tau)     = gt;
    device_hydro_scalar.access(is, ia, hydro_info::dwdsT1)        = dwdsT1;
    device_hydro_scalar.access(is, ia, hydro_info::sigl)          = dsigma_dt/sigma - 1/t;
    device_hydro_scalar.access(is, ia, hydro_info::bigPI)         = bigPI;
    device_hydro_scalar.access(is, ia, hydro_info::C)             = C;
    device_hydro_scalar.access(is, ia, hydro_info::eta_o_tau)     = eta_o_tau;
    device_hydro_scalar.access(is, ia, hydro_info::Agam)          = Agam;
    device_hydro_scalar.access(is, ia, hydro_info::Agam2)         = ( Agam - eta_o_tau/3.0 - dwdsT1*shv00 ) / gamma;
    device_hydro_scalar.access(is, ia, hydro_info::Ctot)          = C + eta_o_tau*(1.0/g2-1.0);
  };

  Cabana::simd_parallel_for(*(sysPtr->simd_policy), fill_auxiliary_variables, "fill_auxiliary_variables");
  Kokkos::fence();

  auto fill_Btot = KOKKOS_LAMBDA(const int is, const int ia){
    auto bsub = 0;
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double Agam = device_hydro_scalar.access(is, ia, hydro_info::Agam);
    double Agam2 = device_hydro_scalar.access(is, ia, hydro_info::Agam2);
    double eta_o_tau = device_hydro_scalar.access(is, ia, hydro_info::eta_o_tau);
    double sigl = device_hydro_scalar.access(is, ia, hydro_info::sigl);
    double bigPI = device_hydro_scalar.access(is, ia, hydro_info::bigPI);
    double tauRelax = device_hydro_scalar.access(is, ia, hydro_info::tauRelax);
    double dwdsT1 = device_hydro_scalar.access(is, ia, hydro_info::dwdsT1);
    double gt = device_hydro_scalar.access(is, ia, hydro_info::gamma_tau);
    double shv33 = device_hydro_scalar.access(is, ia, hydro_info::shv33);
    if (using_shear){
        
      milne::Vector<double,D> v, u; //Already filled
      milne::Matrix<double,D+1,D+1> shv; //Already filled 
      milne::Matrix<double,D,D> gradV, uu, pimin; //Already filled
      milne::Matrix<double,D,D> piu, piutot, gradU; //Filled only if shear is enabled
        
      for(int idir=0; idir<D; ++idir){
        v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
        u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
        for(int jdir=0; jdir<D; ++jdir){
          pimin(idir, jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::pimin, idir, jdir );
          uu(idir, jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::uu, idir, jdir );
          gradV(idir, jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, jdir );
        }
      }
      for(int idir=0; idir<D+1; ++idir)
      for(int jdir=0; jdir<D+1; ++jdir)
        shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir );
      double g2 = device_hydro_scalar.access(is, ia, hydro_info::gamma_squared);
      double g3 = device_hydro_scalar.access(is, ia, hydro_info::gamma_cube);


      piu         = milne::rowp1(0,shv)*u;
      piutot      = piu + milne::transpose(piu);
      double pig  = shv(0,0)/g2;

      gradU        = gamma*gradV+g3*(v*(v*gradV));
      for (int i=0; i<=1; i++)
      for (int j=0; j<=1; j++)
        bsub += gradU(i,j) * ( pimin(i,j) + pig*uu(j,i)
                          - ( piu(i,j) + piu(j,i) ) /gamma );
      //Fill cache arrays
      for(int idir=0; idir<D; ++idir)
      for(int jdir=0; jdir<D; ++jdir){
        device_hydro_space_matrix.access(is, ia, hydro_info::gradU, idir, jdir) = gradU(idir, jdir);
        device_hydro_space_matrix.access(is, ia, hydro_info::piu, idir, jdir) = piu(idir, jdir);
        device_hydro_space_matrix.access(is, ia, hydro_info::piutot, idir, jdir) = piutot(idir, jdir);
      }
    }
    double dwdsT = 1-dwdsT1;
    device_hydro_scalar.access(is, ia, hydro_info::Btot) = 
      ( Agam*gamma + 2.0*eta_o_tau/3.0*gamma )*sigl + bigPI/tauRelax + dwdsT*( gt*shv33 + bsub );
  };
  
  Cabana::simd_parallel_for(*(sysPtr->simd_policy), fill_Btot, "fill_Btot");
  Kokkos::fence();

  auto compute_velocity_derivative = KOKKOS_LAMBDA(const int is, const int ia){
          // THIS IS THE ORIGINAL PART IN TIME DERIVATIVES
      
      double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);

      double Agam2 = device_hydro_scalar.access(is, ia, hydro_info::Agam2);
      double Ctot = device_hydro_scalar.access(is, ia, hydro_info::Ctot);
      double g2 = device_hydro_scalar.access(is, ia, hydro_info::gamma_squared);
      double dwdsT1 = device_hydro_scalar.access(is, ia, hydro_info::dwdsT1);
      double Btot = device_hydro_scalar.access(is, ia, hydro_info::Btot);
      milne::Vector<double,D> u, gradshear, gradP, gradBulk, divshear;
      milne::Matrix<double,D,D> gradU,uu,piu, pimin;
      milne::Matrix<double,D+1,D+1> shv;
      for(int idir=0; idir<D; ++idir){
        for(int jdir=0; jdir<D; ++jdir){
          gradU(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradU, idir, jdir);
          uu(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::uu, idir, jdir);
          piu(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::piu, idir, jdir);
          pimin(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::pimin, idir, jdir);
          shv(idir+1,jdir+1) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir+1, jdir+1);
        }
        u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
        gradP(idir) = device_hydro_vector.access(is, ia, hydro_info::gradP, idir);
        gradBulk(idir) = device_hydro_vector.access(is, ia, hydro_info::gradBulk, idir);
        divshear(idir) = device_hydro_vector.access(is, ia, hydro_info::divshear, idir);
        gradshear(idir) = device_hydro_vector.access(is, ia, hydro_info::gradshear, idir);
        shv(0,idir+1) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, 0, idir+1);
        shv(idir+1,0) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir+1, 0);
      }
      shv(0,0) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, 0, 0);
      

      // set the Mass and the Force
      milne::Matrix <double,D,D> M = Agam2*uu - (1+4/3./g2)*piu
              + dwdsT1*transpose(piu) + gamma*pimin;
      for(int idir=0; idir<D; ++idir) M(idir, idir) += Ctot*gamma;

      milne::Vector<double,D> F    = Btot*u + gradshear
                              - ( gradP + gradBulk + divshear );

      if ( using_shear )
      {
        double stauRelax = device_hydro_scalar.access(is, ia, hydro_info::stauRelax);
        double eta_o_tau = device_hydro_scalar.access(is, ia, hydro_info::eta_o_tau);
        double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
        double dsigma_dt = device_hydro_scalar.access(is, ia, hydro_info::dsigma_dt);
        milne::Vector<double, D> v;
        for(int idir=0; idir<D; ++idir)
          v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
        double gamt = 1.0/gamma/stauRelax;
        double pre  = eta_o_tau/gamma;
        double p1   = gamt - 4.0/3.0/sigma*dsigma_dt + 1.0/t/3.0;
        milne::Matrix <double,D,D> partU = gradU + milne::transpose( gradU );
        milne::Vector<double,D> minshv   = milne::rowp1(0, shv);
        F += pre*v*partU + p1*minshv;
      }
      
      milne::Matrix<double,D,D> MI = milne::inverse(M);
      milne::Vector<double,D> du_dt = F*MI;
      for(int idir=0; idir<D; ++idir) device_hydro_vector.access(is,ia,hydro_info::du_dt, idir) = du_dt(idir);
  };
  Cabana::simd_parallel_for(*(sysPtr->simd_policy), compute_velocity_derivative, "compute_velocity_derivative");
  Kokkos::fence();
  #ifdef DEBUG
  auto particles_host =
        Cabana::create_mirror_view_and_copy( Kokkos::HostSpace(), sysPtr->cabana_particles);
  auto h_matrix_view =  Cabana::slice<particle_info::hydro_space_matrix_info>(particles_host);
  auto h_vec_view =  Cabana::slice<particle_info::hydro_vector_info>(particles_host);
  auto h_scalar_view =  Cabana::slice<particle_info::hydro_scalar_info>(particles_host);
  auto h_spacetime_matrix_view =  Cabana::slice<particle_info::hydro_spacetime_matrix_info>(particles_host);
  auto h_ID =  Cabana::slice<particle_info::id>(particles_host);
  auto h_pos =  Cabana::slice<particle_info::position>(particles_host);
  Cabana::deep_copy(h_M,d_M);

  //Print hydro scalars for a particle in the midfle of the event
  int c=11000;
  //Scalars
  cout << "---- Hydro Info Scalars ----" << std::endl;
  cout << "particle_id..: " << h_ID(c) << std::endl;
  cout << "position.....: " << h_pos(c,0) << " " << h_pos(c,1) << std::endl;
  cout << "dsigma_dt....: " << h_scalar_view(c,hydro_info::dsigma_dt) << std::endl;
  cout << "gamma........: " << h_scalar_view(c,hydro_info::gamma)  << std::endl;
  cout << "gamma squared: " << h_scalar_view(c,hydro_info::gamma_squared)  << std::endl;
  cout << "gamma cube...: " << h_scalar_view(c,hydro_info::gamma_cube)  << std::endl;
  cout << "gamma tau....: " << h_scalar_view(c,hydro_info::gamma_tau)  << std::endl;
  cout << "dwdsT1.......: " << h_scalar_view(c,hydro_info::dwdsT1) << std::endl;
  cout << "sigl.........: " << h_scalar_view(c,hydro_info::sigl) << std::endl;
  cout << "bigPI........: " << h_scalar_view(c,hydro_info::bigPI) << std::endl;
  cout << "eta_o_tau....: " << h_scalar_view(c,hydro_info::eta_o_tau) << std::endl;
  cout << "Agam.........: " << h_scalar_view(c,hydro_info::Agam)  << std::endl;
  cout << "Agam2........: " << h_scalar_view(c,hydro_info::Agam2) << std::endl;
  cout << "C............: " << h_scalar_view(c,hydro_info::C) << std::endl;
  cout << "Ctot.........: " << h_scalar_view(c,hydro_info::Ctot) << std::endl;
  cout << "Btot.........: " << h_scalar_view(c,hydro_info::Btot) << std::endl;
  cout << "u............: ";
  for(int j=0;  j<D; ++j)
    cout << h_vec_view(c,hydro_info::u, j) << " ";
  cout << endl;
  cout << "du_dt........: ";
  for(int j=0;  j<D; ++j)
    cout << h_vec_view(c,hydro_info::du_dt, j) << " ";
  cout << endl;
  cout << "shv........: " << endl;
  for(int j=0;  j<D+1; ++j){
    for(int i=0;  i<D+1; ++i)
      cout << h_spacetime_matrix_view(c,hydro_info::shv, i, j) << " ";
    cout << endl;
  }
  cout << "gradU........: " << endl;
  for(int j=0;  j<D; ++j){
    for(int i=0;  i<D; ++i)
      cout << h_matrix_view(c,hydro_info::gradU, i, j) << " ";
    cout << endl;
  }
  cout << "gradV........: " << endl;
  for(int j=0;  j<D; ++j){
    for(int i=0;  i<D; ++i)
      cout << h_matrix_view(c,hydro_info::gradV, i, j) << " ";
    cout << endl;
  }
  cout << "piu........: " << endl;
  for(int j=0;  j<D; ++j){
    for(int i=0;  i<D; ++i)
      cout << h_matrix_view(c,hydro_info::piu, i, j) << " ";
    cout << endl;
  }
  cout << "piutot.....: " << endl;
  for(int j=0;  j<D; ++j){
    for(int i=0;  i<D; ++i)
      cout << h_matrix_view(c,hydro_info::piu, i, j) << " ";
    cout << endl;
  }

  #endif
}
/*
    double Bsub    = 0;
    if ( using_shear )
    {
      milne::Matrix<double,D+1,D+1> shv;
      for (int i=0; i<=1; i++)
      for (int j=0; j<=1; j++){
        
      }

      milne::Matrix<double,D,D> gradU = gamma*gradV+g3*(v*(v*gradV));
      // these quantities will all be zero without shear
      milne::Matrix<double,D,D> piu  = milne::rowp1(0,shv)*u;
      milne::Matrix<double,D,D> piutot = piu+milne::transpose(piu);
      double pig  = shv(0,0)/g2;

      for (int i=0; i<=1; i++)
      for (int j=0; j<=1; j++)
        Bsub += gradU(i,j) * ( pimin(i,j) + pig*uu(j,i)
             - ( piu(i,j) + piu(j,i) ) / gamma );
        //Initialize shear dependent variables, if necessary

      }

      double Btot    = ( Agam*gamma + 2.0*eta_o_tau/3.0*gamma )*sigl
                          + bigPI/tauRelax
                          + dwdsT*( gt*shv33 + Bsub );
  }

  auto update_velocity = KOKKOS_LAMBDA(int const i, int const is){
    
    
    //Load caches
    double gamma, g2;
    gamma = device_hydro_scalar(is, ia, hydro_info::gamma);
    g2 = gamma*gamma;

    milne::Matrix<double,D,D> pimin, uu;
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir){
      pimin(idir, jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::pimin, idir, jdir);
      uu(idir, jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::uu, idir, jdir);
    }
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir){
      pimin(idir, jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::pimin, idir, jdir);
      uu(idir, jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::uu, idir, jdir);
    }


    milne::Matrix<double,D,D> piu = rowp1(0,shv)*u;
    milne::Matrix <double,D,D> M = Agam2*uu - (1+4/3./g2)*piu
              + dwdsT1*milne::transpose(piu) + gamma*pimin;
    for(int idir=0, idir<D; ++idir) M(idir,idir) += Ctot*gamma;

    milne::Matrix<double,D,D> MInverse = milne::inverse(M);
    milne::Vector<double,D> du_dt;
    du_dt = F*MInverse;

    //TODO:
    //- Compute:
    // F
    // MInverse - Ok
    // M - Ok
    // Agam2
    // uu - Ok
    // g2 - Ok
    // piu
    // dwdsT1
    // gamma - Ok
    // pimin - Ok
    // shv
    // u


    
    for(int idir=0; idir<D; ++idir)
      device_hydro_vector.access(is, ia, hydro_info::du_dt, idir) = du_dt(idir);
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir){
      device_hydro_space_matrix.access(is, ia, hydro_info::piu, idir, jdir) = piu(idir, jdir);
    }
    
  }
  #ifdef DEBUG
  //Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace> simd_policy(0, particles.size());
  //Cabana::simd_parallel_for(simd_policy,single_particle_update,"single_particle_update");
  int c = 110000;
  auto particles_host =
        Cabana::create_mirror_view_and_copy( Kokkos::HostSpace(), particles);
  auto h_matrix_view =  Cabana::slice<particle_info::hydro_space_matrix_info>(particles_host);
  auto h_vec_view =  Cabana::slice<particle_info::hydro_vector_info>(particles_host);
  auto h_scalar_view =  Cabana::slice<particle_info::hydro_scalar_info>(particles_host);
  auto h_spacetime_matrix_view =  Cabana::slice<particle_info::hydro_spacetime_matrix_info>(particles_host);
  Cabana::deep_copy(h_M,d_M);

  cout << "du_dt: " << h_vec_view(c,device_hydro::du_dt,0) <<" " <<h_vec_view(c,device_hydro::du_dt,1) << endl;
  
  #endif
  

  auto single_particle_update = KOKKOS_LAMBDA(int const iparticle)
  {

      //1) Cache locally the particle data
      double sigma = device_hydro_scalar(iparticle, ccake::hydro_info::sigma);
      double gamma = device_hydro_scalar(iparticle, ccake::hydro_info::gamma);
      double t_squared = t*t;
      double Temperature = device_thermo(iparticle, ccake::thermo_info::T);
      double dwds = device_thermo(iparticle, ccake::thermo_info::dwds);
      double Bulk = device_hydro_scalar(iparticle, ccake::hydro_info::Bulk);
      double w = device_thermo(iparticle, ccake::thermo_info::w);
      double setas = device_hydro_scalar(iparticle, ccake::hydro_info::setas);
      double stauRelax = device_hydro_scalar(iparticle, ccake::hydro_info::stauRelax);
      double s = device_thermo(iparticle, ccake::thermo_info::s);
      double zeta = device_hydro_scalar(iparticle, ccake::hydro_info::zeta);
      double tauRelax = device_hydro_scalar(iparticle, ccake::hydro_info::tauRelax);
      double dwdB = device_thermo(iparticle, ccake::thermo_info::dwdB);
      double rhoB = device_thermo(iparticle, ccake::thermo_info::rhoB);
      double dwdS = device_thermo(iparticle, ccake::thermo_info::dwdS);
      double rhoS = device_thermo(iparticle, ccake::thermo_info::rhoS);
      double dwdQ = device_thermo(iparticle, ccake::thermo_info::dwdQ);
      double rhoQ = device_thermo(iparticle, ccake::thermo_info::rhoQ);
      double shv33 = device_hydro_scalar(iparticle, ccake::hydro_info::shv33);

      milne::Vector<double,D> v, u, gradshear, gradP, gradBulk, divshear, du_dt;
      milne::Matrix<double, D, D> gradV,uu,pimin;
      milne::Matrix<double, D+1, D+1> shv;

      for (unsigned int idir=0; idir<D; idir++){
        v(idir) = device_hydro_vector(iparticle, ccake::hydro_info::v, idir);
        u(idir) = device_hydro_vector(iparticle, ccake::hydro_info::u, idir);
        //minshv(idir) = device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, 0, idir+1);
        gradshear(idir) = device_hydro_vector(iparticle, ccake::hydro_info::gradshear, idir);
        gradP(idir) = device_hydro_vector(iparticle, ccake::hydro_info::gradP, idir);
        gradBulk(idir) = device_hydro_vector(iparticle, ccake::hydro_info::gradBulk, idir);
        divshear(idir) = device_hydro_vector(iparticle, ccake::hydro_info::divshear, idir);
        du_dt(idir) = 0.;
        for (unsigned int jdir=0; jdir<D; jdir++){
          gradV(idir,jdir) = device_hydro_space_matrix(iparticle, ccake::hydro_info::gradV, idir, jdir);
          //Load quiantities computed in reset pi tensor 
          uu(idir,jdir)= device_hydro_space_matrix(iparticle, ccake::hydro_info::uu, idir, jdir);
          pimin(idir,jdir)= device_hydro_space_matrix(iparticle, ccake::hydro_info::pimin, idir, jdir);
        }
      }
      for (unsigned int idir=0; idir<D+1; idir++)
      for (unsigned int jdir=0; jdir<D+1; jdir++)
          shv(idir,jdir) = device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, idir, jdir);



      //2) Perform computations - Minimal changes when compared to CCAKE-v1 //For now, 2D only
      double dsigma_dt = -sigma * ( gradV(0,0) + gradV(1,1) );
      //double dsigma_dt = 0;
      //for(int idir=0;idir<D;++idir) dsigma_dt += gradV[idir][idir];
      //dsigma_dt *= -sigma; //3D version?
      double gamma_squared = gamma*gamma;
      double gamma_cube = gamma*gamma_squared;
      double gamma_tau = gamma*t;
      double dwdsT    = dwds/Temperature;
      double dwdsT1 = 1.0 - dwdsT;
      double sigl = dsigma_dt/sigma - 1/t;
      milne::Matrix<double, D, D> gradU = gamma*gradV+gamma_cube*(v*(v*gradV));
      double bigPi = Bulk*sigma/gamma_tau;
      double C = w + bigPi;

      double eta_o_tau    = using_shear ? setas/stauRelax : 0.0;

      double Agam = w - dwds*(s+bigPi/Temperature) - zeta/tauRelax
                    - dwdB*rhoB - dwdS*rhoS - dwdQ*rhoQ;
      double Agam2 = (Agam - eta_o_tau/3.0 - dwdsT1*shv(0,0) ) / gamma;
      double Ctot  = C + eta_o_tau*(1./gamma_squared-1.);

      //Previously in bsub function
      double bsub = 0;
      milne::Vector<double,D> minshv   = milne::rowp1(0, shv);
      milne::Matrix<double,D,D> piu, piutot;
      if (using_shear){
        piu    = minshv*u;
        piutot = piu+milne::transpose(piu);
        double pig  = shv(0,0)/gamma_squared;
        for (int i=0; i<D; i++)
        for (int j=0; j<D; j++)
          bsub += gradU(i,j) * ( pimin(i,j) + pig*uu(j,i)
                - piutot(i,j)  / gamma );
      }

      double Btot = ( Agam*gamma + 2.0*eta_o_tau/3.0*gamma )*sigl
                          + bigPi/tauRelax
                          + dwdsT*( gamma_tau*shv33 + bsub );

      double gamt = 0.0, pre = 0.0, p1 = 0.0;
      if ( using_shear )
      {
        gamt = 1.0/gamma/stauRelax;
        pre  = eta_o_tau/gamma;
        p1   = gamt - 4.0/3.0/sigma*dsigma_dt + 1.0/t/3.0;
      }

      milne::Matrix <double,D,D> partU = gradU + milne::transpose( gradU );

      milne::Matrix<double, D, D> M = Agam2*uu - (1+4/3./gamma_squared)*piu
              + dwdsT1*transpose(piu) + gamma*pimin;
      
      for (int idir=0; idir<D; idir++) M(idir,idir) += Ctot*gamma;

      milne::Vector<double,D> F  = Btot*u - gradP 
                                  //+ gradshear - ( gradBulk + divshear )
                                  ;
      if ( using_shear )
        F += pre*v*partU + p1*minshv;

      
      #ifdef DEBUG
      for (int idir=0; idir<D; idir++) 
      for (int jdir=0; jdir<D; jdir++) dM(iparticle,idir,jdir) = M(idir,jdir);
      #endif
      milne::Matrix<double, D, D> MInverse = milne::inverse(M);
      du_dt = F*MInverse;
      /*dot(v, gradV, t_squared, aux_vector);
      for (int idir=0; idir<D; idir++)
      for (int jdir=0; jdir<D; jdir++){
        M[idir][jdir] = 0;
        Ipi[idir][jdir] = 0;
        uu[idir][jdir] = u[idir]*u[jdir];
        gradU[idir][jdir] = gamma*gradV[jdir][idir] - gamma_cube*aux_vector[idir]*v[jdir];
        pi_u[idir][jdir] = shv[0][idir+1]*u[jdir];
        pimin[idir][jdir] = shv[idir+1][jdir+1];
      }
      for (int idir=0; idir<D; idir++)
      for (int jdir=0; jdir<D; jdir++){
        piutot[idir][jdir] = pi_u[idir][jdir] + pi_u[jdir][idir];
        partU[idir][jdir] = gradU[idir][jdir] + gradU[jdir][idir];
      }

      double bsub = 0;
      for (int idir=0; idir<D; idir++)
      for (int jdir=0; jdir<D; jdir++)
        bsub -= gradU[idir][jdir]*( pimin[idir][jdir]
                                    + uu[idir][jdir]*shv[0][0]/gamma_squared
                                    - piutot[idir][jdir]/gamma );


      double Btot = ( Agam*gamma + 2.0*eta_o_tau/3.0*gamma )*sigl ///TODO: 2/3 or 1/3. See Jaki's (Eq. 274)?
                      + bigPi/tauRelax
                      + dwdsT*( gamma_tau*shv33 + bsub );

      // set the Mass and the Force
      for (int idir=0; idir<D; idir++){
        M[idir][idir] = Ctot;
        for (int jdir=0; jdir<D; jdir++){
          M[idir][jdir] += Agam2*uu[idir][jdir]
                          -(1+4./3./gamma_squared)*pi_u[idir][jdir]
                          + dwdsT1*pi_u[jdir][idir] + gamma*pimin[idir][jdir]; 
        }
        F[idir] = Btot*u[idir] + gradshear[idir] - ( gradP[idir] + gradBulk[idir] + divshear[idir] );
      }

      //===============
      // shear contribution
      double gamt = 0.0, pre = 0.0, p1 = 0.0;
      ///TODO: Wrap on if for shear
      //if ( this->settingsPtr->using_shear )
        gamt = 1.0/gamma/stauRelax;
        pre  = eta_o_tau/gamma;
        p1   = gamt - 4.0/3.0/sigma*dsigma_dt + 1.0/t/3.0;
      //}
        dot(v, partU, t_squared, aux_vector);
        for (int idir=0; idir<D; idir++)
          F[idir] += p1*minshv[idir]-pre*aux_vector[idir];


      Utilities::inverse<D>(&M,&M_inverse);

      //===============
      // compute acceleration
      for (int idir=0; idir<D; idir++){
        for (int jdir=0; jdir<D; jdir++){
          du_dt[idir] += F[jdir]*M_inverse[idir][jdir];
        }
      }

      //===============
      // "coordinate" divergence
      double div_u = (1./ gamma)*dot(u, du_dt, t_squared)
                      - ( gamma/ sigma ) * dsigma_dt;

      //===============
      // "covariant" divergence
      double bigtheta = div_u*t + gamma;

      //===============
      for (int idir=0; idir<D; idir++){
        aux_vector[idir] = shv[0][0]*v[idir] - minshv[idir];
        for (int jdir=0; jdir<D; jdir++)
          sub[idir][jdir] = pimin[idir][jdir] + (shv[0][0]/gamma_squared)*uu[idir][jdir]
                            - piutot[idir][jdir]/gamma;
      }
      //TODO: Fill inside only if shear is on
      double aux = 0.;

      double inside = t*( dot( aux_vector, du_dt, t_squared ) - gamma*t*shv33+full_contraction(sub,gradU,t_squared) );
      double d_dt_specific_s = 1./sigma/Temperature*( -bigPi*bigtheta + inside );

      //formulating simple setup for Beta_Bulk derivative  //WMS: I don't know what is this
      //hi.finite_diff_cs2 = (ti.cs2 - hi.prev_cs2)/0.05; // Asadek
	    //hi.finite_diff_T   = (Temperature - hi.prev_T)/0.05;     // Asadek
	    //hi.finite_diff_w   = (ti.w - hi.prev_w)/0.05;     // Asadek
	    //hi.dBeta_dt        = 0.5*((-hi.finite_diff_T/(Temperature*Temperature))*(1/ti.w)*(1/((1/3-ti.cs2)*(1/3-ti.cs2))))
	    //                   + 0.5*((-hi.finite_diff_w/(ti.w*ti.w))*(1/Temperature)*(1/((1/3-ti.cs2)*(1/3-ti.cs2))))
			//	          	     + 0.5*((4*ti.cs2*hi.finite_diff_cs2*(1/((1/3-ti.cs2)*(1/3-ti.cs2)*(1/3-ti.cs2))))*(1/Temperature)*(1/ti.w));//Asadek

      double dBulk_dt = -(zeta*bigtheta/sigma + Bulk/gamma)/tauRelax;

      //===============
      // define auxiliary variables
      ///TODO: Write in terms of contractions
      double vduk = dot(v, du_dt, t_squared);
      double ulpi[D][D], ududt[D][D], vsub[D][D];
      for (int idir=0;idir<D;++idir){
        Ipi[idir][idir] = -2.*eta_o_tau/3.;
        for (int jdir=0;jdir<D;++jdir){

          ulpi[idir][jdir] = u[idir]*shv[0][jdir];
          Ipi[idir][jdir] += - 2.*eta_o_tau*uu[idir][jdir]/3.
                            + 4.*pimin[idir][jdir]/3.;
          ududt[idir][jdir] = u[idir]*du_dt[jdir];
          vsub[idir][jdir] = 0;
          for (int kdir=0;kdir<D;++kdir)
            vsub[idir][jdir] += (u[idir]*pimin[jdir][kdir]
                             + u[jdir]*pimin[idir][kdir])*du_dt[kdir];
        }
      }

      //===============
      ///TODO: - ADD READABLE TERM NAMES
      ///TODO: Wrap on if for shear
      //if ( this->settingsPtr->using_shear )
      for (int idir=0;idir<D;++idir)
      for (int jdir=0;jdir<D;++jdir)
        dshv_dt[idir][jdir] = - gamt*( pimin[idir][jdir] + setas*partU[idir][jdir] )
                              - eta_o_tau*( ududt[idir][jdir] + ududt[jdir][idir] )
                              + vsub[idir][jdir] + sigl*Ipi[idir][jdir]
                              - vduk*( ulpi[idir][jdir] + ulpi[jdir][idir] + Ipi[idir][jdir]/gamma );

    //3) Update the particle data
    ///TODO: Do I really need to update all of these?
    /
    device_hydro_scalar(iparticle, ccake::hydro_info::dsigma_dt) = dsigma_dt;
    device_hydro_scalar(iparticle, ccake::hydro_info::gamma_squared) = gamma_squared;
    device_hydro_scalar(iparticle, ccake::hydro_info::gamma_cube) = gamma_cube;
    device_hydro_scalar(iparticle, ccake::hydro_info::gamma_tau) = gamma_tau;
    device_hydro_scalar(iparticle, ccake::hydro_info::dwdsT1) = dwdsT1;
    device_hydro_scalar(iparticle, ccake::hydro_info::sigl) = sigl;
    device_hydro_scalar(iparticle, ccake::hydro_info::bigPI) = bigPi;
    device_hydro_scalar(iparticle, ccake::hydro_info::C) = C;
    device_hydro_scalar(iparticle, ccake::hydro_info::eta_o_tau) = eta_o_tau;
    device_hydro_scalar(iparticle, ccake::hydro_info::Agam) = Agam;
    device_hydro_scalar(iparticle, ccake::hydro_info::Agam2) = Agam2;
    device_hydro_scalar(iparticle, ccake::hydro_info::Ctot) = Ctot;
    device_hydro_scalar(iparticle, ccake::hydro_info::Btot) = Btot;
    // device_hydro_scalar(iparticle, ccake::hydro_info::div_u) = div_u;
    // device_hydro_scalar(iparticle, ccake::hydro_info::bigtheta) = bigtheta;
    // device_hydro_scalar(iparticle, ccake::hydro_info::inside) = inside;
    // device_hydro_scalar(iparticle, ccake::hydro_info::dBulk_dt) = dBulk_dt;
    // device_d_dt_spec(iparticle, ccake::densities_info::s) = d_dt_specific_s;


    
    for (int idir=0; idir<D; idir++){
      device_hydro_vector(iparticle, ccake::hydro_info::du_dt, idir) = du_dt(idir);
      for (int jdir=0; jdir<D; jdir++){
        device_hydro_space_matrix(iparticle, ccake::hydro_info::gradU, idir, jdir) = gradU(idir,jdir);
        device_hydro_space_matrix(iparticle, ccake::hydro_info::dshv_dt, idir, jdir) = 0;//dshv_dt[idir][jdir];
        device_hydro_space_matrix(iparticle, ccake::hydro_info::piu, idir, jdir) = piu(idir,jdir);
        device_hydro_space_matrix(iparticle, ccake::hydro_info::piutot, idir, jdir) = piutot(idir,jdir);
      }
    }
    
  };

  formatted_output::report("Start computing the derivatives");
  Kokkos::parallel_for( "single_particle_update", particles.size(),
                        single_particle_update );
  Kokkos::fence();
  formatted_output::report("Finish computing the derivatives");
  #ifdef DEBUG_VERBOSE
  //Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace> simd_policy(0, particles.size());
  //Cabana::simd_parallel_for(simd_policy,single_particle_update,"single_particle_update");
  auto particles_host =
        Cabana::create_mirror_view_and_copy( Kokkos::HostSpace(), particles);
  auto h_matrix_view =  Cabana::slice<particle_info::hydro_space_matrix_info>(particles_host);
  auto h_vec_view =  Cabana::slice<particle_info::hydro_vector_info>(particles_host);
  auto h_scalar_view =  Cabana::slice<particle_info::hydro_scalar_info>(particles_host);
  auto h_spacetime_matrix_view =  Cabana::slice<particle_info::hydro_spacetime_matrix_info>(particles_host);
  Cabana::deep_copy(h_M,d_M);

  //Print hydro scalars for a particle in the midfle of the event
  int c=19000;
  //Scalars
  cout << "---- Hydro Info Scalars ----" << std::endl;
  cout << "gamma........: " << h_scalar_view(c,hydro_info::gamma)  << std::endl;
  cout << "Agam.........: " << h_scalar_view(c,hydro_info::Agam)  << std::endl;
  cout << "Agam2........: " << h_scalar_view(c,hydro_info::Agam2) << std::endl;
  cout << "Btot.........: " << h_scalar_view(c,hydro_info::Btot) << std::endl;
  cout << "Ctot.........: " << h_scalar_view(c,hydro_info::Ctot) << std::endl;
  cout << "C............: " << h_scalar_view(c,hydro_info::C) << std::endl;
  cout << "dsigma_dt....: " << h_scalar_view(c,hydro_info::dsigma_dt) << std::endl;
  cout << "gamma_squared: " << h_scalar_view(c,hydro_info::gamma_squared) << std::endl;
  cout << "gamma_cube...: " << h_scalar_view(c,hydro_info::gamma_cube) << std::endl;
  cout << "sigl.........: " << h_scalar_view(c,hydro_info::sigl) << std::endl;
  cout << "bigPI........: " << h_scalar_view(c,hydro_info::bigPI) << std::endl;
  cout << "eta_o_tau....: " << h_scalar_view(c,hydro_info::eta_o_tau) << std::endl;
  cout << "dwdsT1.......: " << h_scalar_view(c,hydro_info::dwdsT1) << std::endl;
  cout << "gamma_tau....: " << h_scalar_view(c,hydro_info::gamma_tau) << std::endl;
  cout << "sigma........: " << h_scalar_view(c,hydro_info::sigma) << std::endl;
  cout << "shv33........: " << h_scalar_view(c,hydro_info::shv33) << std::endl;
  cout << "---- Hydro Info Vectors ----" << std::endl;
  cout << "u............: ";
  for(int i=0; i<D; ++i) std::cout << h_vec_view(c,hydro_info::u, i) << " ";
  std::cout << std::endl;
  cout << "du_dt........: ";
  for(int i=0; i<D; ++i) std::cout << h_vec_view(c,hydro_info::du_dt, i) << " ";
  std::cout << std::endl;

  cout << "---- Hydro Info Matrices ----" << std::endl;
  cout << "gradV........: ";
  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j) std::cout << h_matrix_view(c,hydro_info::gradV, i, j) << " ";
  std::cout << std::endl;
  cout << "gradU........: ";
  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j) std::cout << h_matrix_view(c,hydro_info::gradU, i, j) << " ";
  std::cout << std::endl;
  cout << "piu........: ";
  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j) std::cout << h_matrix_view(c,hydro_info::piu, i, j) << " ";
  std::cout << std::endl;
  cout << "pimin........: ";
  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j) std::cout << h_matrix_view(c,hydro_info::pimin, i, j) << " ";
  std::cout << std::endl;
  cout << "uu..........: ";
  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j) std::cout << h_matrix_view(c,hydro_info::uu, i, j) << " ";
  std::cout << std::endl;
  cout << "piu.........: ";
  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j) std::cout << h_matrix_view(c,hydro_info::piu, i, j) << " ";
  std::cout << std::endl;
    cout << "piutot....: ";
  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j) std::cout << h_matrix_view(c,hydro_info::piutot, i, j) << " ";
  std::cout << std::endl;
  cout << "---- Shear Matrix ----" << std::endl;
  for(int i=0; i<D+1; ++i)
  for(int j=0; j<D+1; ++j) std::cout << h_spacetime_matrix_view(c,hydro_info::shv, i, j) << " ";
  std::cout << std::endl;
  

  std::cout << "+++++++++++++Mass Matrix++++++++++" << endl;
  for(int i=0; i<D; ++i){
    for(int j=0; j<D; ++j)
      std::cout <<  hM(c,i,j) << " ";
  }
  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++" << endl;
  #endif
}*/
}