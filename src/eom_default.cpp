#include "eom_default.h"
// #include <fstream>
// #include <iostream>

/// @file eom_default.cpp
/// @brief Implementation of the default equations of motion for the 
/// hydrodynamic evolution.

using namespace ccake;

template class EoM_default<1>;
template class EoM_default<2>;
template class EoM_default<3>;

/// @brief Update gamma and v
/// @tparam D The number of spatial dimensions.
template<unsigned int D>
void EoM_default<D>::update_velocity(std::shared_ptr<SystemState<D>> sysPtr)
{
  CREATE_VIEW(device_, sysPtr->cabana_particles)
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH, ExecutionSpace>(0, sysPtr->cabana_particles.size());
  double t =  (sysPtr->t);
  double t2 = t*t;
  auto update_gammas = KOKKOS_LAMBDA(const int is, int ia) 
  {
    milne::Vector<double,D> u;
    for(int idir=0; idir<D; ++idir) u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
    milne::Vector<double,D> u_cov = u;
    //Updates gamma and velocities, and sigma 
    u_cov.make_covariant(t2); //Transforms u^i to -u_i
    double gamma = Kokkos::sqrt(1-milne::contract(u_cov,u)); //Calculates gamma = \sqrt{1+u^i u_i}
    device_hydro_scalar.access(is,ia, hydro_info::gamma) = gamma;
    for(int idir=0; idir<D; ++idir) device_hydro_vector.access(is, ia, hydro_info::v, idir) = u(idir)/gamma;
  };

  Cabana::simd_parallel_for(simd_policy, update_gammas, "update_gamma_sigma_kernel");
  Kokkos::fence();
}


/// @brief Enforces the constraints for the shear viscous tensor \f$ \pi^{\mu\nu} \f$.
/// @details Enforces the constraints  \f$\pi^{\mu\nu}u_\nu = 0 \f$,
/// \f$\pi^\mu_\mu = 0 \f$ and \f$ \pi^{\mu\nu} = \pi^{\nu\mu} \f$.
/// We assume that the components \f$ \pi^{xx} \f$, \f$ \pi^{xy} \f$,
/// \f$ \pi^{x\eta} \f$, \f$ \pi^{yy} \f$ and \f$ \pi^{y\eta} \f$ are computed
/// during evolution (for 2+1D simulations, \f$ \pi^{x\eta} \f$ and 
/// \f$ \pi^{y\eta} \f$ are assumed to be zero).
/// \f[\begin{align*}
/// \pi^{0i} & = -\pi^{ij}u_j/\gamma \\
/// \pi^{00} & = u_i u_j \pi^{ij}/\gamma^2 \\
/// \end{align*}\f]
///
/// We also take the opportunity to update the gamma factor and the velocity
/// vector, as well the matrices \f$u^i u^j\f$ and cache variable
/// `pimin` = \f$\pi^{ij}\f$.
/// @tparam D The number of spatial dimensions.
/// @param sysPtr A pointer to an object of class SystemState.
template<unsigned int D>
void EoM_default<D>::reset_pi_tensor(std::shared_ptr<SystemState<D>> sysPtr)
{
  double t =  (sysPtr->t);
  double t2 = t*t;
  CREATE_VIEW(device_,sysPtr->cabana_particles)
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());

  auto kokkos_ensure_consistency = KOKKOS_LAMBDA(const int is, int ia) 
  { 
    milne::Vector<double,D> u;
    for(int idir=0; idir<D; ++idir) u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
    milne::Vector<double,D> u_cov = u;
    //Updates gamma and velocities, and sigma 
    u_cov.make_covariant(t2); //Transforms u^i to -u_i
    double gamma = device_hydro_scalar.access(is,ia, hydro_info::gamma);
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    //compute bulk from big bulk 
    double big_bulk = device_hydro_scalar.access(is, ia, hydro_info::bigBulk);
    device_hydro_scalar.access(is, ia, hydro_info::bulk) = big_bulk*sigma;
    //Declare caches
    milne::Matrix<double,4,4> shv;


    //compute shear from big shear
    for(int idir=0; idir<2; ++idir)
    for(int jdir=idir; jdir<3; ++jdir)
      shv(idir+1,jdir+1) = device_hydro_shear_aux_vector.access(is, ia, hydro_info::bigshv, idir, jdir)*sigma;



    //Symmetrizes space part
    for( int i=1; i<4; i++ )
    for( int j=i+1; j<4; j++ )
      shv(j,i) = shv(i,j);

    //pi^{0i} = -\pi^{ij}u_j/gamma 
    for(int idir=1; idir<4; ++idir){
      milne::Vector<double,3> shv_i;
      for(int jdir=1; jdir<4; ++jdir) shv_i(jdir-1) = shv(idir,jdir);
      shv(0,idir) = -1.*milne::contract(u_cov,shv_i)/gamma; 
    }

    //Symmetrizes time part
    for( int i=1; i<4; i++ )
      shv(i,0) = shv(0,i);


    //pi^00 = u_i u_j pi^{ij}/gamma^2
    shv(0,0) =( u_cov(0)*u_cov(0)*shv(1,1) + u_cov(1)*u_cov(1)*shv(2,2)
                + 2.*u_cov(1)*u_cov(0)*shv(2,1) + 2.*u_cov(2)*u_cov(0)*shv(3,1) 
                + 2.*u_cov(2)*u_cov(1)*shv(3,2)
                - u_cov(2)*u_cov(2)*(shv(1,1)+shv(2,2))/t2
              )/(gamma*gamma-u_cov(2)*u_cov(2)/(t2));


    //pi^33 = (pi^00 - pi^11 - pi^22)/t^2
    shv(3,3) = (shv(0,0) - shv(1,1) - shv(2,2))/t2;

    //Return shv
    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir) = shv(idir,jdir) ;

 
  };
  Cabana::simd_parallel_for(simd_policy,kokkos_ensure_consistency,"kokkos_ensure_consistency");
  Kokkos::fence();
}


///@brief Specialization for (2+1)D 
template<>
void EoM_default<2>::reset_pi_tensor(std::shared_ptr<SystemState<2>> sysPtr)
{
  CREATE_VIEW(device_,sysPtr->cabana_particles)
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());
  double t =  (sysPtr->t);
  double t2 = t*t;
  auto kokkos_ensure_consistency = KOKKOS_LAMBDA(const int is, int ia) 
  {
    milne::Vector<double,2> u;
    for(int idir=0; idir<2; ++idir) u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
    milne::Vector<double,2> u_cov = u;
    u_cov.make_covariant(t2);

    double gamma = device_hydro_scalar.access(is,ia, hydro_info::gamma);
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);

    //compute bulk from big bulk 
    double big_bulk = device_hydro_scalar.access(is, ia, hydro_info::bigBulk);
    device_hydro_scalar.access(is, ia, hydro_info::bulk) = big_bulk*sigma;
    //Declare caches
    milne::Matrix<double,4,4> shv;

    //compute shear from big shear
    for(int idir=0; idir<2; ++idir){
    for(int jdir=idir; jdir<3; ++jdir){
      shv(idir+1,jdir+1) = device_hydro_shear_aux_vector.access(is, ia, hydro_info::bigshv, idir, jdir)*sigma;
     // std::cout << "bigshv: " << device_hydro_spacetime_matrix.access(is, ia, hydro_info::bigshv, idir, jdir) << std::endl;
      }}

    //Symmetrizes space part
    for( int i=1; i<4; i++ )
    for( int j=i+1; j<4; j++ )
      shv(j,i) = shv(i,j);

    //pi^{0i} = -\pi^{ij}u_j/gamma 
    for(int idir=1; idir<3; ++idir){
      milne::Vector<double,2> shv_i;
      for(int jdir=1; jdir<3; ++jdir) shv_i(jdir-1) = shv(idir,jdir);
      shv(0,idir) = -1.*milne::contract(u_cov,shv_i)/gamma; 
    }
    //enforce 0 for other components
    for(int idir=3; idir<4; ++idir) shv(0,idir) = 0.0;

    //Symmetrizes time part
    for( int i=1; i<4; i++ )
      shv(i,0) = shv(0,i);


    //pi^00 = u_i u_j pi^{ij}/gamma^2
    shv(0,0) =( u_cov(0)*u_cov(0)*shv(1,1) + u_cov(1)*u_cov(1)*shv(2,2)
                + 2.*u_cov(1)*u_cov(0)*shv(1,2))/(gamma*gamma);


    //pi^33 = (pi^00 - pi^11 - pi^22)/t^2
    shv(3,3) = (shv(0,0) - shv(1,1) - shv(2,2))/(t2);

    //Return shv
    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir) = shv(idir,jdir) ;

  };
  Cabana::simd_parallel_for(simd_policy,kokkos_ensure_consistency,"kokkos_ensure_consistency");
  Kokkos::fence();
}


/// @brief Specialization for (1+1)D
template<>
void EoM_default<1>::reset_pi_tensor(std::shared_ptr<SystemState<1>> sysPtr)
{
  double t =  (sysPtr->t);
  double t2 = t*t;
  CREATE_VIEW(device_,sysPtr->cabana_particles)
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());

  auto kokkos_ensure_consistency = KOKKOS_LAMBDA(const int is, int ia) 
  {
    milne::Vector<double,1> u;
    for(int idir=0; idir<1; ++idir) u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
    milne::Vector<double,1> u_cov = u;
    u_cov.make_covariant(t2); //Transforms u^i to -u_i
    double gamma = device_hydro_scalar.access(is,ia, hydro_info::gamma);
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    //compute bulk from big bulk 
    double big_bulk = device_hydro_scalar.access(is, ia, hydro_info::bigBulk);
    device_hydro_scalar.access(is, ia, hydro_info::bulk) = big_bulk*sigma;
    //Declare caches

    milne::Matrix<double,4,4> shv;


    //compute shear from big shear
    for(int idir=0; idir<2; ++idir){
    for(int jdir=idir; jdir<3; ++jdir){
      shv(idir+1,jdir+1) = device_hydro_shear_aux_vector.access(is, ia, hydro_info::bigshv, idir, jdir)*sigma;
      //std::cout << "bigshv: " << device_hydro_spacetime_matrix.access(is, ia, hydro_info::bigshv, idir, jdir) << std::endl;
    }
    }

    //Symmetrizes space part
    for( int i=1; i<4; i++ ){
    for( int j=i+1; j<4; j++ ){
      shv(j,i) = shv(i,j);
    } 
    }

    //pi^{0i} = -\pi^{ij}u_j/gamma 
    for(int idir=1; idir<2; ++idir){
      milne::Vector<double,1> shv_i;
      for(int jdir=1; jdir<2; ++jdir) shv_i(jdir-1) = shv(idir,jdir);
      shv(0,idir) = -1.*milne::contract(u_cov,shv_i)/gamma; 
    }
    //enforce 0 for other components
    for(int idir=2; idir<4; ++idir){
      shv(0,idir) = 0.0;
    } 

    //Symmetrizes time part
    for( int i=1; i<4; i++ ){
      shv(i,0) = shv(0,i);
    }


    //pi^00 = u_i u_j pi^{ij}/gamma^2
    shv(0,0) =(- u_cov(0)*u_cov(0)*(shv(1,1)+shv(2,2)))/(gamma*gamma*t2-u(0)*u(0));


    //pi^33 = (pi^00 - pi^11 - pi^22)/t^2
    shv(3,3) = (shv(0,0) - shv(1,1) - shv(2,2))/(t2);

    //Return shv
    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir) = shv(idir,jdir) ;


  };
  Cabana::simd_parallel_for(simd_policy,kokkos_ensure_consistency,"kokkos_ensure_consistency");
  Kokkos::fence();
}




/// @brief Calculates the Lorentz contraction factor \f$ \gamma = u^0\f$.
/// @details Calculates the Lorentz contraction factor \f$ \gamma \f$ in the
/// general case. This assumes that the last component of the velocity vector is
/// the longitudinal velocity. The expression is
/// \f[ \gamma = \sqrt{1 + \sum_{i=1}^{D-1} u_i^2 + \tau^2 u_D^2} \f]
/// @param u The velocity vector (space components only).
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return The value of \f$ \gamma \f$.
template<unsigned int D> KOKKOS_FUNCTION
double EoM_default<D>::gamma_calc(double u[D], const double &time_squared)
{
    double dot_u = 0;
    for (unsigned int i=0; i<D-1; i++)
      dot_u+= u[i]*u[i];
    dot_u += time_squared*u[D-1]*u[D-1];
    return sqrt(1.0+dot_u);
}

/// @brief Calculates the Lorentz contraction factor \f$ \gamma = u^0\f$ for 2D 
/// simulations.
/// @details Calculates the Lorenz contraction factor \f$ \gamma \f$. The 
/// expression is
/// \f[ \gamma = \sqrt{1 + \sum_{i=1}^{D-1} u_i^2 \f]
/// @param u The velocity vector (space components only).
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return The value of \f$ \gamma \f$.
template<> KOKKOS_FUNCTION
double EoM_default<2>::gamma_calc(double u[2], const double &time_squared)
{
    double dot_u = 0;
    for (unsigned int i=0; i<2; i++)
      dot_u += u[i]*u[i];
    return sqrt(1.0+dot_u);
}
/// @brief Transforms a scalar from the lab frame to the Local Rest Frame (LRF).
/// @details The expression implemented here is \f[ \frac{lab}{\gamma \tau} \f]
/// @tparam D The number of spatial dimensions.
/// @param lab The quantity in the lab (computational) frame.
/// @param gamma Lorentz contraction factor.
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return The quantity in the fluid LRF (local rest frame).
template<unsigned int D> KOKKOS_FUNCTION
double EoM_default<D>::get_LRF(const double &lab, const double &gamma,
                               const double &t)
{
                                return lab/gamma/t;
}
template<unsigned int D>
void EoM_default<D>::evaluate_time_derivatives( std::shared_ptr<SystemState<D>> sysPtr, std::shared_ptr<Settings> settingsPtr)
{
  // #ifdef DEBUG
  // ofstream outfile;
  // outfile.open("gradients.dat");
  // #endif
  double t = (sysPtr->t);
  double t2 = t*t;
  milne::Vector<double,D> delta_i_eta = milne::delta_i_eta<D>();
  CREATE_VIEW(device_,sysPtr->cabana_particles);
  bool using_shear = settingsPtr->using_shear; 
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());
  if(using_shear){
    calculate_MRF_shear(sysPtr);
    //std::cout << "Hi shear" << std::endl;
  }
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
    double muB = device_thermo.access(is, ia, thermo_info::muB);
    double muQ = device_thermo.access(is, ia, thermo_info::muQ);
    double muS = device_thermo.access(is, ia, thermo_info::muS);
    double sigma_star = device_hydro_scalar.access(is, ia, hydro_info::sigma_star);
    double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);
    double zeta = device_hydro_scalar.access(is, ia, hydro_info::zeta_Pi);
    double tau_Pi = device_hydro_scalar.access(is, ia, hydro_info::tau_Pi);
    double tau_pi = device_hydro_scalar.access(is, ia, hydro_info::tau_pi );
    double delta_PiPi = device_hydro_scalar.access(is, ia, hydro_info::delta_PiPi);
    double a = device_hydro_scalar.access(is, ia, hydro_info::a);
    double F_big_bulk = device_hydro_scalar.access(is, ia, hydro_info::F_big_bulk);
    double F_shv_nabla_u = device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u);
    double F_big_entropy = device_hydro_scalar.access(is, ia, hydro_info::F_big_entropy);
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    double phi1 = device_hydro_scalar.access(is, ia, hydro_info::phi1);
    double phi3 = device_hydro_scalar.access(is, ia, hydro_info::phi3);
    double lambda_Pipi = device_hydro_scalar.access(is, ia, hydro_info::lambda_Pipi);
    //auxiliary zeta tilde to control IR or DNMR
    double zeta_tilde  = zeta +  a*(delta_PiPi - tau_Pi)*bulk;
    
    double rhoQ_ext = 0.;
    double rhoS_ext = 0.;
    double rhoB_ext = 0.;
    //auxiliary array to store muB, muS, muQ
    milne::Vector<double,3> mu_vec = {muB, muS, muQ};
    //auxiliary array to store rho_ext
    milne::Vector<double,3> rho_ext = {rhoB_ext, rhoS_ext, rhoQ_ext};

    
    double j0_ext =0.;
    milne::Vector<double,D> j_ext;
    for(int idir=0; idir<D; ++idir){
      j_ext(idir) = 0.;
    }

    milne::Vector<double,D> v, u, grad_u0;
    milne::Vector<double,D> M_big_bulk_aux, M_shv_nabla_u, M_big_entropy;
    milne::Vector<double,3> R_big_entropy, R_big_bulk, F_big_N;
    milne::Matrix<double,D,D> gradV, grad_uj;
    milne::Matrix<double,3,3> R_big_N;
    milne::Matrix<double,4,4> shv;
    for(int idir=0; idir<D; ++idir){  
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
    }
    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
    milne::Vector<double,D> u_cov = u; 
    u_cov.make_covariant(t2); 
    double geometric_factor = delta_i_eta(D-1)*u(D-1)*u(D-1)*t/gamma;

    milne::Matrix<double,4,4> shv_hybrid = shv;
    //shv_hybrid = \pi^\mu_\nu
    shv_hybrid.make_covariant(1, t2); 

    milne::Matrix<double,4,4> shv_cov = shv_hybrid; 
    //shv_cov = \pi_{\mu\nu}
    shv_cov.make_covariant(0, t2);  

    double divV = 0;
    for(int idir=0; idir<D; ++idir){
      divV += device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, idir);
    for(int jdir=0; jdir<D; ++jdir){
      gradV(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, jdir);
      }
    }
    
    //fill M matrices
    for(int idir=0; idir<D; ++idir){
      M_big_bulk_aux(idir) = (zeta_tilde*u_cov(idir)/gamma)/(sigma*gamma*tau_Pi);
      M_big_entropy(idir) = (bulk*u_cov(idir)/gamma)/(sigma*gamma*T);
    };
    //fill R matrices
    for(int icharge=0; icharge<3; ++icharge){
      R_big_entropy(icharge) = -gamma*sigma*mu_vec(icharge);
      
      R_big_bulk(icharge) = 0.0;
      for(int jcharge=0; jcharge<3; ++jcharge){
        R_big_N(icharge,jcharge) = 0.0;
      }
      R_big_N(icharge,icharge) += gamma*sigma;
    }
    
    
  
    //fill vectors
    F_big_bulk = -(zeta_tilde*(gamma*divV + gamma/t + geometric_factor)
                  +bulk  - a*phi1*bulk*bulk);
    F_big_entropy = ( -bulk*(gamma*divV + gamma/t + geometric_factor)
                    + gamma*j0_ext
                    +milne::contract(u_cov,j_ext))/(sigma*gamma*T);
    
    if(using_shear){ 
      //std::cout << "Hi shear" << std::endl;
      double F_shv_nabla_u = device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u);
      F_big_entropy += F_shv_nabla_u/(sigma*gamma*T);
      F_big_bulk += -(- a*lambda_Pipi*F_shv_nabla_u
                      - a*phi3*milne::contract(shv_cov,shv))/(sigma*gamma*tau_Pi);
      milne::Vector<double,D> M_shv_nabla_u;
      for(int idir=0; idir<D; ++idir){
        M_shv_nabla_u(idir) =  device_hydro_vector.access(is, ia, hydro_info::M_shv_nabla_u,idir);
        M_big_bulk_aux(idir) += a*lambda_Pipi*M_shv_nabla_u(idir)/(sigma*gamma*tau_Pi);
        M_big_entropy(idir) += M_shv_nabla_u(idir)/(sigma*gamma*T);
        
      }

    }

    for(int icharge=0; icharge<3; ++icharge){
      F_big_N(icharge) = rho_ext(icharge);
    }
    //stores the results
    device_hydro_scalar.access(is, ia, hydro_info::F_big_bulk) = F_big_bulk;
    device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u) = F_shv_nabla_u;
    device_hydro_scalar.access(is, ia, hydro_info::F_big_entropy) = F_big_entropy;
    for(int idir=0; idir<D; ++idir){
      device_hydro_vector.access(is, ia, hydro_info::M_big_bulk, idir) = M_big_bulk_aux(idir);
      device_hydro_vector.access(is, ia, hydro_info::M_big_entropy, idir) = M_big_entropy(idir);
    }
    for(int icharge=0; icharge<3; ++icharge){
      device_hydro_vector.access(is, ia, hydro_info::R_big_entropy, icharge) = R_big_entropy(icharge);
      device_hydro_vector.access(is, ia, hydro_info::F_big_N, icharge) = F_big_N(icharge);
      device_hydro_vector.access(is, ia, hydro_info::R_big_bulk,icharge) = R_big_bulk(icharge);
      for(int jcharge=0; jcharge<3; ++jcharge){
        device_hydro_space_matrix.access(is, ia, hydro_info::R_big_N, icharge, jcharge) = R_big_N(icharge,jcharge);
      }
    }
   
  };

  Cabana::simd_parallel_for(simd_policy, fill_auxiliary_variables, "fill_auxiliary_variables");
  Kokkos::fence();
  //calculate du/dt
  auto compute_velocity_derivative = KOKKOS_LAMBDA(const int is, const int ia)
  { 

      double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
      double bigBulk = device_hydro_scalar.access(is, ia,hydro_info::bigBulk );
      double bulk = device_hydro_scalar.access(is, ia,hydro_info::bulk );
      double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
      double w = device_thermo.access(is, ia, thermo_info::w);
      double s = device_thermo.access(is, ia, thermo_info::s);
      double dwds = device_thermo.access(is, ia, thermo_info::dwds);
      double rhob = device_thermo.access(is, ia, thermo_info::rhoB);
      double rhos = device_thermo.access(is, ia, thermo_info::rhoS);
      double rhoq = device_thermo.access(is, ia, thermo_info::rhoQ);
      double dwdrhoB = device_thermo.access(is, ia, thermo_info::dwdB);
      double dwdrhoS = device_thermo.access(is, ia, thermo_info::dwdS); 
      double dwdrhoQ = device_thermo.access(is, ia, thermo_info::dwdQ);
      double F_big_bulk = device_hydro_scalar.access(is, ia, hydro_info::F_big_bulk);
      double F_big_entropy = device_hydro_scalar.access(is, ia, hydro_info::F_big_entropy);
      double tau_Pi = device_hydro_scalar.access(is, ia, hydro_info::tau_Pi);

      double divV = 0;
      milne::Vector<double,D> u, gradshear, gradP, gradBulk, divshear;
      milne::Vector<double,D> M_bulk, M_S, M_big_bulk;
      milne::Vector<double,D> F_0i_shear, j_ext; 
      milne::Matrix<double,D,D>  M_0i_shear;
      milne::Vector<double,3> R_S,R_big_bulk, rhoVec, dwdrhoVec, F_big_N;
      milne::Matrix<double,D,3> R_0i_shear,M_big_N;
      milne::Matrix<double,4,4> shv;
      for(int idir=0; idir<D; ++idir){
        u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
        gradP(idir) = device_hydro_vector.access(is, ia, hydro_info::gradP, idir);
        gradBulk(idir) = device_hydro_vector.access(is, ia, hydro_info::gradBulk, idir);
        divshear(idir) = device_hydro_vector.access(is, ia, hydro_info::divshear, idir);
        gradshear(idir) = device_hydro_vector.access(is, ia, hydro_info::gradshear, idir);
        M_big_bulk(idir) = device_hydro_vector.access(is, ia, hydro_info::M_big_bulk, idir);
        M_S(idir) = device_hydro_vector.access(is, ia, hydro_info::M_big_entropy, idir);
        F_0i_shear(idir) = device_hydro_vector.access(is, ia, hydro_info::F_0i_shear, idir);
        j_ext(idir) = device_hydro_vector.access(is, ia, hydro_info::j_ext, idir);
        divV += device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, idir);
        for(int  jdir=0; jdir<D; ++jdir){
          M_0i_shear(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::M_0i_shear, idir, jdir);
        }
      }
      for(int idir=0; idir<4; ++idir)
      for(int jdir=0; jdir<4; ++jdir)
        shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
      for(int icharge=0; icharge<3; ++icharge){
        R_S(icharge) = device_hydro_vector.access(is, ia, hydro_info::R_big_entropy, icharge);
        R_big_bulk(icharge) = device_hydro_vector.access(is, ia, hydro_info::R_big_bulk, icharge);
        F_big_N(icharge) = device_hydro_vector.access(is, ia, hydro_info::F_big_N, icharge);
        for(int idir=0; idir<D; ++idir){
          R_0i_shear(idir,icharge) = device_hydro_space_matrix.access(is, ia, hydro_info::R_0i_shear, idir, icharge);
          M_big_N(idir,icharge) = device_hydro_space_matrix.access(is, ia, hydro_info::M_big_N, idir, icharge);
        }
      }
      rhoVec = {rhob, rhos, rhoq};
      dwdrhoVec = {dwdrhoB, dwdrhoS, dwdrhoQ};
      milne:: Vector<double,D> gradP_contra, gradBulk_contra, u_cov;
      u_cov = u;
      gradP_contra = gradP;
      gradBulk_contra = gradBulk;
      gradP_contra.make_contravariant(t2);
      gradBulk_contra.make_contravariant(t2);
      u_cov.make_covariant(t2);

      double geometric_factor = delta_i_eta(D-1)*u(D-1)*u(D-1)*t/gamma;   
      
      milne::Vector<double,D> M_w = (s*dwds/gamma
                                    +milne::contract(rhoVec,dwdrhoVec)/gamma)*u_cov/gamma
                                    +sigma*dwds*M_S;
      milne::Vector<double,3> R_w = sigma*dwdrhoVec 
                                    +sigma*dwds*R_S;
      double F_w = -(s*dwds/gamma
                    +milne::contract(rhoVec,dwdrhoVec))*(geometric_factor+ gamma/t + gamma*divV)
                    +sigma*dwds*F_big_entropy;
      double F_bulk = F_big_bulk - bulk*(gamma*divV + gamma/t + geometric_factor)/gamma;
      M_bulk = M_big_bulk + bulk*u_cov/(gamma*gamma*tau_Pi);                         



      // set the Mass and the Force
      milne::Matrix <double,D,D> M_u;
      for(int idir=0; idir<D; ++idir){
        for(int jdir=0; jdir<D; ++jdir){
          M_u(idir,jdir) = u(idir)*gamma*(M_big_bulk(jdir) + M_w(jdir))
                          -u(idir)*u_cov(jdir)*(w+bulk)/gamma;
        }
        M_u(idir,idir) += gamma*(w+bulk);
      }
      milne::Vector<double,D> F_u = j_ext + (gradP_contra + gradBulk_contra)
                                  -gamma*(F_w+F_bulk)*u
                                  -(w+bulk)*(gamma*divV + gamma/t + geometric_factor)*u
                                  -2.*(w+bulk)*gamma*u(D-1)*delta_i_eta/t;
      
                                 
      milne::Matrix<double,D,3> R_u;
      for(int idir=0; idir<D; ++idir){
        for(int icharge=0; icharge<3; ++icharge){
          R_u(idir,icharge) = -u(idir)*gamma*(R_w(icharge) + R_big_bulk(icharge));
        }
      }

      if ( using_shear )
      { //std::cout << "Hi shear" << std::endl;
        for(int idir=0; idir<D; ++idir){
          F_u(idir) +=  -shv(idir+1,0)/t - 2.*delta_i_eta(idir)*shv(0,3)/t 
                        -divshear(idir) + gradshear(idir) -F_0i_shear(idir);
          for(int icharge=0; icharge<3; ++icharge){
            R_u(idir,icharge) += -R_0i_shear(idir,icharge);
          }
          for(int jdir=0; jdir<D; ++jdir){
            M_u(idir,jdir) += M_0i_shear(idir,jdir);
            //stores M_u in the space matrix
            device_hydro_space_matrix.access(is, ia, hydro_info::M_u, idir, jdir) = M_u(idir,jdir);
          }
          //stores F_u in the space vector
          device_hydro_vector.access(is, ia, hydro_info::F_u, idir) = F_u(idir);
        }
      }
    

     milne::Matrix<double,D,D> aux_MR;
     milne::Vector<double,D> aux_FR;
     milne::Matrix<double,3,3> R_big_N;
     for(int icharge=0; icharge<3; ++icharge){
       for(int jcharge=0; jcharge<3; ++jcharge){
         R_big_N(icharge,jcharge) = device_hydro_space_matrix.access(is, ia, hydro_info::R_big_N, icharge, jcharge);
       }
     }
     milne::Matrix<double,3,3> RI = milne::inverse(R_big_N);

     for(int idir=0; idir<D; ++idir){
       aux_FR(idir) = 0.0;
       for(int jdir=0; jdir<D; ++jdir){
          aux_MR(idir,jdir) = 0.0;
          for(int icharge=0; icharge<3; ++icharge){
            for(int jcharge=0; jcharge<3; ++jcharge){
              aux_MR(idir,jdir) += R_u(idir,icharge)*RI(icharge,jcharge)*M_big_N(jdir,jcharge);
            }
          }
          for(int icharge=0; icharge<3; ++icharge){
            for(int jcharge=0; jcharge<3; ++jcharge){
              aux_FR(idir) += R_u(idir,icharge)*RI(icharge,jcharge)*F_big_N(jcharge);
            }
          }
        }
      };

      /// @brief 
      milne::Matrix<double,D,D> MI = milne::inverse(M_u-aux_MR);
      milne::Vector<double,D> du_dt = MI*(F_u+
      aux_FR);
      for(int idir=0; idir<D; ++idir) device_hydro_vector.access(is,ia,hydro_info::du_dt, idir) = du_dt(idir);
  };
  Cabana::simd_parallel_for(simd_policy, compute_velocity_derivative, "compute_velocity_derivative");
  Kokkos::fence();
  //all derivatives and quantities that depends on du/dt
  auto compute_derivatives = KOKKOS_LAMBDA(const int is, const int ia)
  {


    double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double divV = 0;
    
    milne::Vector<double,D> v, du_dt, u,u_cov;
    milne::Vector<double,D> M_big_entropy, M_big_bulk;
    milne::Vector<double,3> F_big_N, R_big_entropy, R_big_bulk;
    for(int idir=0; idir<D; ++idir){
      v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
      du_dt(idir) = device_hydro_vector.access(is, ia, hydro_info::du_dt, idir);
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      divV += device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, idir);
      M_big_entropy(idir) = device_hydro_vector.access(is, ia, hydro_info::M_big_entropy, idir);
      M_big_bulk(idir) = device_hydro_vector.access(is, ia, hydro_info::M_big_bulk, idir);
    }
    for(int icharge=0; icharge<3; ++icharge){
      F_big_N(icharge) = device_hydro_vector.access(is, ia, hydro_info::F_big_N, icharge);
      R_big_entropy(icharge) = device_hydro_vector.access(is, ia, hydro_info::R_big_entropy, icharge);
      R_big_bulk(icharge) = device_hydro_vector.access(is, ia, hydro_info::R_big_bulk, icharge);
    }
    milne::Matrix<double,4,4> shv;
    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      shv(idir, jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
    double geometric_factor = delta_i_eta(D-1)*u(D-1)*u(D-1)*t/gamma;
    u_cov = u;
    u_cov.make_covariant(t2);

    milne::Vector<double,3> MU_aux;
    for(int icharge=0; icharge<3; ++icharge){
      MU_aux(icharge) = 0.;
      for(int idir=0; idir<D; ++idir){
        MU_aux(icharge) +=device_hydro_space_matrix.access(is, ia, hydro_info::M_big_N, idir, icharge)*du_dt(idir);
      }
    }

    milne::Matrix<double,3,3> R_N_inv;
    for(int icharge=0; icharge<3; ++icharge){
      for(int jcharge=0; jcharge<3; ++jcharge){
        R_N_inv(icharge,jcharge) = device_hydro_space_matrix.access(is, ia, hydro_info::R_big_N, icharge, jcharge);
      }
    }
    R_N_inv = milne::inverse(R_N_inv); 
    milne::Vector<double,3> dN_dt = R_N_inv*(F_big_N + MU_aux);

    // time derivative of ``specific entropy density per particle"
    double d_dt_specific_s = milne::contract(M_big_entropy,du_dt)
                            +device_hydro_scalar.access(is, ia, hydro_info::F_big_entropy)
                            +milne::contract(R_big_entropy,dN_dt);
	  //Bulk evolution equation
    double dbigBulk_dt = milne::contract(M_big_bulk,du_dt)
                    +device_hydro_scalar.access(is, ia, hydro_info::F_big_bulk)
                    +milne::contract(R_big_bulk,dN_dt);
    //shv_nabla_u evolution equation
    milne::Vector<double,D> M_shv_nabla_u; 
    for(int idir=0; idir<D; ++idir){
      M_shv_nabla_u(idir) = device_hydro_vector.access(is, ia, hydro_info::M_shv_nabla_u, idir);
    }
    double shv_nabla_u = milne::contract(M_shv_nabla_u,du_dt)
                        +device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u);
    //theta equation
    double theta = -1.*milne::contract(u_cov,du_dt)/gamma + gamma*divV + gamma/t + geometric_factor;

    device_hydro_scalar.access(is, ia, hydro_info::dbigBulk_dt) = dbigBulk_dt;
    device_hydro_scalar.access(is, ia, hydro_info::theta) = theta;
    device_hydro_scalar.access(is, ia, hydro_info::shv_nabla_u) = shv_nabla_u;
    device_d_dt_spec.access(is, ia, densities_info::s) = d_dt_specific_s;
    device_d_dt_spec.access(is, ia, densities_info::rhoB) = dN_dt(0);
    device_d_dt_spec.access(is, ia, densities_info::rhoS) = dN_dt(1);
    device_d_dt_spec.access(is, ia, densities_info::rhoQ) = dN_dt(2);
    
    if(using_shear){
      //std::cout << "Hi shear" << std::endl;
      milne::Matrix<double,2,3> d_bigshv_dt;
      for(int idir=0; idir<2; ++idir){
        for(int jdir=idir; jdir<3; ++jdir){
          double M_du_aux = 0.;
          double R_dn_aux = 0.;
          for(int kdir=0; kdir<D; ++kdir){
            int linear_index_M = jdir * D + kdir;
            M_du_aux += device_hydro_shear_aux_matrix.access(is, ia, hydro_info::M_big_shear, idir,linear_index_M)*du_dt(kdir);
          }
          for(int icharge=0; icharge<3; ++icharge){
            int linear_index_R = jdir * 3 + icharge;
            R_dn_aux += device_hydro_shear_aux_matrix.access(is, ia, hydro_info::R_big_shear, idir, linear_index_R)*dN_dt(icharge);
          }
          d_bigshv_dt(idir,jdir) = M_du_aux 
                                  + device_hydro_shear_aux_vector.access(is, ia, hydro_info::F_big_shear, idir, jdir);
           device_hydro_shear_aux_vector.access(is, ia, hydro_info::d_bigshv_dt, idir,jdir) = d_bigshv_dt(idir,jdir);
        }
      }
    }
    else{
      for(int idir=0; idir<2; ++idir){
        for(int jdir=idir; jdir<3; ++jdir){
          device_hydro_shear_aux_vector.access(is, ia, hydro_info::d_bigshv_dt, idir,jdir) = 0.0;
        }
      }
    }

	  //formulating simple setup for Beta_Bulk derivative
    ///TODO: Implement this
    //hi.finite_diff_cs2   =  (ti.cs2 - hi.prev_cs2)/0.05; // Asadek
	  //hi.finite_diff_T   =  (ti.T - hi.prev_T)/0.05; // Asadek
	  //hi.finite_diff_w   =  (ti.w - hi.prev_w)/0.05; // Asadek
	  //hi.dBeta_dt      = 0.5*((-hi.finite_diff_T/(ti.T*ti.T))*(1/ti.w)*(1/((1/3-ti.cs2)*(1/3-ti.cs2))))
	  //                 + 0.5*((-hi.finite_diff_w/(ti.w*ti.w))*(1/ti.T)*(1/((1/3-ti.cs2)*(1/3-ti.cs2))))
		//			   + 0.5*((4*ti.cs2*hi.finite_diff_cs2*(1/((1/3-ti.cs2)*(1/3-ti.cs2)*(1/3-ti.cs2))))*(1/ti.T)*(1/ti.w));//Asadek 

  };
  Cabana::simd_parallel_for(simd_policy, compute_derivatives, "compute_derivatives");
  Kokkos::fence();

      // computing dEz_dt
  auto compute_Ez_derivative = KOKKOS_LAMBDA(const int is, const int ia)
    {
      double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);
      double shv33 = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv,3,3);
      double sigma_star = device_hydro_scalar.access(is, ia, hydro_info::sigma_star);
      double norm_spec = device_norm_spec.access(is, ia, densities_info::s);
      double p = device_thermo.access(is, ia, thermo_info::p);  
      double e = device_smoothed.access(is, ia, densities_info::e);
      milne::Vector<double,D> u;
      for(int idir=0; idir<D; ++idir){
        u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      }
      double u3 = 0;
      if (D==1){
        u3 = u(0); // eta flow velocity
      }
      else if (D == 3){
        u3 = u(2);  
      }
      else if (D == 2){
        u3 = 0;
      }

      double dEz = 0;
      dEz = ((e + p + bulk)*u3*u3 + (p + bulk)/t2 + shv33 ) * norm_spec * t2 / sigma_star;
      device_contribution_to_total_dEz.access(is, ia) = dEz;
    };
    Cabana::simd_parallel_for(simd_policy, compute_Ez_derivative, "compute_Ez_derivative");
    Kokkos::fence();


  };

/// @brief Calculate the M and R matrices and F vectors for the 
/// shear viscosity. Needs specialization for each dimension.
/// @details This function calculates the M and R matrices and the F vectors
/// for the shear viscosity, which is used in the shear viscosity derivative 
/// calculation as well as in the acceleration equation.
/// @tparam D The number of spatial dimensions.
/// @param sysPtr A pointer to the object of class SystemState
template<unsigned int D>
void EoM_default<D>::calculate_MRF_shear(std::shared_ptr<SystemState<D>> sysPtr){};

template<>
void EoM_default<3>::calculate_MRF_shear(std::shared_ptr<SystemState<3>> sysPtr)
{
  double t = (sysPtr->t);
  double t2 = t*t;
  milne::Vector<double,3> delta_i_eta = milne::delta_i_eta<3>();
  CREATE_VIEW(device_,sysPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());

  auto compute_MRF_shear = KOKKOS_LAMBDA(const int is, const int ia)
  {
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double eta_pi = device_hydro_scalar.access(is, ia, hydro_info::eta_pi);
    double phi6 = device_hydro_scalar.access(is, ia, hydro_info::phi6);
    double phi7 = device_hydro_scalar.access(is, ia, hydro_info::phi7);
    double lambda_piPi = device_hydro_scalar.access(is, ia, hydro_info::lambda_piPi);
    double tau_pipi = device_hydro_scalar.access(is, ia, hydro_info::tau_pipi);
    double a= device_hydro_scalar.access(is, ia, hydro_info::a);
    double delta_pipi = device_hydro_scalar.access(is, ia, hydro_info::delta_pipi);
    double tau_pi = device_hydro_scalar.access(is, ia, hydro_info::tau_pi);
    double tilde_delta = delta_pipi - tau_pi;
    double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);

    milne::Vector<double,3> u, u_cov, M_shv_nabla_u ,F_0i_shear;
    milne::Vector<double,3> v, grad_u0;
    milne::Vector<double,4> four_u;
    milne::Matrix<double,2,3> F_big_shear;
    milne::Matrix<double,3,3> gradV, grad_uj, M_0i_shear;
    milne::Matrix<double,4,4> shv_hybrid, shv, shv_cov;
    milne::Matrix3D<double,2,3,3> M_big_shear;
    for(int idir=0; idir<3; ++idir){
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      four_u(idir+1) = u(idir);
    }
    four_u(0) = gamma;
    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
    u_cov.make_covariant(t2);
    double geometric_factor = delta_i_eta(2)*u(2)*u(2)*t/gamma;
    shv_hybrid = shv;
    shv_hybrid.make_covariant(1, t2);
    shv_cov = shv_hybrid;
    shv_cov.make_covariant(0, t2);

    milne::Vector<double,3> F_i0_sigma;
    milne::Vector<double,3> F_i0_D;
    milne::Vector<double,3> F_i0_domega;
    milne::Vector<double,3> F_i0_dsigma;
    milne::Vector<double,3> F_i0_dd;
    milne::Matrix<double,3,3> M_i0_sigma;
    milne::Matrix<double,3,3> M_i0_D;
    milne::Matrix<double,3,3> M_i0_domega;
    milne::Matrix<double,3,3> M_i0_dsigma;
    milne::Matrix<double,3,3> M_i0_dd;



    double divV = 0;
    for(int idir=0; idir<3; ++idir){
      divV += device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, idir);
    for(int jdir=0; jdir<3; ++jdir){
      gradV(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, jdir);
      }
    }

    grad_uj = gamma*gradV+pow(gamma,3.)*(v*(v*gradV));
    grad_u0 = -1.*milne::contract(grad_uj,u_cov,milne::SecondIndex())/gamma;
    milne::Vector<double,3> grad_u0_contra = grad_u0;
    grad_u0_contra.make_contravariant(t2);

    milne::Matrix<double,3,3> contra_grad_uj;
    contra_grad_uj = grad_uj;
    contra_grad_uj.make_contravariant(0, t2);
    //calculate by hand due to shear dimensions being fixed,
    //and to avoid unnecessary specializations
    //milne::Vector<double,3> shv_hyb_i0  = milne::colp1(0, shv_hybrid);
    //milne::Matrix<double,3,3> vi_shv_hyb_0j = milne::outer(v,milne::rowp1(0,shv_hybrid));
    //calculate by hand due to shear dimensions being fixed,
    //and to avoid unnecessary specializations
    milne::Vector<double,3> shv_hyb_i0;
    for(int idir=0; idir<3; ++idir){
      shv_hyb_i0(idir) = shv_hybrid(idir+1,0);
    }
    
    milne::Matrix<double,3,3> vi_shv_hyb_0j;
    for(int idir=0; idir<3; ++idir){
      for(int jdir=0; jdir<3; ++jdir){
        vi_shv_hyb_0j(idir,jdir) =shv_hybrid(idir+1,jdir+1) -  v(idir)*shv_hybrid(0,jdir+1);
      }
    }    

    milne::Matrix<double,3,3> R_0i_shear;
    for(int idir=0; idir<3; ++idir){
      for(int icharge=0; icharge<3; ++icharge){
        R_0i_shear(idir,icharge) = 0.0;
      }
    }


    //aux vectors and matrices
    for(int idir=0; idir<3; ++idir){
      F_i0_sigma(idir) = -(milne::contract(grad_uj,v,milne::FirstIndex())(idir)
                          -grad_u0_contra(idir))/2.
                          +u(idir)*gamma*(geometric_factor
                          +gamma/t + gamma*divV)/3.
                          -gamma*geometric_factor*u(idir)/2.;     
      F_i0_D(idir) = gamma*(u(idir)*shv_hybrid(0,0)+gamma*shv_hybrid(idir+1,0))*geometric_factor
                     +t*u(2)*shv(3,idir+1);
      F_i0_domega(idir) = 0.0;
      F_i0_dsigma(idir) = 0.0;
      F_i0_dd(idir) = shv(idir+1,0)*(geometric_factor+gamma/t+gamma*divV);
      for(int jdir=0; jdir<3; ++jdir){
        M_i0_sigma(idir,jdir) = u(idir)*u_cov(jdir)/3.;
        M_i0_D(idir,jdir) = gamma*u(idir)*(shv_hybrid(0,jdir+1)-shv_hybrid(0,0)*u_cov(jdir)/gamma)
                            +gamma*gamma*(shv_hybrid(idir+1,jdir+1)-shv_hybrid(idir+1,0)*u_cov(jdir)/gamma);
        M_i0_domega(idir,jdir) = 0.0;
        M_i0_dsigma(idir,jdir) = 0.0;
        M_i0_dd(idir,jdir) = -shv(idir+1,0)*u_cov(jdir)/3.;
      }
      //diagonal terms
      M_i0_sigma(idir,idir) += (1.-gamma*gamma);

    }
    F_i0_sigma(2) += -(-u(2)*t-u(2)/t+gamma*u(2)*u(2)*t
                    +2.*gamma*gamma*u(2)/t)/2.;
    F_i0_D(2) += u(2)*shv(0,0)/t + gamma*shv(3,0)/t;

    //fill F vectors for shear viscosity
    milne::Vector<double,4> shear0mu_aux;
    for(int idir=0; idir<4; ++idir){
      shear0mu_aux(idir) = shv_hybrid(0,idir);
    }

    double F_shv_nabla_u = shv_hybrid(0,0)*geometric_factor
            +milne::contract(shv_hyb_i0 - shv_hybrid(0,0)*v,grad_u0)
            +milne::contract(vi_shv_hyb_0j,grad_uj)
            +shv_hybrid(3,3)*gamma/t
            + (t*shv_hybrid(3,0)+shv_hybrid(0,3)/t)*u(2)*delta_i_eta(2);
    //stores the results
    device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u) = F_shv_nabla_u;

    for(int idir=0; idir<3; ++idir){
      F_0i_shear(idir) = (-tilde_delta*F_i0_dd(idir)-tau_pi*F_i0_D(idir)
                          -a*2.*tau_pi*F_i0_domega(idir) 
                          -a*tau_pipi*F_i0_dsigma(idir)
                          +(2.*eta_pi+a*lambda_piPi*bulk)*F_i0_sigma(idir)
                          +a*phi6*bulk*shv(idir+1,0)
                          +a*phi7*(milne::contract(shv, shear0mu_aux, milne::SecondIndex())(idir+1)
                          +gamma*u(idir)*milne::contract(shv_cov,shv)/3.)/(tau_pi*gamma)
                          -shv(idir+1,0)*(geometric_factor+gamma/t+gamma*divV)/(gamma));
      //stores the results 
      device_hydro_vector.access(is, ia, hydro_info::F_0i_shear, idir) = F_0i_shear(idir);
    }    
    

    //fill M matrices for shear viscosity
    for(int idir=0; idir<3; ++idir){
      M_shv_nabla_u(idir) = shv_hybrid(0,idir+1) - shv_hybrid(0,0)*u_cov(idir)/gamma;
      for(int jdir=0; jdir<3; ++jdir){
        M_0i_shear(idir,jdir) = -(tau_pi*(M_i0_D(idir,jdir)
                                +tilde_delta*M_i0_dd(idir,jdir)/tau_pi
                                +2.*a*M_i0_domega(idir,jdir))
                                +a*tau_pipi*M_i0_dsigma(idir,jdir)
                                -(2.*eta_pi+a*lambda_piPi*bulk)*M_i0_sigma(idir,jdir))/(tau_pi*gamma)
                                +shv(idir+1,0)*u_cov(jdir)/(gamma*gamma);
        //stores the results
        device_hydro_space_matrix.access(is, ia, hydro_info::M_0i_shear, idir, jdir) = M_0i_shear(idir,jdir);
      }
      //stores the results
      device_hydro_vector.access(is, ia, hydro_info::M_shv_nabla_u, idir) = M_shv_nabla_u(idir);
    }


    
    milne::Vector<double,4> contra_metric_diag = milne::get_contra_metric_diagonal<3>(t);

    milne::Matrix3D<double,2,3,3> M_D;
    milne::Matrix3D<double,2,3,3> M_sigma;
    milne::Matrix3D<double,2,3,3> M_delta;
    milne::Matrix3D<double,2,3,3> M_domega;
    milne::Matrix3D<double,2,3,3> M_dsigma;

    milne::Matrix<double,2,3> F_D;
    milne::Matrix<double,2,3> F_sigma;
    milne::Matrix<double,2,3> F_delta;
    milne::Matrix<double,2,3> F_domega;
    milne::Matrix<double,2,3> F_dsigma;

    for(int idir=0; idir<2; ++idir){

      for(int jdir=idir; jdir<3; ++jdir){
        F_D(idir,jdir) = gamma*geometric_factor*(u(idir)*shv_hybrid(jdir+1,0)
                  +u(jdir)*shv_hybrid(idir+1,0));
                
        F_sigma(idir,jdir) = contra_grad_uj(idir,jdir)/2. + contra_grad_uj(jdir,idir)/2.
                  -u(idir)*u(jdir)*(geometric_factor+gamma/t+gamma*divV)/3.;         
        F_delta(idir,jdir) = shv(idir+1,jdir+1)*(geometric_factor+gamma/t+gamma*divV);
        F_domega(idir,jdir) = 0.0;
        F_dsigma(idir,jdir) = 0.0;                

        for(int kdir=0; kdir<3; ++kdir){
          M_D(idir,jdir,kdir) = gamma*(u(idir)*shv_hybrid(jdir+1,kdir+1)
                                      -u(idir)*shv_hybrid(jdir+1,0)*u_cov(kdir)/gamma
                                      +u(jdir)*shv_hybrid(idir+1,kdir+1)
                                      -u(jdir)*shv_hybrid(idir+1,0)*u_cov(kdir)/gamma);
          M_sigma(idir,jdir,kdir) = -2.*u(idir)*u(jdir)*u_cov(kdir)/(3.*gamma);
          M_delta(idir,jdir,kdir) = -shv(idir+1,jdir+1)*u_cov(kdir)/gamma;
          M_domega(idir,jdir,kdir) = 0.0;
          M_dsigma(idir,jdir,kdir) = 0.0;
        }
        //diagonal k terms
        M_sigma(idir,jdir,jdir) += -gamma*u(idir); 
        M_sigma(idir,jdir,idir) += -gamma*u(jdir);
      }
      //diagonal j terms
      for(int kdir=0; kdir<3; ++kdir){
        M_sigma(idir,idir,kdir) += 2.*contra_metric_diag(idir+1)*u_cov(kdir)/(3.*gamma);
      }

      //diagonal and metric terms for F
      F_D(idir,2) +=u(2)*shv(0,idir+1)/t + gamma*shv(3,idir+1)/t;
      F_sigma(idir,idir) += -contra_metric_diag(idir+1)*(geometric_factor+gamma/t+gamma*divV)/3.;
      F_sigma(idir,2) += -u(idir)*u(2)*u(2)*t/2. -u(idir)*u(2)*gamma/t;

    }

    for(int idir=0; idir<2; ++idir){
      for(int jdir=idir; jdir<3; ++jdir){
        F_big_shear(idir,jdir) = (-shv(idir+1,jdir+1) 
                                -tilde_delta*F_delta(idir,jdir)
                                -tau_pi*F_D(idir,jdir)
                                -a*2.*tau_pi*F_domega(idir,jdir)
                                -a*tau_pipi*F_dsigma(idir,jdir)
                                +(2.*eta_pi+a*lambda_piPi*bulk)*F_sigma(idir,jdir)
                                +a*phi6*bulk*shv(idir+1,jdir+1)
                                +a*phi7*(milne::contract(shv,shv_hybrid,milne::SecondIndex(),milne::SecondIndex())(idir,jdir)
                                +four_u(idir+1)*four_u(jdir+1)*milne::contract(shv,shv_cov)/3.))/(sigma*tau_pi*gamma);

        for(int kdir=0; kdir<3; ++kdir){
          
          M_big_shear(idir,jdir,kdir) = -(tau_pi*(M_D(idir,jdir,kdir)
                                      +tilde_delta*M_delta(idir,jdir,kdir)/tau_pi
                                      +2.*a*M_domega(idir,jdir,kdir))
                                      +a*tau_pipi*M_dsigma(idir,jdir,kdir)
                                      -(2.*eta_pi+a*lambda_piPi*bulk)*M_sigma(idir,jdir,kdir))/(sigma*tau_pi*gamma);
          //saves the results 
          int linear_index = jdir * 3 + kdir;
          device_hydro_shear_aux_matrix.access(is, ia, hydro_info::M_big_shear, idir, linear_index) = M_big_shear(idir,jdir,kdir);
        }
      }
      F_big_shear(idir,idir) += a*phi7*(-contra_metric_diag(idir+1)*milne::contract(shv,shv_cov))/(3.*sigma*tau_pi*gamma);
      //saves the results
      for(int jdir=idir; jdir<3; ++jdir){
        device_hydro_shear_aux_vector.access(is, ia, hydro_info::F_big_shear, idir, jdir) = F_big_shear(idir,jdir);
      }

                          
    }   
    //stores R matrix
    for(int idir=0; idir<3; ++idir){
      for(int icharge=0; icharge<3; ++icharge){
        device_hydro_space_matrix.access(is, ia, hydro_info::R_0i_shear, idir, icharge) = R_0i_shear(idir,icharge);
      }
    }
    for(int idir=0; idir<2; ++idir){
      for(int jdir=idir; jdir<3; ++jdir){
        for(int icharge=0; icharge<3; ++icharge){
          int linear_index = jdir * 3 + icharge;
          device_hydro_shear_aux_matrix.access(is, ia, hydro_info::R_big_shear, idir, linear_index) = 0.;
        }
      }
    } 
    
  };
  Cabana::simd_parallel_for(simd_policy, compute_MRF_shear, "compute_MRF_shear");
  Kokkos::fence();
};
    
template<>
void EoM_default<2>::calculate_MRF_shear(std::shared_ptr<SystemState<2>> sysPtr)
{
  double t = (sysPtr->t);
  double t2 = t*t;
  milne::Vector<double,2> delta_i_eta = milne::delta_i_eta<2>();
  CREATE_VIEW(device_,sysPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());

  auto compute_MRF_shear = KOKKOS_LAMBDA(const int is, const int ia)
  {
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double eta_pi = device_hydro_scalar.access(is, ia, hydro_info::eta_pi);
    double phi6 = device_hydro_scalar.access(is, ia, hydro_info::phi6);
    double phi7 = device_hydro_scalar.access(is, ia, hydro_info::phi7);
    double lambda_piPi = device_hydro_scalar.access(is, ia, hydro_info::lambda_piPi);
    double tau_pipi = device_hydro_scalar.access(is, ia, hydro_info::tau_pipi);
    //double a= device_hydro_scalar.access(is, ia, hydro_info::a);
    double a = 0.;
    double delta_pipi = device_hydro_scalar.access(is, ia, hydro_info::delta_pipi);
    double tau_pi = device_hydro_scalar.access(is, ia, hydro_info::tau_pi);
    double tilde_delta = delta_pipi - tau_pi;
    double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);

    milne::Vector<double,2> u, u_cov, M_shv_nabla_u ,F_0i_shear;
    milne::Vector<double,2> v, grad_u0;
    milne::Vector<double,4> four_u, four_u_cov;
    milne::Matrix<double,2,3> F_big_shear;
    milne::Matrix<double,2,2> gradV, grad_uj, M_0i_shear;
    milne::Matrix<double,4,4> shv_hybrid, shv, shv_cov;
    milne::Matrix3D<double,2,3,2> M_big_shear;
    for(int idir=0; idir<2; ++idir){
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      four_u(idir+1) = u(idir);
    }
    four_u(0) = gamma;
    //remaining dimensions
    four_u(3) = 0.0;
    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
    u_cov = u;
    u_cov.make_covariant(t2);
    for(int idir=0; idir<2; ++idir){
      four_u_cov(idir) = u_cov(idir);
    }
    four_u_cov(0) = gamma;
    four_u_cov(3) = 0.0;
    double geometric_factor = 0.;
    shv_hybrid = shv;
    shv_hybrid.make_covariant(1, t2);
    shv_cov = shv_hybrid;
    shv_cov.make_covariant(0, t2);

    milne::Vector<double,2> F_i0_sigma;
    milne::Vector<double,2> F_i0_D;
    milne::Vector<double,2> F_i0_domega;
    milne::Vector<double,2> F_i0_dsigma;
    milne::Vector<double,2> F_i0_dd;
    milne::Matrix<double,2,2> M_i0_sigma;
    milne::Matrix<double,2,2> M_i0_D;
    milne::Matrix<double,2,2> M_i0_domega;
    milne::Matrix<double,2,2> M_i0_dsigma;
    milne::Matrix<double,2,2> M_i0_dd;



    double divV = 0;
    for(int idir=0; idir<2; ++idir){
      divV += device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, idir);
    for(int jdir=0; jdir<2; ++jdir){
      gradV(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, jdir);
      }
    }

    grad_uj = gamma*gradV+pow(gamma,3.)*(v*(v*gradV));
    //grad_u0 = -1.*milne::contract(grad_uj,u_cov,milne::SecondIndex())/gamma;
    for (int idir=0; idir<2; ++idir){
      grad_u0(idir) = 0.0;
      for (int jdir=0; jdir<2; ++jdir){
        grad_u0(idir) += -1.*grad_uj(idir,jdir)*u_cov(jdir)/gamma;
      }
    }
    milne::Vector<double,2> grad_u0_contra = grad_u0;
    grad_u0_contra.make_contravariant(t2);

    milne::Matrix<double,2,2> contra_grad_uj;
    contra_grad_uj = grad_uj;
    contra_grad_uj.make_contravariant(0, t2);
    milne::Matrix<double,3,3> contra_grad_uj_3d;
    for(int idir=0; idir<2; ++idir){
      for(int jdir=0; jdir<2; ++jdir){
        contra_grad_uj_3d(idir,jdir) = contra_grad_uj(idir,jdir);
      }
    }
    //fill remaning dimensions
    for(int idir=0; idir<2; ++idir){
      contra_grad_uj_3d(idir,2) = 0.0;
      contra_grad_uj_3d(2,idir) = 0.0;
    }
    contra_grad_uj_3d(2,2) = 0.0;

    //calculate by hand due to shear dimensions being fixed,
    //and to avoid unnecessary specializations
    milne::Vector<double,2> shv_hyb_i0;
    for(int idir=0; idir<2; ++idir){
      shv_hyb_i0(idir) = shv_hybrid(idir+1,0);
    }
    
    milne::Matrix<double,2,2> vi_shv_hyb_0j;
    for(int idir=0; idir<2; ++idir){
      for(int jdir=0; jdir<2; ++jdir){
        vi_shv_hyb_0j(idir,jdir) =shv_hybrid(idir+1,jdir+1) -  v(idir)*shv_hybrid(0,jdir+1);
      }
    }

    milne::Matrix<double,2,3> R_0i_shear;
    for(int idir=0; idir<2; ++idir){
      for(int icharge=0; icharge<3; ++icharge){
        R_0i_shear(idir,icharge) = 0.0;
      }
    }


    //aux vectors and matrices
    for(int idir=0; idir<2; ++idir){
      F_i0_sigma(idir) = -((v*grad_uj)(idir)
                          -grad_u0_contra(idir))/2.
                          +u(idir)*gamma*(gamma/t + gamma*divV)/3.
                          -gamma*geometric_factor*u(idir)/2.;
      F_i0_D(idir) = 0.0;
      F_i0_domega(idir) = 0.0;
      F_i0_dsigma(idir) = 0.0;
      F_i0_dd(idir) = shv(idir+1,0)*(geometric_factor+gamma/t+gamma*divV);
      for(int jdir=0; jdir<2; ++jdir){
        M_i0_sigma(idir,jdir) = u(idir)*u_cov(jdir)/(2.*3.);
        M_i0_D(idir,jdir) = gamma*u(idir)*(shv_hybrid(0,jdir+1)-shv_hybrid(0,0)*u_cov(jdir)/gamma)
                            +gamma*gamma*(shv_hybrid(idir+1,jdir+1)-shv_hybrid(idir+1,0)*u_cov(jdir)/gamma);
        M_i0_domega(idir,jdir) = 0.0;
        M_i0_dsigma(idir,jdir) = 0.0;
        M_i0_dd(idir,jdir) = -shv(idir+1,0)*u_cov(jdir)/gamma;
      }
      //diagonal terms
      M_i0_sigma(idir,idir) += (1.-gamma*gamma)/2.;

    }


    //fill F vectors for shear viscosity
    milne::Vector<double,4> shear0mu_aux;
    for(int idir=0; idir<4; ++idir){
      shear0mu_aux(idir) = shv_hybrid(0,idir);
    }

    double F_shv_nabla_u = shv_hybrid(0,0)*geometric_factor
            +milne::contract(shv_hyb_i0 - shv_hybrid(0,0)*v,grad_u0)
            +milne::contract(vi_shv_hyb_0j,grad_uj)
            +shv_hybrid(3,3)*gamma/t;
    //stores the results
    device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u) = F_shv_nabla_u;

    for(int idir=0; idir<2; ++idir){
      F_0i_shear(idir) = (-shv(idir+1,0)-tilde_delta*F_i0_dd(idir)-tau_pi*F_i0_D(idir)
                          -a*2.*tau_pi*F_i0_domega(idir) 
                          -a*tau_pipi*F_i0_dsigma(idir)
                          +(2.*eta_pi+a*lambda_piPi*bulk)*F_i0_sigma(idir)
                          +a*phi6*bulk*shv(idir+1,0)
                          +a*phi7*(milne::contract(shv, shear0mu_aux, milne::SecondIndex())(idir+1)
                          +gamma*u(idir)*milne::contract(shv_cov,shv)/3.))/(tau_pi*gamma)
                          -shv(idir+1,0)*(geometric_factor+gamma/t+gamma*divV)/(gamma);
      //stores the results 
      device_hydro_vector.access(is, ia, hydro_info::F_0i_shear, idir) = F_0i_shear(idir);
    }    
    

    //fill M matrices for shear viscosity
    for(int idir=0; idir<2; ++idir){
      M_shv_nabla_u(idir) = shv_hybrid(0,idir+1) - shv_hybrid(0,0)*u_cov(idir)/gamma;
      for(int jdir=0; jdir<2; ++jdir){
        M_0i_shear(idir,jdir) = -(tau_pi*(M_i0_D(idir,jdir)
                                +tilde_delta*M_i0_dd(idir,jdir)/tau_pi
                                +2.*a*M_i0_domega(idir,jdir))
                                +a*tau_pipi*M_i0_dsigma(idir,jdir)
                                -(2.*eta_pi+a*lambda_piPi*bulk)*M_i0_sigma(idir,jdir))/(tau_pi*gamma)
                                +shv(idir+1,0)*u_cov(jdir)/(gamma*gamma);
        //stores the results
        device_hydro_space_matrix.access(is, ia, hydro_info::M_0i_shear, idir, jdir) = M_0i_shear(idir,jdir);
      }
      //stores the results
      device_hydro_vector.access(is, ia, hydro_info::M_shv_nabla_u, idir) = M_shv_nabla_u(idir);
    }


    
    milne::Vector<double,4> contra_metric_diag = milne::get_contra_metric_diagonal<3>(t);

    milne::Matrix3D<double,2,3,2> M_D;
    milne::Matrix3D<double,2,3,2> M_sigma;
    milne::Matrix3D<double,2,3,2> M_delta;
    milne::Matrix3D<double,2,3,2> M_domega;
    milne::Matrix3D<double,2,3,2> M_dsigma;

    milne::Matrix<double,2,3> F_D;
    milne::Matrix<double,2,3> F_sigma;
    milne::Matrix<double,2,3> F_delta;
    milne::Matrix<double,2,3> F_domega;
    milne::Matrix<double,2,3> F_dsigma;

    for(int idir=0; idir<2; ++idir){

      for(int jdir=idir; jdir<3; ++jdir){
        F_D(idir,jdir) = gamma*geometric_factor*(four_u(idir+1)*shv_hybrid(jdir+1,0)
                  +four_u(jdir+1)*shv_hybrid(idir+1,0));
                
        F_sigma(idir,jdir) = contra_grad_uj_3d(idir,jdir)/2. + contra_grad_uj_3d(jdir,idir)/2.
                  +four_u(idir+1)*four_u(jdir+1)*(geometric_factor+gamma/t+gamma*divV)/3.;         
        F_delta(idir,jdir) = shv(idir+1,jdir+1)*(geometric_factor+gamma/t+gamma*divV);
        F_domega(idir,jdir) = 0.0;
        F_dsigma(idir,jdir) = 0.0;                

        for(int kdir=0; kdir<2; ++kdir){
          M_D(idir,jdir,kdir) = gamma*(four_u(idir+1)*shv_hybrid(jdir+1,kdir+1)
                                      -four_u(idir+1)*shv_hybrid(jdir+1,0)*four_u_cov(kdir+1)/gamma
                                      +four_u(jdir+1)*shv_hybrid(idir+1,kdir+1)
                                      -four_u(jdir+1)*shv_hybrid(idir+1,0)*four_u_cov(kdir+1)/gamma);
          M_sigma(idir,jdir,kdir) = -four_u(idir+1)*four_u(jdir+1)*four_u_cov(kdir+1)/(3.*gamma);
          M_delta(idir,jdir,kdir) = -shv(idir+1,jdir+1)*four_u_cov(kdir+1)/gamma;
          M_domega(idir,jdir,kdir) = 0.0;
          M_dsigma(idir,jdir,kdir) = 0.0;
        }
        //diagonal i = k terms
        M_sigma(idir,jdir,idir) += -gamma*four_u(jdir+1)/2.;
      }
      for(int jdir = idir; jdir<2; ++jdir){
        //j=k, skipping j=3 
        M_sigma(idir,jdir,jdir) += -gamma*four_u(idir+1)/2.; 
      }

      //diagonal j terms
      for(int kdir=0; kdir<2; ++kdir){
        //j=i
        M_sigma(idir,idir,kdir) += contra_metric_diag(idir+1)*four_u_cov(kdir+1)/(gamma*3.);
      }
      
      //diagonal and metric terms for F
      F_D(idir,2) += gamma*shv(3,idir+1)/t;
      F_sigma(idir,idir) += -contra_metric_diag(idir+1)*(geometric_factor+gamma/t+gamma*divV)/3.;

    }

    for(int idir=0; idir<2; ++idir){
      for(int jdir=idir; jdir<3; ++jdir){
        F_big_shear(idir,jdir) = (-shv(idir+1,jdir+1) 
                                -tilde_delta*F_delta(idir,jdir)
                                -tau_pi*F_D(idir,jdir)
                                -a*2.*tau_pi*F_domega(idir,jdir)
                                -a*tau_pipi*F_dsigma(idir,jdir)
                                +(2.*eta_pi+a*lambda_piPi*bulk)*F_sigma(idir,jdir)
                                +a*phi6*bulk*shv(idir+1,jdir+1)
                                +a*phi7*(milne::contract(shv,shv_hybrid,milne::SecondIndex(),milne::SecondIndex())(idir,jdir)
                                +four_u(idir+1)*four_u(jdir+1)*milne::contract(shv,shv_cov)/3.))/(sigma*tau_pi*gamma);
        //std::cout << "F_big_shear(" << idir << "," << jdir << ") = " << F_big_shear(idir,jdir) << std::endl;


        for(int kdir=0; kdir<2; ++kdir){
          M_big_shear(idir,jdir,kdir) = -(tau_pi*(M_D(idir,jdir,kdir)
                                      +tilde_delta*M_delta(idir,jdir,kdir)/tau_pi
                                      +2.*a*M_domega(idir,jdir,kdir))
                                      +a*tau_pipi*M_dsigma(idir,jdir,kdir)
                                      -(2.*eta_pi+a*lambda_piPi*bulk)*M_sigma(idir,jdir,kdir))/(sigma*tau_pi*gamma);
          //saves the results 
          int linear_index = jdir * 2 + kdir;
          device_hydro_shear_aux_matrix.access(is, ia, hydro_info::M_big_shear, idir, linear_index) = M_big_shear(idir,jdir,kdir);
        }
      }
      F_big_shear(idir,idir) += a*phi7*(-contra_metric_diag(idir+1)*milne::contract(shv,shv_cov))/(3.*sigma*tau_pi*gamma);
      //saves the results
      for(int jdir=idir; jdir<3; ++jdir){
        device_hydro_shear_aux_vector.access(is, ia, hydro_info::F_big_shear, idir, jdir) = F_big_shear(idir,jdir);
      }

                          
    }   
    //stores R matrix
    for(int idir=0; idir<2; ++idir){
      for(int icharge=0; icharge<3; ++icharge){
        device_hydro_space_matrix.access(is, ia, hydro_info::R_0i_shear, idir, icharge) = R_0i_shear(idir,icharge);
      }
    } 
    for(int idir=0; idir<2; ++idir){
      for(int jdir=idir; jdir<3; ++jdir){
        for(int icharge=0; icharge<3; ++icharge){
          int linear_index = jdir * 3 + icharge;
          device_hydro_shear_aux_matrix.access(is, ia, hydro_info::R_big_shear, idir, linear_index) =0.;
        }
      }
    }
    
  };
  Cabana::simd_parallel_for(simd_policy, compute_MRF_shear, "compute_MRF_shear");
  Kokkos::fence();
};
    
    
    
    
    
    
template<>
void EoM_default<1>::calculate_MRF_shear(std::shared_ptr<SystemState<1>> sysPtr)
{
  double t = (sysPtr->t);
  double t2 = t*t;
  milne::Vector<double,1> delta_i_eta = milne::delta_i_eta<1>();
  CREATE_VIEW(device_,sysPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());

  auto compute_MRF_shear = KOKKOS_LAMBDA(const int is, const int ia)
  {
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double eta_pi = device_hydro_scalar.access(is, ia, hydro_info::eta_pi);
    double phi6 = device_hydro_scalar.access(is, ia, hydro_info::phi6);
    double phi7 = device_hydro_scalar.access(is, ia, hydro_info::phi7);
    double lambda_piPi = device_hydro_scalar.access(is, ia, hydro_info::lambda_piPi);
    double tau_pipi = device_hydro_scalar.access(is, ia, hydro_info::tau_pipi);
    double a= device_hydro_scalar.access(is, ia, hydro_info::a);
    double delta_pipi = device_hydro_scalar.access(is, ia, hydro_info::delta_pipi);
    double tau_pi = device_hydro_scalar.access(is, ia, hydro_info::tau_pi);
    double tilde_delta = delta_pipi - tau_pi;
    double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);

    milne::Vector<double,1> u, u_cov, M_shv_nabla_u ,F_0i_shear;
    milne::Vector<double,1> v, grad_u0;
    milne::Vector<double,4> four_u,four_u_cov;
    milne::Matrix<double,2,3> F_big_shear;
    milne::Matrix<double,1,1> gradV, grad_uj, M_0i_shear;
    milne::Matrix<double,4,4> shv_hybrid, shv, shv_cov;
    milne::Matrix3D<double,2,3,1> M_big_shear;
    for(int idir=0; idir<1; ++idir){
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
    }
    four_u(0) = gamma;
    four_u(1) = 0.0;
    four_u(2) = 0.0;
    four_u(3) = u(0);

    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
    u_cov.make_covariant(t2);
    four_u_cov(0) = gamma;
    four_u_cov(1) = 0.0;
    four_u_cov(2) = 0.0;
    four_u_cov(3) = u_cov(0);
    double geometric_factor = u(0)*u(0)*t/gamma;
    shv_hybrid = shv;
    shv_hybrid.make_covariant(1, t2);
    shv_cov = shv_hybrid;
    shv_cov.make_covariant(0, t2);

    milne::Vector<double,1> F_i0_sigma;
    milne::Vector<double,1> F_i0_D;
    milne::Vector<double,1> F_i0_domega;
    milne::Vector<double,1> F_i0_dsigma;
    milne::Vector<double,1> F_i0_dd;
    milne::Matrix<double,1,1> M_i0_sigma;
    milne::Matrix<double,1,1> M_i0_D;
    milne::Matrix<double,1,1> M_i0_domega;
    milne::Matrix<double,1,1> M_i0_dsigma;
    milne::Matrix<double,1,1> M_i0_dd;



    double divV = 0;
    for(int idir=0; idir<1; ++idir){
      divV += device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, idir);
    for(int jdir=0; jdir<1; ++jdir){
      gradV(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, jdir);
      }
    }

    grad_uj = gamma*gradV+pow(gamma,3.)*(v*(v*gradV));
    grad_u0 = -1.*milne::contract(grad_uj,u_cov,milne::SecondIndex())/gamma;
    milne::Vector<double,1> grad_u0_contra = grad_u0;
    grad_u0_contra.make_contravariant(t2);

    milne::Matrix<double,1,1> contra_grad_uj;
    contra_grad_uj = grad_uj;
    contra_grad_uj.make_contravariant(0, t2);
    milne::Matrix<double,3,3> contra_grad_uj_3d;
    contra_grad_uj_3d(2,2) = contra_grad_uj(0,0);
    //fill remaning dimensions
    for(int idir=0; idir<2; ++idir){
      for(int jdir=0; jdir<2; ++jdir){
        contra_grad_uj_3d(idir,jdir) = 0.0;
      }
      contra_grad_uj_3d(idir,2) = 0.0;
      contra_grad_uj_3d(2,idir) = 0.0;
    }
    //calculate by hand due to shear dimensions being fixed,
    //and to avoid unnecessary specializations
    //calculate by hand due to shear dimensions being fixed,
    //and to avoid unnecessary specializations
    milne::Vector<double,1> shv_hyb_i0;
    for(int idir=0; idir<1; ++idir){
      shv_hyb_i0(idir) = shv_hybrid(idir+1,0);
    }
    
    milne::Matrix<double,1,1> vi_shv_hyb_0j;
    for(int idir=0; idir<1; ++idir){
      for(int jdir=0; jdir<1; ++jdir){
        vi_shv_hyb_0j(idir,jdir) =shv_hybrid(idir+1,jdir+1) -  v(idir)*shv_hybrid(0,jdir+1);
      }
    }

    milne::Matrix<double,1,3> R_0i_shear;
    for(int idir=0; idir<1; ++idir){
      for(int icharge=0; icharge<3; ++icharge){
        R_0i_shear(idir,icharge) = 0.0;
      }
    }


    //aux vectors and matrices
    for(int idir=0; idir<1; ++idir){
      F_i0_sigma(idir) = -(milne::contract(grad_uj,v,milne::FirstIndex())(idir)
                          -grad_u0_contra(idir))/2.
                          +u(idir)*gamma*(geometric_factor
                          +gamma/t + gamma*divV)/3.
                          -gamma*geometric_factor*u(idir)/2.;     
      F_i0_D(idir) = gamma*(u(idir)*shv_hybrid(0,0)+gamma*shv_hybrid(idir+1,0))*geometric_factor
                     +t*u(2)*shv(3,idir+1);
      F_i0_domega(idir) = 0.0;
      F_i0_dsigma(idir) = 0.0;
      F_i0_dd(idir) = shv(idir+1,0)*(geometric_factor+gamma/t+gamma*divV);
      for(int jdir=0; jdir<1; ++jdir){
        M_i0_sigma(idir,jdir) = u(idir)*u_cov(jdir)/3.;
        M_i0_D(idir,jdir) = gamma*u(idir)*(shv_hybrid(0,jdir+1)-shv_hybrid(0,0)*u_cov(jdir)/gamma)
                            +gamma*gamma*(shv_hybrid(idir+1,jdir+1)-shv_hybrid(idir+1,0)*u_cov(jdir)/gamma);
        M_i0_domega(idir,jdir) = 0.0;
        M_i0_dsigma(idir,jdir) = 0.0;
        M_i0_dd(idir,jdir) = -shv(idir+1,0)*u_cov(jdir)/3.;
      }
      //diagonal terms
      M_i0_sigma(idir,idir) += (1.-gamma*gamma);

    }
    F_i0_sigma(0) += -(-u(0)*t-u(0)/t+gamma*u(0)*u(0)*t
                    +2.*gamma*gamma*u(0)/t)/2.;
    F_i0_D(0) += u(0)*shv(0,0)/t + gamma*shv(3,0)/t;

    //fill F vectors for shear viscosity
    milne::Vector<double,4> shear0mu_aux;
    for(int idir=0; idir<4; ++idir){
      shear0mu_aux(idir) = shv_hybrid(0,idir);
    }

    double F_shv_nabla_u = shv_hybrid(0,0)*geometric_factor
            +milne::contract(shv_hyb_i0 - shv_hybrid(0,0)*v,grad_u0)
            +milne::contract(vi_shv_hyb_0j,grad_uj)
            +shv_hybrid(3,3)*gamma/t
            + (t*shv_hybrid(3,0)+shv_hybrid(0,3)/t)*u(0);
    //stores the results
    device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u) = F_shv_nabla_u;

    for(int idir=0; idir<1; ++idir){
      F_0i_shear(idir) = (-tilde_delta*F_i0_dd(idir)-tau_pi*F_i0_D(idir)
                          -a*2.*tau_pi*F_i0_domega(idir) 
                          -a*tau_pipi*F_i0_dsigma(idir)
                          +(2.*eta_pi+a*lambda_piPi*bulk)*F_i0_sigma(idir)
                          +a*phi6*bulk*shv(idir+1,0)
                          +a*phi7*(milne::contract(shv, shear0mu_aux, milne::SecondIndex())(idir+1)
                          +gamma*u(idir)*milne::contract(shv_cov,shv)/3.)/(tau_pi*gamma)
                          -shv(idir+1,0)*(geometric_factor+gamma/t+gamma*divV)/(gamma));
      //stores the results 
      device_hydro_vector.access(is, ia, hydro_info::F_0i_shear, idir) = F_0i_shear(idir);
    }    
    

    //fill M matrices for shear viscosity
    for(int idir=0; idir<1; ++idir){
      M_shv_nabla_u(idir) = shv_hybrid(0,idir+1) - shv_hybrid(0,0)*u_cov(idir)/gamma;
      for(int jdir=0; jdir<1; ++jdir){
        M_0i_shear(idir,jdir) = -(tau_pi*(M_i0_D(idir,jdir)
                                +tilde_delta*M_i0_dd(idir,jdir)/tau_pi
                                +2.*a*M_i0_domega(idir,jdir))
                                +a*tau_pipi*M_i0_dsigma(idir,jdir)
                                -(2.*eta_pi+a*lambda_piPi*bulk)*M_i0_sigma(idir,jdir))/(tau_pi*gamma)
                                +shv(idir+1,0)*u_cov(jdir)/(gamma*gamma);
        //stores the results
        device_hydro_space_matrix.access(is, ia, hydro_info::M_0i_shear, idir, jdir) = M_0i_shear(idir,jdir);
      }
      //stores the results
      device_hydro_vector.access(is, ia, hydro_info::M_shv_nabla_u, idir) = M_shv_nabla_u(idir);
    }


    
    milne::Vector<double,4> contra_metric_diag = milne::get_contra_metric_diagonal<3>(t);

    milne::Matrix3D<double,2,3,1> M_D;
    milne::Matrix3D<double,2,3,1> M_sigma;
    milne::Matrix3D<double,2,3,1> M_delta;
    milne::Matrix3D<double,2,3,1> M_domega;
    milne::Matrix3D<double,2,3,1> M_dsigma;

    milne::Matrix<double,2,3> F_D;
    milne::Matrix<double,2,3> F_sigma;
    milne::Matrix<double,2,3> F_delta;
    milne::Matrix<double,2,3> F_domega;
    milne::Matrix<double,2,3> F_dsigma;

    for(int idir=0; idir<2; ++idir){

      for(int jdir=idir; jdir<3; ++jdir){
        F_D(idir,jdir) = gamma*geometric_factor*(four_u(idir+1)*shv_hybrid(jdir+1,0)
                  +four_u(jdir+1)*shv_hybrid(idir+1,0));
                
        F_sigma(idir,jdir) = contra_grad_uj(idir,jdir)/2. + contra_grad_uj(jdir,idir)/2.
                  -four_u(idir+1)*four_u(jdir+1)*(geometric_factor+gamma/t+gamma*divV)/3.;         
        F_delta(idir,jdir) = shv(idir+1,jdir+1)*(geometric_factor+gamma/t+gamma*divV);
        F_domega(idir,jdir) = 0.0;
        F_dsigma(idir,jdir) = 0.0;                

        for(int kdir=0; kdir<1; ++kdir){
          M_D(idir,jdir,kdir) = gamma*(four_u(idir+1)*shv_hybrid(jdir+1,kdir+1)
                                      -four_u(idir+1)*shv_hybrid(jdir+1,0)*four_u_cov(kdir+1)/gamma
                                      +four_u(jdir+1)*shv_hybrid(idir+1,kdir+1)
                                      -four_u(jdir+1)*shv_hybrid(idir+1,0)*four_u_cov(kdir+1)/gamma);
          M_sigma(idir,jdir,kdir) = -2.*four_u(idir+1)*four_u(jdir+1)*four_u_cov(kdir+1)/(3.*gamma);
          M_delta(idir,jdir,kdir) = -shv(idir+1,jdir+1)*four_u_cov(kdir+1)/gamma;
          M_domega(idir,jdir,kdir) = 0.0;
          M_dsigma(idir,jdir,kdir) = 0.0;
        }
        //diagonal k terms
        M_sigma(idir,jdir,jdir) += -gamma*four_u(idir+1); 
        M_sigma(idir,jdir,idir) += -gamma*four_u(jdir+1);
      }
      //diagonal j terms
      for(int kdir=0; kdir<1; ++kdir){
        M_sigma(idir,idir,kdir) += 2.*contra_metric_diag(idir+1)*four_u_cov(kdir+1)/3.;
      }

      //diagonal and metric terms for F
      F_D(idir,2) +=u(0)*shv(0,idir+1)/t + gamma*shv(3,idir+1)/t;
      F_sigma(idir,idir) += -contra_metric_diag(idir+1)*(geometric_factor+gamma/t+gamma*divV)/3.;
      F_sigma(idir,2) += four_u(idir+1)*u(0)*u(0)*t + 2.*four_u(idir+1)*u(0)*gamma/t;

    }

    for(int idir=0; idir<2; ++idir){
      for(int jdir=idir; jdir<3; ++jdir){
        F_big_shear(idir,jdir) = (-shv(idir+1,jdir+1) 
                                -tilde_delta*F_delta(idir,jdir)
                                -tau_pi*F_D(idir,jdir)
                                -a*2.*tau_pi*F_domega(idir,jdir)
                                -a*tau_pipi*F_dsigma(idir,jdir)
                                +(2.*eta_pi+a*lambda_piPi*bulk)*F_sigma(idir,jdir)
                                +a*phi6*bulk*shv(idir+1,jdir+1)
                                +a*phi7*(milne::contract(shv,shv_hybrid,milne::SecondIndex(),milne::SecondIndex())(idir,jdir)
                                +four_u(idir)*four_u(jdir)*milne::contract(shv,shv_cov)/3.))/(sigma*tau_pi*gamma);

        for(int kdir=0; kdir<1; ++kdir){
          
          M_big_shear(idir,jdir,kdir) = -(tau_pi*(M_D(idir,jdir,kdir)
                                      +tilde_delta*M_delta(idir,jdir,kdir)/tau_pi
                                      +2.*a*M_domega(idir,jdir,kdir))
                                      +a*tau_pipi*M_dsigma(idir,jdir,kdir)
                                      -(2.*eta_pi+a*lambda_piPi*bulk)*M_sigma(idir,jdir,kdir))/(sigma*tau_pi*gamma);
          //saves the results 
          int linear_index = jdir * 1 + kdir;
          device_hydro_shear_aux_matrix.access(is, ia, hydro_info::M_big_shear, idir, linear_index) = M_big_shear(idir,jdir,kdir);
        }
      }
      F_big_shear(idir,idir) += a*phi7*(-contra_metric_diag(idir+1)*milne::contract(shv,shv_cov))/(3.*sigma*tau_pi*gamma);
      //saves the results
      for(int jdir=idir; jdir<3; ++jdir){
        device_hydro_shear_aux_vector.access(is, ia, hydro_info::F_big_shear, idir, jdir) = F_big_shear(idir,jdir);
      }

                          
    }   
    //stores R matrix
    for(int idir=0; idir<1; ++idir){
      for(int icharge=0; icharge<3; ++icharge){
        device_hydro_space_matrix.access(is, ia, hydro_info::R_0i_shear, idir, icharge) = R_0i_shear(idir,icharge);
      }
    }
    for(int idir=0; idir<2; ++idir){
      for(int jdir=idir; jdir<3; ++jdir){
        for(int icharge=0; icharge<3; ++icharge){
          int linear_index = jdir * 3 + icharge;
          device_hydro_shear_aux_matrix.access(is, ia, hydro_info::R_big_shear, idir, linear_index) = 0.;
        }
      }
    }
    
  };
  Cabana::simd_parallel_for(simd_policy, compute_MRF_shear, "compute_MRF_shear");
  Kokkos::fence();
};
   
   
   
   
   
   
   
   