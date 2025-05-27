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
    double gamma = Kokkos::sqrt(1-milne::contract(u_cov,u)); //Calculates gamma = \sqrt{1-u^i u_i}
    device_hydro_scalar.access(is,ia, hydro_info::gamma) = gamma;
    for(int idir=0; idir<D; ++idir) device_hydro_vector.access(is, ia, hydro_info::v, idir) = u(idir)/gamma;
  };

  Cabana::simd_parallel_for(simd_policy, update_gammas, "update_gamma_sigma_kernel");
  Kokkos::fence();
}


/// @brief Enforces the constraints for the shear viscous tensor \f$ \pi^{\mu\nu} \f$.
/// @details Calculate the shear viscous tensor \f$ \pi^{\mu\nu} \f$ and the bulk viscous pressure
/// \f$ \Pi \f$ from the extensive (called extensive) shear tensor \f$ \pi^{ij} \f$
//  and the extensive(extensive) bulk pressure \f$ \Pi \f$
/// It also enforces the constraints  \f$\pi^{\mu\nu}u_\nu = 0 \f$,
/// \f$\pi^\mu_\mu = 0 \f$ and \f$ \pi^{\mu\nu} = \pi^{\nu\mu} \f$.
/// We assume that the components \f$ \pi^{xx} \f$, \f$ \pi^{xy} \f$,
/// \f$ \pi^{x\eta} \f$, \f$ \pi^{yy} \f$ and \f$ \pi^{y\eta} \f$ are computed
/// during evolution
/// \f[\begin{align*}
/// \pi^{0i} & = -\pi^{ij}u_j/\gamma \\
/// \pi^{00} & = u_i u_j \pi^{ij}/\gamma^2 \\
/// \end{align*}\f]
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
    //read relevant quantities
    double gamma = device_hydro_scalar.access(is,ia, hydro_info::gamma);
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    double extensive_bulk = device_hydro_scalar.access(is, ia, hydro_info::extensive_bulk);
    //Declare caches
    milne::Matrix<double,4,4> shv;
    milne::Vector<double,D> u;
    milne::Vector<double,3> fixed_size_u, fixed_size_u_cov;
    //fill caches
    for(int idir=0; idir<D; ++idir){
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      fixed_size_u(idir) = u(idir);
    }
    for(int idir=D; idir<3; ++idir){
      fixed_size_u(idir) = 0.0;
    }
    milne::Vector<double,D> u_cov = u;
    u_cov.make_covariant(t2); //Transforms u^i to u_i
    fixed_size_u_cov = fixed_size_u;
    fixed_size_u_cov.make_covariant(t2); //Transforms fu^i to fu_i


    //compute bulk from extensive bulk
    device_hydro_scalar.access(is, ia, hydro_info::bulk) = extensive_bulk*sigma;

    //compute shear from extensive shear
    for(int idir=0; idir<2; ++idir)
    for(int jdir=idir; jdir<3; ++jdir)
      shv(idir+1,jdir+1) = device_hydro_shear_aux_vector.access(is, ia, hydro_info::extensive_shv, idir, jdir)*sigma;



    //Symmetrizes space part
    for( int i=1; i<4; i++ )
    for( int j=i+1; j<4; j++ )
      shv(j,i) = shv(i,j);

    //pi^{0i} = -\pi^{ij}u_j/gamma
    for(int idir=1; idir<D+1; ++idir){
      milne::Vector<double,D> shv_i;
      for(int jdir=1; jdir<D+1; ++jdir) shv_i(jdir-1) = shv(idir,jdir);
      shv(0,idir) = -1.*milne::contract(u_cov,shv_i)/gamma;
    }

    //Symmetrizes time part
    for( int i=1; i<4; i++ )
      shv(i,0) = shv(0,i);


    //pi^00 = u_i u_j pi^{ij}/gamma^2
    shv(0,0) =( fixed_size_u_cov(0)*fixed_size_u_cov(0)*shv(1,1)
                + fixed_size_u_cov(1)*fixed_size_u_cov(1)*shv(2,2)
                + 2.*fixed_size_u_cov(1)*fixed_size_u_cov(0)*shv(2,1)
                + 2.*fixed_size_u_cov(2)*fixed_size_u_cov(0)*shv(3,1)
                + 2.*fixed_size_u_cov(2)*fixed_size_u_cov(1)*shv(3,2)
                - fixed_size_u_cov(2)*fixed_size_u_cov(2)*(shv(1,1)+shv(2,2))/t2
              )/(gamma*gamma-fixed_size_u_cov(2)*fixed_size_u_cov(2)/(t2));


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

/// @brief Calculates the Lorentz contraction factor \f$ \gamma = u^0\f$.
/// @details Calculates the Lorentz contraction factor \f$ \gamma \f$ in the
/// general case.
/// @param u The velocity vector (space components only).
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return The value of \f$ \gamma \f$.
template<unsigned int D> KOKKOS_FUNCTION
double EoM_default<D>::gamma_calc(double u_vec[D], const double &time_squared)
{
    milne::Vector<double,D> u,u_cov;
    for (unsigned int i=0; i<D; i++) u(i) = u_vec[i];
    u_cov = u;
    u_cov.make_covariant(time_squared);
    double gamma = Kokkos::sqrt(1-milne::contract(u_cov,u));
    return gamma;
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

/// @brief Calculates the time derivatives of the hydrodynamic variables.
/// @details Calculates the time derivatives of the hydrodynamic variables
/// using the MRF formalism.
/// @tparam D The number of spatial dimensions.
/// @param sysPtr A pointer to an object of class SystemState
/// @param settingsPtr A pointer to an object of class Settings
template<unsigned int D>
void EoM_default<D>::evaluate_time_derivatives( std::shared_ptr<SystemState<D>> sysPtr, std::shared_ptr<Settings> settingsPtr)
{
  // #ifdef DEBUG
  // ofstream outfile;
  // outfile.open("gradients.dat");
  // #endif

  std::cout << "===========================================================================\n";
  std::cout << "Particle #0 at " << __FUNCTION__ << "::" << __LINE__ << ":\n";
  std::cout << sysPtr->particles[0] << std::endl;
  sysPtr->print_neighbors(0);

  double t = (sysPtr->t);
  double t2 = t*t;
  milne::Vector<double,D> delta_i_eta = milne::delta_i_eta<D>();
  CREATE_VIEW(device_,sysPtr->cabana_particles);
  bool using_shear = settingsPtr->using_shear;
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());
  //calculate the M,R,F matrices due to shear, when using shear
  if(using_shear){
    calculate_MRF_shear(sysPtr);
  }
   auto fill_auxiliary_variables = KOKKOS_LAMBDA(int const is, int const ia){
    //read relevant quantities
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
    double sigma_lab = device_hydro_scalar.access(is, ia, hydro_info::sigma_lab);
    double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);
    double zeta = device_hydro_scalar.access(is, ia, hydro_info::zeta_Pi);
    double tau_Pi = device_hydro_scalar.access(is, ia, hydro_info::tau_Pi);
    double tau_pi = device_hydro_scalar.access(is, ia, hydro_info::tau_pi );
    double delta_PiPi = device_hydro_scalar.access(is, ia, hydro_info::delta_PiPi);
    double a = device_hydro_scalar.access(is, ia, hydro_info::a);
    double F_extensive_bulk = device_hydro_scalar.access(is, ia, hydro_info::F_extensive_bulk);
    double F_shv_nabla_u = device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u);
    double F_extensive_entropy = device_hydro_scalar.access(is, ia, hydro_info::F_extensive_entropy);
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    double phi1 = device_hydro_scalar.access(is, ia, hydro_info::phi1);
    double phi3 = device_hydro_scalar.access(is, ia, hydro_info::phi3);
    double lambda_Pipi = device_hydro_scalar.access(is, ia, hydro_info::lambda_Pipi);
    double j0_ext = device_hydro_scalar.access(is, ia, hydro_info::j0_ext);
    double rhoQ_ext = device_hydro_scalar.access(is, ia, hydro_info::rho_Q_ext);
    double rhoS_ext = device_hydro_scalar.access(is, ia, hydro_info::rho_S_ext);
    double rhoB_ext = device_hydro_scalar.access(is, ia, hydro_info::rho_B_ext);
    //auxiliary zeta tilde to control IR or DNMR
    double zeta_tilde  = zeta +  a*(delta_PiPi - tau_Pi)*bulk;
    //declare caches
    milne::Vector<double,D> v, u, grad_u0, u_cov, j_ext;
    milne::Vector<double,D> M_extensive_bulk_aux, M_shv_nabla_u, M_extensive_entropy;
    milne::Vector<double,3> R_extensive_entropy, R_extensive_bulk, F_extensive_N;
    milne::Matrix<double,D,D> gradV, grad_uj;
    milne::Matrix<double,3,3> R_extensive_N;
    milne::Matrix<double,4,4> shv, shv_hybrid, shv_cov;
    //fill caches
    for(int idir=0; idir<D; ++idir){
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
      j_ext(idir) = device_hydro_vector.access(is, ia, hydro_info::j_ext, idir);
    }
    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
    u_cov = u;
    u_cov.make_covariant(t2);
    shv_hybrid = shv;
    //shv_hybrid = \pi^\mu_\nu
    shv_hybrid.make_covariant(1, t2);
    shv_cov = shv_hybrid;
    //shv_cov = \pi_{\mu\nu}
    shv_cov.make_covariant(0, t2);

    double divV = 0;
    for(int idir=0; idir<D; ++idir){
      divV += device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, idir);
    for(int jdir=0; jdir<D; ++jdir){
      gradV(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, jdir);
      }
    }
    double geometric_factor = delta_i_eta(D-1)*u(D-1)*u(D-1)*t/gamma;

    //auxiliary array to store muB, muS, muQ
    milne::Vector<double,3> mu_vec = {muB, muS, muQ};
    //auxiliary array to store rho_ext
    milne::Vector<double,3> rho_ext = {rhoB_ext, rhoS_ext, rhoQ_ext};


    //double j0_ext =0.;
    //milne::Vector<double,D> j_ext;
    //for(int idir=0; idir<D; ++idir){
    //  j_ext(idir) = 0.;
    //}

    //fill M matrices
    for(int idir=0; idir<D; ++idir){
      M_extensive_bulk_aux(idir) = (zeta_tilde*u_cov(idir)/gamma)/(sigma*gamma*tau_Pi);
      M_extensive_entropy(idir) = (bulk*u_cov(idir)/gamma)/(sigma*gamma*T);
      //if(M_extensive_entropy(idir) > 1e-10){
      //  std::cout << "M_extensive_entropy: " << M_extensive_entropy(idir) << std::endl;
      //}
    };
    //fill R matrices
    for(int icharge=0; icharge<3; ++icharge){
      R_extensive_entropy(icharge) = -gamma*sigma*mu_vec(icharge);

      R_extensive_bulk(icharge) = 0.0;
      for(int jcharge=0; jcharge<3; ++jcharge){
        R_extensive_N(icharge,jcharge) = 0.0;
      }
      R_extensive_N(icharge,icharge) += gamma*sigma;
    }



    //fill F matrices
    F_extensive_bulk = -(zeta_tilde*(gamma*divV + gamma/t + geometric_factor)
                  +bulk  - a*phi1*bulk*bulk);
    F_extensive_bulk =0.0;
    F_extensive_entropy = ( -bulk*(gamma*divV + gamma/t + geometric_factor)
                    + gamma*j0_ext
                    +milne::contract(u_cov,j_ext))/(sigma*gamma*T);


    for(int icharge=0; icharge<3; ++icharge){
      F_extensive_N(icharge) = rho_ext(icharge);
    }
    if(using_shear){
      //add contribution from shear MRFs
      double F_shv_nabla_u = device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u);
      F_extensive_entropy += F_shv_nabla_u/(sigma*gamma*T);
      F_extensive_bulk += -(- a*lambda_Pipi*F_shv_nabla_u
                      - a*phi3*milne::contract(shv_cov,shv))/(sigma*gamma*tau_Pi);
      milne::Vector<double,D> M_shv_nabla_u;
      for(int idir=0; idir<D; ++idir){
        M_shv_nabla_u(idir) =  device_hydro_vector.access(is, ia, hydro_info::M_shv_nabla_u,idir);
        M_extensive_bulk_aux(idir) += a*lambda_Pipi*M_shv_nabla_u(idir)/(sigma*gamma*tau_Pi);
        M_extensive_entropy(idir) += M_shv_nabla_u(idir)/(sigma*gamma*T);
      }
    }

    //stores the results
    device_hydro_scalar.access(is, ia, hydro_info::F_extensive_bulk) = F_extensive_bulk;
    device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u) = F_shv_nabla_u;
    device_hydro_scalar.access(is, ia, hydro_info::F_extensive_entropy) = F_extensive_entropy;
    for(int idir=0; idir<D; ++idir){
      device_hydro_vector.access(is, ia, hydro_info::M_extensive_bulk, idir) = M_extensive_bulk_aux(idir);
      device_hydro_vector.access(is, ia, hydro_info::M_extensive_entropy, idir) = M_extensive_entropy(idir);
    }
    for(int icharge=0; icharge<3; ++icharge){
      device_hydro_vector.access(is, ia, hydro_info::R_extensive_entropy, icharge) = R_extensive_entropy(icharge);
      device_hydro_vector.access(is, ia, hydro_info::F_extensive_N, icharge) = F_extensive_N(icharge);
      device_hydro_vector.access(is, ia, hydro_info::R_extensive_bulk,icharge) = R_extensive_bulk(icharge);
      for(int jcharge=0; jcharge<3; ++jcharge){
        device_hydro_space_matrix.access(is, ia, hydro_info::R_extensive_N, icharge, jcharge) = R_extensive_N(icharge,jcharge);
      }
    }

  };

  // std::cout << "===========================================================================\n";
  // std::cout << "Particle #0 at " << __FUNCTION__ << "::" << __LINE__ << ":\n";
  // std::cout << sysPtr->particles[0] << std::endl;
  // sysPtr->print_neighbors(0);


  Cabana::simd_parallel_for(simd_policy, fill_auxiliary_variables, "fill_auxiliary_variables");
  Kokkos::fence();
  //calculate du/dt
  auto compute_velocity_derivative = KOKKOS_LAMBDA(const int is, const int ia)
  {
      //read relevant quantities
      double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
      double extensive_bulk = device_hydro_scalar.access(is, ia,hydro_info::extensive_bulk );
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
      double F_extensive_bulk = device_hydro_scalar.access(is, ia, hydro_info::F_extensive_bulk);
      double F_extensive_entropy = device_hydro_scalar.access(is, ia, hydro_info::F_extensive_entropy);
      double tau_Pi = device_hydro_scalar.access(is, ia, hydro_info::tau_Pi);
      double j0_ext = device_hydro_scalar.access(is, ia, hydro_info::j_ext);

      //declare caches
      double divV = 0;
      milne::Vector<double,D> u, gradshear, gradP, gradBulk, divshear;
      milne::Vector<double,D> M_bulk, M_S, M_extensive_bulk;
      milne::Vector<double,D> F_0i_shear, j_ext;
      milne:: Vector<double,D> gradP_contra, gradBulk_contra, u_cov;
      milne::Matrix<double,D,D>  M_0i_shear, gradV;
      milne::Vector<double,3> R_S,R_extensive_bulk, rhoVec, dwdrhoVec, F_extensive_N;
      milne::Matrix<double,D,3> R_0i_shear,M_extensive_N;
      milne::Matrix<double,4,4> shv;
      milne::Matrix<double,3,3> R_extensive_N;
      //fill caches
      for(int idir=0; idir<D; ++idir){
        u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
        gradP(idir) = device_hydro_vector.access(is, ia, hydro_info::gradP, idir);
        gradBulk(idir) = device_hydro_vector.access(is, ia, hydro_info::gradBulk, idir);
        divshear(idir) = device_hydro_vector.access(is, ia, hydro_info::divshear, idir);
        gradshear(idir) = device_hydro_vector.access(is, ia, hydro_info::gradshear, idir);
        M_extensive_bulk(idir) = device_hydro_vector.access(is, ia, hydro_info::M_extensive_bulk, idir);
        M_S(idir) = device_hydro_vector.access(is, ia, hydro_info::M_extensive_entropy, idir);
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
        R_S(icharge) = device_hydro_vector.access(is, ia, hydro_info::R_extensive_entropy, icharge);
        R_extensive_bulk(icharge) = device_hydro_vector.access(is, ia, hydro_info::R_extensive_bulk, icharge);
        F_extensive_N(icharge) = device_hydro_vector.access(is, ia, hydro_info::F_extensive_N, icharge);
        for(int idir=0; idir<D; ++idir){
          R_0i_shear(idir,icharge) = device_hydro_space_matrix.access(is, ia, hydro_info::R_0i_shear, idir, icharge);
          M_extensive_N(idir,icharge) = device_hydro_space_matrix.access(is, ia, hydro_info::M_extensive_N, idir, icharge);
        }
        for(int jcharge=0; jcharge<3; ++jcharge){
         R_extensive_N(icharge,jcharge) = device_hydro_space_matrix.access(is, ia, hydro_info::R_extensive_N, icharge, jcharge);
       }
      }
      u_cov = u;
      gradP_contra = gradP;
      gradBulk_contra = gradBulk;
      gradP_contra.make_contravariant(t2);
      gradBulk_contra.make_contravariant(t2);
      u_cov.make_covariant(t2);

      //declare auxiliary variables
      double geometric_factor = delta_i_eta(D-1)*u(D-1)*u(D-1)*t/gamma;
      rhoVec = {rhob, rhos, rhoq};
      dwdrhoVec = {dwdrhoB, dwdrhoS, dwdrhoQ};

      //calculate the aux MRFs from the entalphy derivative
      milne::Vector<double,D> M_w = (s*dwds/gamma
                                    +milne::contract(rhoVec,dwdrhoVec)/gamma)*u_cov/gamma
                                    +sigma*dwds*M_S;
      milne::Vector<double,3> R_w = sigma*dwdrhoVec
                                    +sigma*dwds*R_S;
      double F_w = -(s*dwds/gamma
                    +milne::contract(rhoVec,dwdrhoVec)/gamma)*(geometric_factor+ gamma/t + gamma*divV)
                    +sigma*dwds*F_extensive_entropy;

      //convert the extensive MRFs to the intensive ones
      double F_bulk = F_extensive_bulk - bulk*(gamma*divV + gamma/t + geometric_factor)/gamma;
      M_bulk = M_extensive_bulk + bulk*u_cov/(gamma*gamma*tau_Pi);



      // set the MFs for the four-acceleration
      milne::Matrix <double,D,D> M_u;
      for(int idir=0; idir<D; ++idir){
        for(int jdir=0; jdir<D; ++jdir){
          M_u(idir,jdir) = u(idir)*gamma*(M_extensive_bulk(jdir) + M_w(jdir))
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
          R_u(idir,icharge) = -u(idir)*gamma*(R_w(icharge) + R_extensive_bulk(icharge));
        }
      }
      //add the contribution from the shear MRFs
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
          }
        }
      }

    //aux MRFs
     milne::Matrix<double,D,D> aux_MR;
     milne::Vector<double,D> aux_FR;
     milne::Matrix<double,3,3> RI = milne::inverse(R_extensive_N);

     for(int idir=0; idir<D; ++idir){
       aux_FR(idir) = 0.0;
       for(int jdir=0; jdir<D; ++jdir){
          aux_MR(idir,jdir) = 0.0;
          for(int icharge=0; icharge<3; ++icharge){
            for(int jcharge=0; jcharge<3; ++jcharge){
              aux_MR(idir,jdir) += R_u(idir,icharge)*RI(icharge,jcharge)*M_extensive_N(jdir,jcharge);
            }
          }
          for(int icharge=0; icharge<3; ++icharge){
            for(int jcharge=0; jcharge<3; ++jcharge){
              aux_FR(idir) += R_u(idir,icharge)*RI(icharge,jcharge)*F_extensive_N(jcharge);
            }
          }
        }
      };

      /// @brief
      //milne::Matrix<double,D,D> MI = milne::inverse(M_u-aux_MR);
      //milne::Vector<double,D> du_dt = MI*(F_u+aux_FR);
      milne::Matrix<double,D,D> MI = milne::inverse(M_u-aux_MR);
      milne::Vector<double,D> du_dt = MI*(F_u+aux_FR);
if (std::isnan(du_dt(0)))
{
  std::cout << "Matrices:\n";
  for(int idir=0; idir<D; ++idir)
  for(int jdir=0; jdir<D; ++jdir)
    std::cout <<"\t" << idir << " " << jdir << " " << MI(idir,jdir) << " "
              << M_u(idir,jdir) << " " << M_0i_shear(idir,jdir)
              << u(idir)*gamma*(M_extensive_bulk(jdir) << " " << M_w(jdir))
              << " " << -u(idir)*u_cov(jdir)*(w+bulk)/gamma;
  }
  std::cout <<"\t Scalar: " << gamma*(w+bulk) << "\n";

milne::Vector<double,D> F_u_LOCAL = j_ext + (gradP_contra + gradBulk_contra)
                            -gamma*(F_w+F_bulk)*u
                            -(w+bulk)*(gamma*divV + gamma/t + geometric_factor)*u
                            -2.*(w+bulk)*gamma*u(D-1)*delta_i_eta/t;  std::cout << "Vectors:\n";
  for(int idir=0; idir<D; ++idir)
    std::cout << "\t" << idir << " " << F_u(idir) << " " << -shv(idir+1,0)/t << " "
              << - 2.*delta_i_eta(idir)*shv(0,3)/t
              << " " << -divshear(idir) << " " << gradshear(idir) << " " << -F_0i_shear(idir)
              << " " << F_u_LOCAL(idir) << "\n";

  exit(1);
}
      for(int idir=0; idir<D; ++idir) device_hydro_vector.access(is,ia,hydro_info::du_dt, idir) = du_dt(idir);
  };
  Cabana::simd_parallel_for(simd_policy, compute_velocity_derivative, "compute_velocity_derivative");
  Kokkos::fence();

  // std::cout << "===========================================================================\n";
  // std::cout << "Particle #0 at " << __FUNCTION__ << "::" << __LINE__ << ":\n";
  // std::cout << sysPtr->particles[0] << std::endl;
  // sysPtr->print_neighbors(0);


  //Calculate the derivatives and quantities that depends on du/dt:
  //dS/dt, dN_{B,S,Q}/dt, dextensive_bulk/dt,d_extensive_shv/dt, theta and shv_nabla_u
  //using the MRF formalism and the previously calculated du/dt
  auto compute_derivatives = KOKKOS_LAMBDA(const int is, const int ia)
  {

    //read relevant quantities
    double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double divV = 0;

    //declare caches
    milne::Vector<double,D> v, du_dt, u,u_cov;
    milne::Vector<double,D> M_extensive_entropy, M_extensive_bulk;
    milne::Vector<double,3> F_extensive_N, R_extensive_entropy, R_extensive_bulk;
    milne::Matrix<double,4,4> shv;
    //fill caches
    for(int idir=0; idir<D; ++idir){
      v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
      du_dt(idir) = device_hydro_vector.access(is, ia, hydro_info::du_dt, idir);
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      divV += device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, idir);
      M_extensive_entropy(idir) = device_hydro_vector.access(is, ia, hydro_info::M_extensive_entropy, idir);
      M_extensive_bulk(idir) = device_hydro_vector.access(is, ia, hydro_info::M_extensive_bulk, idir);
    }
    for(int icharge=0; icharge<3; ++icharge){
      F_extensive_N(icharge) = device_hydro_vector.access(is, ia, hydro_info::F_extensive_N, icharge);
      R_extensive_entropy(icharge) = device_hydro_vector.access(is, ia, hydro_info::R_extensive_entropy, icharge);
      R_extensive_bulk(icharge) = device_hydro_vector.access(is, ia, hydro_info::R_extensive_bulk, icharge);
    }
    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      shv(idir, jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
    double geometric_factor = delta_i_eta(D-1)*u(D-1)*u(D-1)*t/gamma;
    u_cov = u;
    u_cov.make_covariant(t2);

    //Calculate the extensive charges derivatives
    milne::Vector<double,3> MU_aux;
    for(int icharge=0; icharge<3; ++icharge){
      MU_aux(icharge) = 0.;
      for(int idir=0; idir<D; ++idir){
        MU_aux(icharge) +=device_hydro_space_matrix.access(is, ia, hydro_info::M_extensive_N, idir, icharge)*du_dt(idir);
      }
    }
    milne::Matrix<double,3,3> R_N_inv;
    for(int icharge=0; icharge<3; ++icharge){
      for(int jcharge=0; jcharge<3; ++jcharge){
        R_N_inv(icharge,jcharge) = device_hydro_space_matrix.access(is, ia, hydro_info::R_extensive_N, icharge, jcharge);
      }
    }
    R_N_inv = milne::inverse(R_N_inv);
    milne::Vector<double,3> dN_dt = R_N_inv*(F_extensive_N + MU_aux);

    // time derivative of ``extensive (extensive) entropy density per particle"
    double d_dt_extensive_s = milne::contract(M_extensive_entropy,du_dt)
                            +device_hydro_scalar.access(is, ia, hydro_info::F_extensive_entropy)
                            +milne::contract(R_extensive_entropy,dN_dt);
    //if(d_dt_extensive_s > 1e-10){
    //  std::cout << "d_dt_extensive_s: " << d_dt_extensive_s << std::endl;
    //}
if (std::isnan(d_dt_extensive_s))
{
  std::cout << "===========================================================================\n";
  std::cout << "Particle #0 at " << __FUNCTION__ << "::" << __LINE__ << ":\n";
  std::cout << milne::contract(M_extensive_entropy,du_dt) << std::endl;
  std::cout << device_hydro_scalar.access(is, ia, hydro_info::F_extensive_entropy) << std::endl;
  std::cout << milne::contract(R_extensive_entropy,dN_dt) << std::endl;
  for(int idir=0; idir<D; ++idir)
    std::cout << "\t"
              << M_extensive_entropy(idir) << "  "
              << du_dt(idir) << "  "
              << R_extensive_entropy(idir) << "  "
              << dN_dt(idir) << std::endl;

  std::cout << "Crashing!!\n";
  exit(1);
}

	  //time derivative of the extensive bulk pressure
    double d_extensive_bulk_dt = milne::contract(M_extensive_bulk,du_dt)
                    +device_hydro_scalar.access(is, ia, hydro_info::F_extensive_bulk)
                    +milne::contract(R_extensive_bulk,dN_dt);
    //calculate the shv_nabla_u = \pi^{ij} \nabla_i u_j
    milne::Vector<double,D> M_shv_nabla_u;
    for(int idir=0; idir<D; ++idir){
      M_shv_nabla_u(idir) = device_hydro_vector.access(is, ia, hydro_info::M_shv_nabla_u, idir);
    }
    double shv_nabla_u = milne::contract(M_shv_nabla_u,du_dt)
                        +device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u);

    //theta equation
    double theta = -1.*milne::contract(u_cov,du_dt)/gamma + gamma*divV + gamma/t + geometric_factor;

if (std::isnan(d_dt_extensive_s))
{
  std::cout << "===========================================================================\n";
  std::cout << "Particle #0 at " << __FUNCTION__ << "::" << __LINE__ << ":\n";
  std::cout << sysPtr->particles[0] << std::endl;
  sysPtr->print_neighbors(0);

  std::cout << "Should have already crashed!!\n";
  exit(1);
}

    //stores the results
    device_hydro_scalar.access(is, ia, hydro_info::d_extensive_bulk_dt) = d_extensive_bulk_dt;
    device_hydro_scalar.access(is, ia, hydro_info::theta) = theta;
    device_hydro_scalar.access(is, ia, hydro_info::shv_nabla_u) = shv_nabla_u;
    device_d_dt_extensive.access(is, ia, densities_info::s) = d_dt_extensive_s;
    device_d_dt_extensive.access(is, ia, densities_info::rhoB) = dN_dt(0);
    device_d_dt_extensive.access(is, ia, densities_info::rhoS) = dN_dt(1);
    device_d_dt_extensive.access(is, ia, densities_info::rhoQ) = dN_dt(2);

if (std::isnan(d_dt_extensive_s))
{
  std::cout << "===========================================================================\n";
  std::cout << "Particle #0 at " << __FUNCTION__ << "::" << __LINE__ << ":\n";
  std::cout << sysPtr->particles[0] << std::endl;
  sysPtr->print_neighbors(0);

  std::cout << "Should have already crashed!!\n";
  exit(1);
}
    //calculate the extensive shear tensor derivative
    if(using_shear){
      milne::Matrix<double,2,3> d_extensive_shv_dt;
      milne::Matrix<double,2,3> sigma_tensor;
      for(int idir=0; idir<2; ++idir){
        for(int jdir=idir; jdir<3; ++jdir){
          double M_du_aux = 0.;
          double R_dn_aux = 0.;
          double Msigma_du_aux = 0.;
          for(int kdir=0; kdir<D; ++kdir){
            int linear_index_M = jdir * D + kdir;
            M_du_aux += device_hydro_shear_aux_matrix.access(is, ia, hydro_info::M_extensive_shear, idir,linear_index_M)*du_dt(kdir);
            Msigma_du_aux += device_hydro_shear_aux_matrix.access(is, ia, hydro_info::M_sigma_tensor, idir,linear_index_M)*du_dt(kdir);
          }
          for(int icharge=0; icharge<3; ++icharge){
            int linear_index_R = jdir * 3 + icharge;
            R_dn_aux += device_hydro_shear_aux_matrix.access(is, ia, hydro_info::R_extensive_shear, idir, linear_index_R)*dN_dt(icharge);
          }
          d_extensive_shv_dt(idir,jdir) = M_du_aux
                                  + device_hydro_shear_aux_vector.access(is, ia, hydro_info::F_extensive_shear, idir, jdir);
          device_hydro_shear_aux_vector.access(is, ia, hydro_info::d_extensive_shv_dt, idir,jdir) = d_extensive_shv_dt(idir,jdir);
          sigma_tensor(idir,jdir) = Msigma_du_aux
                                 + device_hydro_shear_aux_vector.access(is, ia, hydro_info::F_sigma_tensor, idir, jdir);
          device_hydro_spacetime_matrix.access(is, ia, hydro_info::sigma_tensor, idir+1,jdir+1) = sigma_tensor(idir,jdir);
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

  std::cout << "===========================================================================\n";
  std::cout << "Particle #0 at " << __FUNCTION__ << "::" << __LINE__ << ":\n";
  std::cout << sysPtr->particles[0] << std::endl;
  sysPtr->print_neighbors(0);


      // computing dEz_dt
  auto compute_Ez_derivative = KOKKOS_LAMBDA(const int is, const int ia)
    {
      double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);
      double shv33 = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv,3,3);
      double sigma_lab = device_hydro_scalar.access(is, ia, hydro_info::sigma_lab);
      double sph_mass = device_sph_mass.access(is, ia, densities_info::s);
      double p = device_thermo.access(is, ia, thermo_info::p);
      double e = device_thermo.access(is, ia, thermo_info::e);
      double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
      milne::Vector<double,D> u;
      for(int idir=0; idir<D; ++idir){
        u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      }

      milne::Vector<double,4> contra_metric_diag = milne::get_contra_metric_diagonal<3>(t);
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
      dEz = ((e + p + bulk)*u3*u3 + (p + bulk)/t2 + shv33 ) * sph_mass * t2 / sigma_lab;
      //double dEz = 0;
      //dEz = (p/t2) * sph_mass* t2  / sigma_lab;
      device_contribution_to_total_dEz.access(is, ia) = dEz;
    };
    Cabana::simd_parallel_for(simd_policy, compute_Ez_derivative, "compute_Ez_derivative");
    Kokkos::fence();

    // if using shear calculate shear knudsen number
    if(using_shear){
    compute_hydro_numbers(sysPtr);
    }

  std::cout << "===========================================================================\n";
  std::cout << "Particle #0 at " << __FUNCTION__ << "::" << __LINE__ << ":\n";
  std::cout << sysPtr->particles[0] << std::endl;
  sysPtr->print_neighbors(0);


  };


/// @brief Calculate the MRs for the shear tensor
/// @details This function calculates all the MRFs necessary
//  when we have shear. This includes the MRFs that appear
/// in the four-acceleration equation  M_0i_shear, R_0i_shear
/// F_0i_shear, the MRFs used for calculating the entropy
/// production (from shv_nabla_u = \pi^{ij} \nabla_i u_j)
/// M_shv_nabla_u, F_shv_nabla_u, and the MRFs used for the
/// extensive shear tensor evolution M_extensive_shear, R_extensive_shear,
/// and F_extensive_shear.
/// @param sysPtr A shared pointer to the system state
template <unsigned int D>
void EoM_default<D>::calculate_MRF_shear(std::shared_ptr<SystemState<D>> sysPtr)
{
  double t = (sysPtr->t);
  double t2 = t*t;
  //auxiliary metric variables
  milne::Vector<double,D> delta_i_eta = milne::delta_i_eta<D>();
  milne::Vector<double,4> contra_metric_diag = milne::get_contra_metric_diagonal<3>(t);

  CREATE_VIEW(device_,sysPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());
  auto compute_MRF_shear = KOKKOS_LAMBDA(const int is, const int ia)
  {
    //read relevant quantities
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

    //declare caches
    double ueta = 0;
    milne::Vector<double,D> u, u_cov, M_shv_nabla_u ,F_0i_shear;
    milne::Vector<double,D> v, grad_u0;
    //use fixed size u (ux, uy, ueta) to avoid unnecessary specializations
    milne::Vector<double,3> fixed_size_u, fixed_size_u_cov;
    milne::Matrix<double,D,3> R_0i_shear;
    milne::Matrix<double,2,3> F_extensive_shear;
    milne::Matrix<double,D,D> gradV, grad_uj, M_0i_shear;
    milne::Matrix<double,4,4> shv_hybrid, shv, shv_cov;
    milne::Matrix3D<double,2,3,D> M_extensive_shear;
    milne::Matrix3D<double, 2, 3, 3> R_extensive_shear;
    //fill caches
    for(int idir=0; idir<D; ++idir){
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      fixed_size_u(idir) = u(idir);
      v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
    }
    //fill remaning dimensions
    for(int idir=D; idir<3; ++idir){
      fixed_size_u(idir) = 0.0;
    }

    for(int idir=0; idir<4; ++idir)
    for(int jdir=0; jdir<4; ++jdir)
      shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
    u_cov = u;
    u_cov.make_covariant(t2);
    fixed_size_u_cov = fixed_size_u;
    fixed_size_u_cov.make_covariant(t2);
    double geometric_factor = delta_i_eta(D-1)*u(D-1)*u(D-1)*t/gamma;
    shv_hybrid = shv;
    shv_hybrid.make_covariant(1, t2);

    shv_cov = shv_hybrid;
    shv_cov.make_covariant(0, t2);

    //initialize the auxiliary Ms and Fs
    milne::Vector<double,D> F_i0_sigma;
    milne::Vector<double,D> F_i0_D;
    milne::Vector<double,D> F_i0_domega;
    milne::Vector<double,D> F_i0_dsigma;
    milne::Vector<double,D> F_i0_dd;
    milne::Matrix<double,D,D> M_i0_sigma;
    milne::Matrix<double,D,D> M_i0_D;
    milne::Matrix<double,D,D> M_i0_domega;
    milne::Matrix<double,D,D> M_i0_dsigma;
    milne::Matrix<double,D,D> M_i0_dd;

    milne::Matrix<double, 2, 3> F_D;
    milne::Matrix<double, 2, 3> F_sigma;
    milne::Matrix<double, 2, 3> F_delta;
    milne::Matrix<double, 2, 3> F_domega;
    milne::Matrix<double, 2, 3> F_dsigma;
    milne::Matrix3D<double, 2, 3, D> M_D;
    milne::Matrix3D<double, 2, 3, D> M_sigma;
    milne::Matrix3D<double, 2, 3, D> M_delta;
    milne::Matrix3D<double, 2, 3, D> M_domega;
    milne::Matrix3D<double, 2, 3, D> M_dsigma;


    //calculate aux quantities
    double divV = 0;
    for(int idir=0; idir<D; ++idir){
      divV += device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, idir);
    for(int jdir=0; jdir<D; ++jdir){
      gradV(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradV, idir, jdir);
      }
    }

    for (int idir=0; idir<D; ++idir){
      for (int jdir=0; jdir<D; ++jdir){
        grad_uj(idir,jdir) = gamma*gradV(idir,jdir);
        for (int kdir=0; kdir<D; ++kdir){
          grad_uj(idir,jdir) += -gamma*u(jdir)*u_cov(kdir)*gradV(idir,kdir);
        }
      }
    }

    grad_u0 = -1.*milne::contract(grad_uj,u_cov,milne::SecondIndex())/gamma;
    milne::Vector<double,D> grad_u0_contra = grad_u0;
    grad_u0_contra.make_contravariant(t2);

    milne::Matrix<double,D,D> contra_grad_uj;
    contra_grad_uj = grad_uj;
    contra_grad_uj.make_contravariant(0, t2);
    // use fixed size grad_uj to avoid unnecessary specializations
    milne::Matrix<double,3,3> contra_grad_uj_3d;
    for(int idir=0; idir<D; ++idir){
      for(int jdir=0; jdir<D; ++jdir){
        contra_grad_uj_3d(idir,jdir) = contra_grad_uj(idir,jdir);
      }
    }
    //fill remaning dimensions
    for(int idir=D; idir<3; ++idir){
      for(int jdir=0; jdir<3; ++jdir){
        contra_grad_uj_3d(idir,jdir) = 0.0;
      }
    }
    milne::Vector<double,4> shear0mu_aux;
    for(int idir=0; idir<4; ++idir){
      shear0mu_aux(idir) = shv_hybrid(0,idir);
    }

    //aux vectors and matrices for the four-acceleration equation
    for(int idir=0; idir<D; ++idir){
      F_i0_sigma(idir) =  -(milne::contract(grad_uj,v,milne::FirstIndex())(idir)
                          -grad_u0_contra(idir))/2.
                          +u(idir)*gamma*(geometric_factor + gamma/t + gamma*divV)/3.
                          -gamma*geometric_factor*u(idir)/2.;
      F_i0_D(idir) = gamma*(u(idir)*shv_hybrid(0,0)+gamma*shv_hybrid(idir+1,0))*geometric_factor
                     +t*u(D-1)*shv(3,idir+1)*delta_i_eta(idir);
      F_i0_domega(idir) = 0.0;
      F_i0_dsigma(idir) = 0.0;
      F_i0_dd(idir) = shv(idir+1,0)*(geometric_factor+gamma/t+gamma*divV);
      for(int jdir=0; jdir<D; ++jdir){
        M_i0_sigma(idir,jdir) = u(idir)*u_cov(jdir)/(2.*3.);
        M_i0_D(idir,jdir) = gamma*u(idir)*(shv_hybrid(0,jdir+1)-shv_hybrid(0,0)*u_cov(jdir)/gamma)
                            +gamma*gamma*(shv_hybrid(idir+1,jdir+1)-shv_hybrid(idir+1,0)*u_cov(jdir)/gamma);
        M_i0_domega(idir,jdir) = 0.0;
        M_i0_dsigma(idir,jdir) = 0.0;
        M_i0_dd(idir,jdir) = -shv(idir+1,0)*u_cov(jdir)/gamma;
      }
      //diagonal terms
      M_i0_sigma(idir,idir) += (1.-gamma*gamma)/2.;
      //metric only terms
      F_i0_sigma(idir) += delta_i_eta(idir)*(-(-u(D-1)*t-u(D-1)/t
                          +gamma*u(D-1)*u(D-1)*t
                          +2.*gamma*gamma*u(D-1)/t)/2.);
      F_i0_D(idir) += delta_i_eta(idir)*(u(D-1)*shv(0,0)/t
                      + gamma*shv(3,0)/t);
    }

    //calculate M and F for shv_nabla_u = \pi^{ij} \nabla_i u_j
    double F_shv_nabla_u = 0.0;
    for (int idir = 0; idir < D; ++idir) {
      M_shv_nabla_u(idir) = shv_hybrid(0,idir+1) - shv_hybrid(0,0)*u_cov(idir)/gamma;
      device_hydro_vector.access(is, ia, hydro_info::M_shv_nabla_u, idir) = M_shv_nabla_u(idir);
      for (int jdir = 0; jdir < D; ++jdir) {
          F_shv_nabla_u += (shv_hybrid(idir + 1, jdir + 1)
                           - v(idir) * shv_hybrid(0, jdir + 1)) * grad_uj(idir, jdir);
      }
      F_shv_nabla_u += (shv_hybrid(idir + 1, 0) - v(idir) * shv_hybrid(0, 0)) * grad_u0(idir);
    }
    F_shv_nabla_u += +shv_hybrid(3,3)*gamma/t + shv_hybrid(0,0)*geometric_factor
            +(t*shv_hybrid(3,0)+shv_hybrid(0,3)/t)*u(D-1)*delta_i_eta(D-1);
    device_hydro_scalar.access(is, ia, hydro_info::F_shv_nabla_u) = F_shv_nabla_u;
    //calculate the MRFs that will be used in the four-acceleration equation
    for(int idir=0; idir<D; ++idir){
      F_0i_shear(idir) = (-shv(idir+1,0)
                    -tau_pi*F_i0_D(idir)
                    -delta_pipi*F_i0_dd(idir)
                    -a*2.*tau_pi*F_i0_domega(idir)
                    -a*tau_pipi*F_i0_dsigma(idir)
                    +(2.*eta_pi+a*lambda_piPi*bulk)*F_i0_sigma(idir)
                    +a*phi6*bulk*shv(idir+1,0)
                    +a*phi7*(milne::contract(shv, shear0mu_aux, milne::SecondIndex())(idir+1)
                    +gamma*u(idir)*milne::contract(shv_cov,shv)/3.))/(tau_pi*gamma);
      //stores the results
      device_hydro_vector.access(is, ia, hydro_info::F_0i_shear, idir) = F_0i_shear(idir);
      for(int jdir=0; jdir<D; ++jdir){
        M_0i_shear(idir,jdir) = (-tau_pi*M_i0_D(idir,jdir)
                                -delta_pipi*M_i0_dd(idir,jdir)
                                -2.*a*M_i0_domega(idir,jdir)
                                +a*tau_pipi*M_i0_dsigma(idir,jdir)
                                +(2.*eta_pi+a*lambda_piPi*bulk)*M_i0_sigma(idir,jdir))/(tau_pi*gamma);
        //stores the results
        device_hydro_space_matrix.access(is, ia, hydro_info::M_0i_shear, idir, jdir) = M_0i_shear(idir,jdir);
      }
      //calculate the Rs
      for(int icharge=0; icharge<3; ++icharge){
        //TODO: change this when using diffusion coupling
        R_0i_shear(idir,icharge) = 0.0;
        device_hydro_space_matrix.access(is, ia, hydro_info::R_0i_shear, idir, icharge) = R_0i_shear(idir,icharge);
      }
    }

    //calculate the aux MRFs that will be used in the evolution of the shear tensor
    for(int idir=0; idir<2; ++idir){
      for(int jdir=idir; jdir<3; ++jdir){
        F_D(idir,jdir) = gamma*geometric_factor*(fixed_size_u(idir)*shv_hybrid(jdir+1,0)
                  +fixed_size_u(jdir)*shv_hybrid(idir+1,0));

        F_sigma(idir,jdir) = contra_grad_uj_3d(idir,jdir)/2. + contra_grad_uj_3d(jdir,idir)/2.
                  +fixed_size_u(idir)*fixed_size_u(jdir)*(geometric_factor+gamma/t+gamma*divV)/3.;
        F_delta(idir,jdir) = shv(idir+1,jdir+1)*(geometric_factor+gamma/t+gamma*divV);
        F_domega(idir,jdir) = 0.0;
        F_dsigma(idir,jdir) = 0.0;

        for(int kdir=0; kdir<D; ++kdir){
          M_D(idir,jdir,kdir) = gamma*(fixed_size_u(idir)*shv_hybrid(jdir+1,kdir+1)
                                      -fixed_size_u(idir)*shv_hybrid(jdir+1,0)*fixed_size_u_cov(kdir)/gamma
                                      +fixed_size_u(jdir)*shv_hybrid(idir+1,kdir+1)
                                      -fixed_size_u(jdir)*shv_hybrid(idir+1,0)*fixed_size_u_cov(kdir)/gamma);
          M_sigma(idir,jdir,kdir) = -fixed_size_u(idir)*fixed_size_u(jdir)*fixed_size_u_cov(kdir)/(3.*gamma);
          M_delta(idir,jdir,kdir) = -shv(idir+1,jdir+1)*fixed_size_u_cov(kdir)/gamma;
          M_domega(idir,jdir,kdir) = 0.0;
          M_dsigma(idir,jdir,kdir) = 0.0;
        }
        //diagonal i = k terms
        //the if avoids the case when i>k (for D=1)
        if (D==1){
          ueta = fixed_size_u(0);
          M_sigma(0,jdir,0) += -gamma*fixed_size_u(jdir)/2.;
        }
        else{
          M_sigma(idir,jdir,idir) += -gamma*fixed_size_u(jdir)/2.;
          ueta = fixed_size_u(2);
        }
      }
      //diagonal j=k terms
      for(int jdir = idir; jdir<D; ++jdir){
        //j=k , since k<D , we only loop until j<D
        M_sigma(idir,jdir,jdir) += -gamma*fixed_size_u(idir)/2.;
      }
      //diagonal i=j terms
      for(int kdir=0; kdir<D; ++kdir){
        //j=i
        M_sigma(idir,idir,kdir) += contra_metric_diag(idir+1)*fixed_size_u_cov(kdir)/(gamma*3.);
      }

      //diagonal and metric terms for F
      F_D(idir,2) += gamma*shv(3,idir+1)/t
                     +ueta*shv(0,idir+1)/t;
      F_sigma(idir,idir) += -contra_metric_diag(idir+1)*(geometric_factor+gamma/t+gamma*divV)/3.;
      F_sigma(idir,2) += -fixed_size_u(idir)*ueta*ueta*t/2.
                         -fixed_size_u(idir)*ueta*gamma/t;

    }

    for(int idir=0; idir<2; ++idir){
      for(int jdir=idir; jdir<3; ++jdir){
        F_extensive_shear(idir,jdir) = (-shv(idir+1,jdir+1)
                                -tilde_delta*F_delta(idir,jdir)
                                -tau_pi*F_D(idir,jdir)
                                -a*2.*tau_pi*F_domega(idir,jdir)
                                -a*tau_pipi*F_dsigma(idir,jdir)
                                +(2.*eta_pi+a*lambda_piPi*bulk)*F_sigma(idir,jdir)
                                +a*phi6*bulk*shv(idir+1,jdir+1)
                                +a*phi7*(milne::contract(shv,shv_hybrid,milne::SecondIndex(),milne::SecondIndex())(idir,jdir)
                                +fixed_size_u(idir)*fixed_size_u(jdir)*milne::contract(shv,shv_cov)/3.))/(sigma*tau_pi*gamma);

        for(int kdir=0; kdir<D; ++kdir){
          M_extensive_shear(idir,jdir,kdir) =(-tau_pi*M_D(idir,jdir,kdir)
                                        -tilde_delta*M_delta(idir,jdir,kdir)
                                        -2.*a*M_domega(idir,jdir,kdir)
                                        -a*tau_pipi*M_dsigma(idir,jdir,kdir)
                                        +(2.*eta_pi+a*lambda_piPi*bulk)*M_sigma(idir,jdir,kdir))/(sigma*tau_pi*gamma);
          //saves the results
          int linear_index = jdir * D + kdir;
          device_hydro_shear_aux_matrix.access(is, ia, hydro_info::M_extensive_shear, idir, linear_index) = M_extensive_shear(idir,jdir,kdir);
          device_hydro_shear_aux_matrix.access(is, ia, hydro_info::M_sigma_tensor, idir, linear_index) = M_sigma(idir,jdir,kdir);
          }
        //calculate the Rs
        for(int icharge=0; icharge<3; ++icharge){
          //TODO: change this when using diffusion coupling
          R_extensive_shear(idir,jdir,icharge) = 0.0;
          int linear_index = jdir * 3 + icharge;
          device_hydro_shear_aux_matrix.access(is, ia, hydro_info::R_extensive_shear, idir, linear_index) = R_extensive_shear(idir,jdir,icharge);
          }
      }
      F_extensive_shear(idir,idir) += a*phi7*(-contra_metric_diag(idir+1)*milne::contract(shv,shv_cov))/(3.*tau_pi*gamma*sigma);
      //saves the results
      for(int jdir=idir; jdir<3; ++jdir){
        device_hydro_shear_aux_vector.access(is, ia, hydro_info::F_extensive_shear, idir, jdir) = F_extensive_shear(idir,jdir);
        device_hydro_shear_aux_vector.access(is, ia, hydro_info::F_sigma_tensor, idir, jdir) = F_sigma(idir,jdir);
      }
    }
  };
  Cabana::simd_parallel_for(simd_policy, compute_MRF_shear, "compute_MRF_shear");
  Kokkos::fence();
};


/// @brief Calculate the knudsen number
/// @details This function calculates the knudsen number
/// for each particle in the system. The knudsen number
/// is defined as for the shear tensor as
/// Kn = tau_pi * abs(sigma_mu_nu sigma^{mu nu})
/// @param sysPtr A shared pointer to the system state
template <unsigned int D>
void EoM_default<D>::compute_hydro_numbers(std::shared_ptr<SystemState<D>> sysPtr)
{
  double t = (sysPtr->t);
  double t2 = t*t;
  //auxiliary metric variables
  milne::Vector<double,D> delta_i_eta = milne::delta_i_eta<D>();
  milne::Vector<double,4> contra_metric_diag = milne::get_contra_metric_diagonal<3>(t);

  CREATE_VIEW(device_,sysPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());
  auto caus = KOKKOS_LAMBDA(const int is, const int ia)
  {
    //read relevant quantities
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double eta_pi = device_hydro_scalar.access(is, ia, hydro_info::eta_pi);
    double phi6 = device_hydro_scalar.access(is, ia, hydro_info::phi6);
    double phi7 = device_hydro_scalar.access(is, ia, hydro_info::phi7);
    double lambda_piPi = device_hydro_scalar.access(is, ia, hydro_info::lambda_piPi);
    double zeta_Pi = device_hydro_scalar.access(is, ia, hydro_info::zeta_Pi);
    double tau_Pi = device_hydro_scalar.access(is, ia, hydro_info::tau_Pi);
    double delta_PiPi = device_hydro_scalar.access(is, ia, hydro_info::delta_PiPi);
    double phi1 = device_hydro_scalar.access(is, ia, hydro_info::phi1);
    double phi3 = device_hydro_scalar.access(is, ia, hydro_info::phi3);
    double lambda_Pipi = device_hydro_scalar.access(is, ia, hydro_info::lambda_Pipi);
    double tau_pipi = device_hydro_scalar.access(is, ia, hydro_info::tau_pipi);
    double a= device_hydro_scalar.access(is, ia, hydro_info::a);
    double delta_pipi = device_hydro_scalar.access(is, ia, hydro_info::delta_pipi);
    double tau_pi = device_hydro_scalar.access(is, ia, hydro_info::tau_pi);
    double delta_Pi = device_hydro_scalar.access(is, ia, hydro_info::delta_PiPi);
    double bulk = device_hydro_scalar.access(is, ia, hydro_info::bulk);
    double theta = device_hydro_scalar.access(is, ia, hydro_info::theta);
    double pressure = device_thermo.access(is, ia, thermo_info::p);
    //declare caches
    milne::Vector<double,D> u, u_cov;
    //use fixed size u (ux, uy, ueta) to avoid unnecessary specializations
    milne::Vector<double,3> fixed_size_u, fixed_size_u_cov;
    milne::Matrix<double,4,4> shv, shv_hybrid, shv_cov;
    //fill caches
    for(int idir=0; idir<D; ++idir){
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
      fixed_size_u(idir) = u(idir);

    }
    //fill remaning dimensions
    for(int idir=D; idir<3; ++idir){
      fixed_size_u(idir) = 0.0;
    }

    for(int idir=0; idir<4; ++idir){
    for(int jdir=0; jdir<4; ++jdir){
      shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
    }}
    u_cov = u;
    u_cov.make_covariant(t2);
    fixed_size_u_cov = fixed_size_u;
    fixed_size_u_cov.make_covariant(t2);
    shv_hybrid = shv;
    shv_hybrid.make_covariant(1, t2);
    shv_cov = shv_hybrid;
    shv_cov.make_covariant(0, t2);


    //calculate the knudsen number
    milne::Matrix<double, 4, 4> sigma_tensor;
    for(int idir=0; idir<2; ++idir){
      for(int jdir=idir; jdir<3; ++jdir){
        sigma_tensor(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::sigma_tensor, idir+1,jdir+1);
      }
    }
    //calculate the other components of the matrix
    //Symmetrizes space part

    for( int i=1; i<4; i++ )
    for( int j=i+1; j<4; j++ )
      sigma_tensor(j,i) = sigma_tensor(i,j);

    //pi^{0i} = -\pi^{ij}u_j/gamma
    for(int idir=1; idir<D+1; ++idir){
      milne::Vector<double,D> sigma_tensor_i;
      for(int jdir=1; jdir<D+1; ++jdir) sigma_tensor_i(jdir-1) = sigma_tensor(idir,jdir);
      sigma_tensor(0,idir) = -1.*milne::contract(u_cov,sigma_tensor_i)/gamma;
    }

    //Symmetrizes time part
    for( int i=1; i<4; i++ )
      sigma_tensor(i,0) = sigma_tensor(0,i);


    //pi^00 = u_i u_j pi^{ij}/gamma^2
    sigma_tensor(0,0) =( fixed_size_u_cov(0)*fixed_size_u_cov(0)*sigma_tensor(1,1)
                + fixed_size_u_cov(1)*fixed_size_u_cov(1)*sigma_tensor(2,2)
                + 2.*fixed_size_u_cov(1)*fixed_size_u_cov(0)*sigma_tensor(2,1)
                + 2.*fixed_size_u_cov(2)*fixed_size_u_cov(0)*sigma_tensor(3,1)
                + 2.*fixed_size_u_cov(2)*fixed_size_u_cov(1)*sigma_tensor(3,2)
                - fixed_size_u_cov(2)*fixed_size_u_cov(2)*(sigma_tensor(1,1)+sigma_tensor(2,2))/t2
              )/(gamma*gamma-fixed_size_u_cov(2)*fixed_size_u_cov(2)/(t2));


    //pi^33 = (pi^00 - pi^11 - pi^22)/t^2
    sigma_tensor(3,3) = (sigma_tensor(0,0) - sigma_tensor(1,1) - sigma_tensor(2,2))/t2;

    milne::Matrix<double, 4, 4> sigma_hybrid = sigma_tensor;
    sigma_hybrid.make_covariant(1, t2);
    milne::Matrix<double, 4, 4> sigma_cov = sigma_hybrid;
    sigma_cov.make_covariant(0, t2);
    double sigma_norm = milne::contract(sigma_tensor, sigma_cov);
    double shear_knudsen_number = tau_pi * sqrt(abs(sigma_norm));
    device_hydro_scalar.access(is, ia, hydro_info::shear_knudsen) = shear_knudsen_number;

    device_hydro_scalar.access(is, ia, hydro_info::bulk_knudsen) = tau_Pi*abs(theta);

    double shv_magnitude = sqrt(abs(milne::contract(shv,shv_cov)));
    device_hydro_scalar.access(is, ia, hydro_info::inverse_reynolds_shear) = shv_magnitude/pressure;
    device_hydro_scalar.access(is, ia, hydro_info::inverse_reynolds_bulk) = abs(bulk)/pressure;

  };
  Cabana::simd_parallel_for(simd_policy, caus, "compute_MRF_shear");
  Kokkos::fence();
};












/// @brief Checks causality conditions in the hydrodynamic evolution.
/// @details Computes the eigenvalues of the shear tensor \(\pi^{\mu\nu}\) and verifies
/// necessary and sufficient conditions for causality using the most general formulation.
/// Stores the results in device_hydro_scalar.
/// @param sysPtr A shared pointer to the system state.
template <unsigned int D>
void EoM_default<D>::check_causality(std::shared_ptr<SystemState<D>> sysPtr)
{
    double t = (sysPtr->t);
    double t2 = t * t;
    CREATE_VIEW(device_, sysPtr->cabana_particles);
    auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH, ExecutionSpace>(0, sysPtr->cabana_particles.size());

    auto get_sorted_eigenvalues_of_pi_mu_nu = [](milne::Matrix<double, 4, 4> &shv_hybrid,
                                                 double &Lambda_0,
                                                 double &Lambda_1, double &Lambda_2, double &Lambda_3,
                                                double &eps, const double &t2) -> bool
    {
        double m[16];
        m[0]  =  shv_hybrid(0, 0); m[1]  =  shv_hybrid(0, 1); m[2]  =  shv_hybrid(0, 2); m[3]  =  shv_hybrid(0, 3);
        m[4]  =  shv_hybrid(1, 0); m[5]  =  shv_hybrid(1, 1); m[6]  =  shv_hybrid(1, 2); m[7]  =  shv_hybrid(1, 3);
        m[8]  =  shv_hybrid(2, 0); m[9]  =  shv_hybrid(2, 1); m[10] =  shv_hybrid(2, 2); m[11] =  shv_hybrid(2, 3);
        m[12] =  shv_hybrid(3, 0); m[13] =  shv_hybrid(3, 1); m[14] =  shv_hybrid(3, 2); m[15] =  shv_hybrid(3, 3);
        //m[0]  =  -shv_hybrid(0, 0); m[1]  =  -shv_hybrid(0, 1); m[2]  =  -shv_hybrid(0, 2); m[3]  =  -shv_hybrid(0, 3);
        //m[4]  =  shv_hybrid(0, 1); m[5]  =  shv_hybrid(1, 1); m[6]  =  shv_hybrid(1, 2); m[7]  =  shv_hybrid(1, 3);
        //m[8]  =  shv_hybrid(0, 2); m[9]  =  shv_hybrid(1, 2); m[10] =  shv_hybrid(2, 2); m[11] =  t2*shv_hybrid(2, 3);
        //m[12] =  t2*shv_hybrid(0, 3); m[13] =  t2*shv_hybrid(1,3); m[14] =  t2*shv_hybrid(2,3); m[15] =  t2*shv_hybrid(3, 3);
      	gsl_vector_complex *eval = gsl_vector_complex_alloc(4);
      	gsl_matrix_complex *evec = gsl_matrix_complex_alloc(4, 4);

      	gsl_matrix_view mat = gsl_matrix_view_array(m, 4, 4);
      	gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(4);
      	int success = gsl_eigen_nonsymmv (&mat.matrix, eval, evec, w);
      	gsl_eigen_nonsymmv_free(w);

      	// sort by magnitude first
      	gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

      	// check eigensystem
      	/*if (true)
      	for (int i = 0; i < 4; i++)
            {
              gsl_complex eval_i
                 = gsl_vector_complex_get (eval, i);
              gsl_vector_complex_view evec_i
                 = gsl_matrix_complex_column (evec, i);

              printf ("eigenvalue = %g + %gi\n",
                      GSL_REAL(eval_i), GSL_IMAG(eval_i));
              printf ("eigenvector = \n");
              for (int j = 0; j < 4; ++j)
                {
                  gsl_complex z =
                    gsl_vector_complex_get(&evec_i.vector, j);
                  printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
                }
            }*/



      	for ( int elem = 0; elem < 4; elem++ )
      		if ( abs(GSL_IMAG(gsl_vector_complex_get(eval, elem)))
      				> 0.01*abs(GSL_REAL(gsl_vector_complex_get(eval, elem))) )
      		{
      			return false;
      		}

      	Lambda_0 = GSL_REAL(gsl_vector_complex_get(eval, 0));
      	double ratio = abs(Lambda_0 / (abs(GSL_REAL(gsl_vector_complex_get(eval, 3)))+eps));

      	/*cout << "Check #1 here: "
      		<< GSL_REAL(gsl_vector_complex_get(eval, 0)) << "   "
      		<< GSL_REAL(gsl_vector_complex_get(eval, 1)) << "   "
      		<< GSL_REAL(gsl_vector_complex_get(eval, 2)) << "   "
      		<< GSL_REAL(gsl_vector_complex_get(eval, 3)) << endl;*/

      	// sort by value next
      	gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

      	/*cout << "Check #2 here: "
      		<< GSL_REAL(gsl_vector_complex_get(eval, 0)) << "   "
      		<< GSL_REAL(gsl_vector_complex_get(eval, 1)) << "   "
      		<< GSL_REAL(gsl_vector_complex_get(eval, 2)) << "   "
      		<< GSL_REAL(gsl_vector_complex_get(eval, 3)) << endl;*/

      	double tmp0 = GSL_REAL(gsl_vector_complex_get(eval, 0));
      	double tmp1 = GSL_REAL(gsl_vector_complex_get(eval, 1));
      	double tmp2 = GSL_REAL(gsl_vector_complex_get(eval, 2));
      	double tmp3 = GSL_REAL(gsl_vector_complex_get(eval, 3));


      	// sort eval by values
      	//Lambda_0 = 0.0;
      	Lambda_1 = tmp0;
      	Lambda_2 = ( abs(tmp1) > abs(tmp2) ) ? tmp1 : tmp2;
      	Lambda_3 = tmp3;

      	if ( ratio > 0.01 )
      	{
      		success++;
      	}
      	//else
      	//	cerr << "Found zero eigenvalue with ratio = " << ratio << endl;

      	gsl_vector_complex_free(eval);
      	gsl_matrix_complex_free(evec);

      	//return ( ( success == 0 ) and ( ratio <= epsilon ) );
      	return ( success == 0 );
    };


    auto causality_check = KOKKOS_LAMBDA(const int is, const int ia)
    {
        double e, p, Pi, eta, zeta, tau_pi, tau_Pi, delta_PiPi, lambda_Pipi, delta_pipi, lambda_piPi, cs2,
               tau_pipi;
        double Lambda_0, Lambda_1, Lambda_2, Lambda_3;

        e  = device_thermo.access(is, ia, thermo_info::e);
        p  = device_thermo.access(is, ia, thermo_info::p);
        Pi = device_hydro_scalar.access(is, ia, hydro_info::bulk);
        eta = device_hydro_scalar.access(is, ia, hydro_info::eta_pi);
        zeta = device_hydro_scalar.access(is, ia, hydro_info::zeta_Pi);
        tau_pi = device_hydro_scalar.access(is, ia, hydro_info::tau_pi);
        tau_pipi = device_hydro_scalar.access(is, ia, hydro_info::tau_pipi);
        tau_Pi = device_hydro_scalar.access(is, ia, hydro_info::tau_Pi);
        delta_PiPi = device_hydro_scalar.access(is, ia, hydro_info::delta_PiPi);
        lambda_Pipi = device_hydro_scalar.access(is, ia, hydro_info::lambda_Pipi);
        delta_pipi = device_hydro_scalar.access(is, ia, hydro_info::delta_pipi);
        lambda_piPi = device_hydro_scalar.access(is, ia, hydro_info::lambda_piPi);
        cs2 = device_thermo.access(is, ia, ccake::thermo_info::cs2);

        milne::Matrix<double, 4, 4> shv;
        for (int idir = 0; idir < 4; ++idir)
        {
            for (int jdir = 0; jdir < 4; ++jdir)
            {
                shv(idir, jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
            }
        }
        milne::Matrix<double, 4, 4> shv_hybrid = shv;
        shv_hybrid.make_covariant(0, t2);
        //shv_mu^nu

        int causality_result;

        bool success = get_sorted_eigenvalues_of_pi_mu_nu(shv, Lambda_0 ,Lambda_1, Lambda_2, Lambda_3, e, t2);
        if (!success)
        {
            std::cout << "Error: Failed to compute eigenvalues of pi_mu_nu for particle " << ia << " in system, setting causality to -100." << std::endl;
            //print shv(mu,3) values
            for (int idir = 0; idir < 4; ++idir)
            {
                std::cout << "shv(" << idir << ",3) = " << shv(idir, 3) << std::endl;
            }
            causality_result = -100; // Error code

        }

        std::array<double, 3> Lambda = {Lambda_1, Lambda_2, Lambda_3};

        //basic conditions
				bool necessary_conditions
						= (tau_Pi>=0) && (tau_pi>=0) && (eta>=0) && (zeta>=0)
							&& (tau_pipi>=0) && (delta_PiPi>=0) && (lambda_Pipi>=0)
							&& (delta_pipi>=0) && (lambda_piPi>=0) && (cs2>=0)
							&& (e>=0) && (p>=0) && (e+p+Pi>0);

        // Evaluate conditions that do not require loops
        bool condition4a = (2. * eta + lambda_piPi * Pi) - (tau_pi * fabs(Lambda_1) / 2.) >= 0;

        bool condition4b = (e + p + Pi - (1.0 / (2. * tau_pi)) * (2. * eta + lambda_piPi * Pi) -
                           (tau_pipi / (4. * tau_pi)) * Lambda_3) >= 0;

        // Check conditions that require a single index "d"
        for (int d = 1; d <= 3; ++d) {
            double Lambda_d = Lambda[d - 1];  // Get Lambda_d
            bool condition4e = (1.0 / (2. * tau_pi)) * (2. * eta + lambda_piPi * Pi) +
                              (tau_pipi / (2. * tau_pi)) * Lambda_d +
                              (1.0 / (6. * tau_pi)) * (2. * eta + lambda_piPi * Pi +
                                                       (6 * delta_pipi - tau_pipi) * Lambda_d) +
                              (zeta + delta_PiPi * Pi + lambda_Pipi * Lambda_d) / tau_Pi +
                              (e + p + Pi + Lambda_d) * cs2 >= 0;

            bool condition4f = (e + p + Pi + Lambda_d - (1.0 / (2. * tau_pi)) * (2. * eta + lambda_piPi * Pi) -
                               (tau_pipi / (2. * tau_pi)) * Lambda_d -
                               (1.0 / (6. * tau_pi)) * (2. * eta + lambda_piPi * Pi +
                                                        (6 * delta_pipi - tau_pipi) * Lambda_d) -
                               (zeta + delta_PiPi * Pi + lambda_Pipi * Lambda_d) / tau_Pi -
                               (e + p + Pi + Lambda_d) * cs2) >= 0;

            necessary_conditions &= (condition4e && condition4e);
        }

        // Check conditions that require two indices "a, d" with a ≠ d
        for (int a = 1; a <= 3; ++a) {
            for (int d = 1; d <= 3; ++d) {
                if (a == d) continue; // Skip cases where a == d

                double Lambda_a = Lambda[a - 1];  // Get Lambda_a
                double Lambda_d = Lambda[d - 1];  // Get Lambda_d

                bool condition4c = (1.0 / (2. * tau_pi)) * (2. * eta + lambda_piPi * Pi) +
                                  (tau_pipi / (4. * tau_pi)) * (Lambda_a + Lambda_d) >= 0;

                bool condition4d = (e + p + Pi + Lambda_a - (1.0 / (2. * tau_pi)) * (2. * eta + lambda_piPi * Pi) -
                                   (tau_pipi / (4. * tau_pi)) * (Lambda_d + Lambda_a)) >= 0;

                necessary_conditions &= (condition4c && condition4d);
            }
        }

        // Combine the results
        necessary_conditions &= (condition4a && condition4b);
        if (!necessary_conditions)
        {
          causality_result = -1; //not causal

        }
        else{
            // Condition 5a
            bool condition5a =
            (e + p + Pi - fabs(Lambda_1)) - (1.0 / (2. * tau_pi)) * (2. * eta + lambda_piPi * Pi) -
            (tau_pipi / (2. * tau_pi)) * Lambda_3 >= 0;

            bool condition5b =
            (2. * eta + lambda_piPi * Pi) - tau_pipi * fabs(Lambda_1) > 0;
            bool condition5c =
            tau_pipi <= 6. * delta_pipi;
            bool condition5d =
            (lambda_Pipi / tau_Pi + cs2 - (tau_pipi / (12. * tau_pi))) >= 0;
            // Condition 5e
            bool condition5e =
            (1.0 / (3. * tau_pi)) * (4. * eta + 2. * lambda_piPi * Pi + (3. * delta_pipi + tau_pipi) * Lambda_3) +
            (zeta + delta_PiPi * Pi + lambda_Pipi * Lambda_3) / tau_Pi +
            fabs(Lambda_1) + Lambda_3 * cs2 +
            ( (12. * delta_pipi - tau_pipi) / (12. * tau_pi) ) *
            ( (lambda_Pipi / tau_Pi + cs2 - tau_pipi / (12. * tau_pi)) * pow(Lambda_3 + fabs(Lambda_1), 2.) ) /
            ( (e + p + Pi - fabs(Lambda_1)) - (1.0 / (2. * tau_pi)) * (2. * eta + lambda_piPi * Pi) - (tau_pipi / (2. * tau_pi)) * Lambda_3 )
            <= (e + p + Pi) * (1. - cs2);

            // Condition 5f
            bool condition5f =
            (1.0 / (6. * tau_pi)) * (2. * eta + lambda_piPi * Pi + (tau_pipi - 6. * delta_pipi) * fabs(Lambda_1)) +
            (zeta + delta_PiPi * Pi - lambda_Pipi * fabs(Lambda_1)) / tau_Pi +
            (e + p + Pi - fabs(Lambda_1)) * cs2 >= 0;

            // Condition 5g
            bool condition5g =
            1. >= ( (12. * delta_pipi - tau_pipi) / (12. * tau_pi) ) *
            ( (lambda_Pipi / tau_Pi + cs2 - tau_pipi / (12. * tau_pi)) * pow(Lambda_3 + fabs(Lambda_1), 2.) ) /
            pow( (1.0 / (2. * tau_pi)) * (2 * eta + lambda_piPi * Pi) - (tau_pipi / (2. * tau_pi)) * fabs(Lambda_1), 2. );


            bool condition5h =
            (1.0 / (3. * tau_pi)) * (4. * eta + 2. * lambda_piPi * Pi - (3 * delta_pipi + tau_pipi) * fabs(Lambda_1)) +
            (zeta + delta_PiPi * Pi - lambda_Pipi * fabs(Lambda_1)) / tau_Pi +
            (e + p + Pi - fabs(Lambda_1)) * cs2 >=
            ( (e + p + Pi + Lambda_2) * (e + p + Pi + Lambda_3) ) /
            ( 3. * (e + p + Pi - fabs(Lambda_1)) ) *
            ( 1. + 2. * ( (1.0 / (2. * tau_pi)) * (2. * eta + lambda_piPi * Pi) + (tau_pipi / (2. * tau_pi)) * Lambda_3 ) /
            ( e + p + Pi - fabs(Lambda_1) ) );

            //condition5h = true; // Temporarily set to true for testing

            bool sufficient_conditions =
                (condition5a && condition5b && condition5c && condition5d && condition5e && condition5f && condition5g && condition5h);
            causality_result = sufficient_conditions ? 1 : 0; // 1: Causal, 0: Not determined
            //std::cout << "Causality result for particle " << ia << ": " << causality_result << std::endl;
            //write which conditions failed
            //if (!sufficient_conditions)
            //{
            //  std::cout << "Failed conditions: " << std::endl;
            //  if (!condition5a) std::cout << "Condition 5a failed." << std::endl;
            //  if (!condition5b) std::cout << "Condition 5b failed." << std::endl;
            //  if (!condition5c) std::cout << "Condition 5c failed." << std::endl;
            //  if (!condition5d) std::cout << "Condition 5d failed." << std::endl;
            //  if (!condition5e) std::cout << "Condition 5e failed." << std::endl;
            //  if (!condition5f) std::cout << "Condition 5f failed." << std::endl;
            //  if (!condition5g) std::cout << "Condition 5g failed." << std::endl;
            //  if (!condition5h) std::cout << "Condition 5h failed." << std::endl;
            //}

        }
        device_hydro_scalar.access(is, ia, hydro_info::causality) = causality_result;
  };

    Cabana::simd_parallel_for(simd_policy, causality_check, "check_causality_conditions");
    Kokkos::fence();
}
