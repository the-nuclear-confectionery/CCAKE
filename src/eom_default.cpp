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
  double t2 = (sysPtr->t)*(sysPtr->t);
  CREATE_VIEW(device_,sysPtr->cabana_particles)
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());

  auto kokkos_ensure_consistency = KOKKOS_LAMBDA(const int is, int ia) 
  {
    //Declare caches
    milne::Vector<double,D> u;
    milne::Matrix<double,D+1,D+1> shv;

    for(int idir=0; idir<D; ++idir) u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
    for(int idir=0; idir<D+1; ++idir)
    for(int jdir=0; jdir<D+1; ++jdir)
      shv(idir,jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);

    //Symmetrizes space part
    for( int i=1; i<D+1; i++ )
    for( int j=i+1; j<D+1; j++ )
      shv(j,i) = shv(i,j);

    //pi^{0i} = -\pi^{ij}u_j/gamma = \pi^{ij}u^j/gamma
    milne::Vector<double,D> u_cov = u;
    u_cov.make_covariant(t2); //Transforms u^i to -u_i
    double gamma = Kokkos::sqrt(1+milne::inner(u_cov,u)); //Calculates gamma = \sqrt{1+u^i u_i}
    for(int idir=1; idir<D+1; ++idir){
      milne::Vector<double,D> colp1_shv;
      for(int jdir=1; jdir<D+1; ++jdir) colp1_shv(jdir-1) = shv(idir,jdir);
      shv(0,idir) = 1./gamma*milne::inner(u_cov,colp1_shv); //Computes \pi^{0i} = -\pi^{ij}u_j/gamma
    }

    //Symmetrizes time part
    for( int i=1; i<D+1; i++ )
      shv(i,0) = shv(0,i);

    //Build pi^ij and u^i*u^j
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

    //Transform vector into covariant form, that is u^i*u_j -> u_i*u_j
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
    device_hydro_scalar.access(is, ia, hydro_info::shv33) = shv33;
    //Updates gamma and velocities
    device_hydro_scalar.access(is,ia, hydro_info::gamma) = gamma;
    for(int idir=0; idir<D; ++idir) device_hydro_vector.access(is, ia, hydro_info::v, idir) = u(idir)/gamma;

  };
  Cabana::simd_parallel_for(simd_policy,kokkos_ensure_consistency,"kokkos_ensure_consistency");
  Kokkos::fence();
}

/// @brief Enforces the constraints for the shear viscous tensor \f$ \pi^{\mu\nu} \f$.
/// @details In the (1+1)D scenario, we have 3 independent components of the 
/// shear tensor, \f$ \pi^{00} \f$, \f$ \pi^{01} \f$ and \f$ \pi^{11} \f$. The 
/// constraints gives us fully
/// determined expressions for the components of the shear tensor. By solving
/// the constraints, we get that either \f$ \pi^{\mu\nu} = 0 \f$ or the fluid
/// is moving at the speed of light. We do not foresee the second case to happen
/// in the simulations, so we enforce the absence of shear viscous tensor.
/// @tparam D The number of spatial dimensions.
/// @param sysPtr A pointer to an object of class SystemState.
template<>
void EoM_default<1>::reset_pi_tensor(std::shared_ptr<SystemState<1>> sysPtr)
{
  double t2 = (sysPtr->t)*(sysPtr->t);
  CREATE_VIEW(device_,sysPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());

  auto set_shear_zero = KOKKOS_LAMBDA(const int is, const int ia){
    device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, 0, 0) = 0.0;
    device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, 0, 1) = 0.0;
    device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, 1, 0) = 0.0;
    device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, 1, 1) = 0.0;
    device_hydro_scalar.access(is, ia, hydro_info::shv33) = 0.0;

    milne::Vector<double,1> u;
    for(int idir=0; idir<1; ++idir) u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
    milne::Vector<double,1> u_cov = u;
    u_cov.make_covariant(t2); //Transforms u^i to -u_i
    double gamma = Kokkos::sqrt(1+milne::inner(u_cov,u)); //Calculates gamma = \sqrt{1+u^i u_i}

    //Build pi^ij and u^i*u^j
    device_hydro_space_matrix.access(is, ia, hydro_info::uu, 0, 0) = u(0)*u(0);
    device_hydro_space_matrix.access(is, ia, hydro_info::pimin, 0, 0) = 0;
    ///Updates gamma and velocities
    device_hydro_scalar.access(is,ia, hydro_info::gamma) = gamma;
    device_hydro_vector.access(is, ia, hydro_info::v, 0) = u(0)/gamma;
  };
  Cabana::simd_parallel_for(simd_policy,set_shear_zero,"set_shear_zero");
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

/// @brief Evaluates the time derivatives of the hydrodynamic variables.
/// @details The expressions that is expected to be computed here are for
/// \f$du^i/d\tau\f$, \f$ds/d\tau\f$, \f$d\pi^{ij}/d\tau\f$ and d\Pi/d\tau.
/// This function is broken down in the following steps:
/// - `fill_auxiliary_variables`
/// - `fill_Btot`
/// - `compute_velocity_derivative`
/// - `compute_bulk_derivative`
/// - `compute_shear_derivative`
/// We give an overview of the expressions that are computed in each step.
/// Some of the quantities are computed only if we enable the shear viscosity
/// tensor. These are set to zero if shear viscosity is not enabled. We will
/// annotate these quantities with a "\f$^*\f$".
/// ## `fill_auxiliary_variables`
///
/// The first lambda function fills the cache variables that are used later.
/// The table below shows the quantities that are computed in this step.
/// | Variable    | Expression | Variable | Expression |
/// |-------------|------------|----------|------------|
/// | `dsigma_dt` | \f[\frac{d\sigma}{d\tau} = -\sigma \partial_i v^i \f] | `eta_o_tau`\f$^*\f$ | \f[\frac{\eta}{\tau_{\pi}}\f] |
/// | `g2 = gamma_squared` | \f[\gamma^2 \f]                              | `Agam`      | \f[\begin{align*} A_\gamma = w - \frac{\partial w}{\partial s}\left(s + \frac{\Pi}{T} \right) - \frac{\zeta}{\tau_{\Pi}} - \frac{dw}{d\rho_B}\rho_B - \frac{dw}{d\rho_S}\rho_S - \frac{dw}{d\rho_Q}\rho_Q \end{align*}\f]
/// | `gamma_cube`| \f[\gamma^3 \f]                                       | `sigl`      | \f[ \frac{1}{\tau}\frac{d\sigma}{d\tau} - \frac{1}{\tau} \f] |
/// | `gt`        | \f[\gamma\tau \f]                                     | `Agam2`     | \f[ A_{\gamma,\,2} = \frac{1}{\gamma}\left[A_\gamma - \frac{\eta}{3\tau_\pi} - \pi^{00}\left(1 - \frac{\partial  w}{\partial s}\frac{1}{T}\right)\right] \f]|
/// | `dwdsT1`    | \f[1 - \frac{\partial  w}{\partial s}\frac{1}{T} \f]  | `C`         | \f[w + \Pi \f] |
/// | `bigPI`     | \f[ \Pi = \frac{\tilde \Pi\sigma}{\gamma\tau} \f]     | `Ctot`      | \f[C_{\text{tot}}  = C + \frac{\eta}{\tau_\pi}\left(\frac{1}{\gamma^2}-1\right)\f]  |
/// \f$^*\f$ = Only computed if shear viscosity is enabled.
/// ## `fill_Btot`
///
/// The second lambda function fills the cache variable `Btot`. If shear is 
/// enabled, the expression fill other cache variables as well
/// | Variable         | Expression                       | Variable | Expression |
/// |------------------|----------------------------------|------------------|------------|
/// | `piu`\f$^*\f$    | \f[\pi^{0i}u^j \f]               | `gradU`\f$^*\f$  |  \f[\partial_i u^j = \gamma \partial_i v^j + \gamma^3 v^i v^k \partial_k v^j \f] |
/// | `piutot`\f$^*\f$ | \f[\pi^{0i}u^j + \pi^{0j}u^i \f] | `bsub`\f$^*\f$   |  \f[b = \sum_{i=1}^D \sum_{j=1}^D \partial_i u^j \left(\pi^{ij} + \frac{\pi^{00}}{\gamma^2} u^i u^j - \frac{\pi^{0i}u^j+\pi^{0j}u^i}{\gamma}\right) \f] |
/// | `pig`\f$^*\f$    | \f[\frac{\pi^{00}}{\gamma^2} \f] | `Btot`           |  \f[B = \left(A\gamma + \frac{2}{3}\frac{\eta}{\tau_\pi}\right)\left(\frac{1}{\tau}\frac{d\sigma}{d\tau}-\frac{1}{\tau}\right) + \frac{\Pi}{\tau_\Pi} +\frac{1}{T}\frac{\partial w}{\partial s} (\gamma \tau \pi^{33}+b) \f]  |
/// \f$^*\f$ = Only computed if shear viscosity is enabled.
///
/// ## `compute_velocity_derivative`
///
/// The third lambda function computes the time derivative of the velocity 
/// vector. To this end, it computes a vector $F^i$ and a matrix $M_{ij}$.
/// The derivative is obtained by solving the equation
/// \f[ M^{ij}\frac{du^j}{d\tau} = F^i \f]
///
/// The expressions being computed are
/// | Variable         | Expression | Variable | Expression    |
/// |------------------|---------------------------------------|--------------|------------|
/// | `gamt`\f$^*\f$   | \f[\frac{1}{\gamma \tau_\pi}\f]       | `p1`\f$^*\f$ | \f[p_1 = \frac{1}{\gamma \tau_\pi} - \frac{4}{3 \sigma}\frac{d\sigma}{d\tau} + \frac{1}{3\tau}\f] |
/// | `pre`\f$^*\f$    | \f[\frac{\eta}{\gamma \tau_\pi}\f]    | `M`          | \f[M^{ij} = A_{\gamma,\,2} u^i u^j - \left[1+\frac{4}{3\gamma^2}\pi^{0i}u^j+\left(1-\frac{1}{T}\frac{\partial w}{\partial s}\right)\right]\pi^{0j}u^i + \gamma \pi^{ij} + \delta^{ij} \gamma C_{\text{tot}} \f] |
/// | `partU`\f$^*\f$  | \f[\partial_i u^j + \partial_j u^i\f] | `F`          | \f[F^i = B u^i + v^j\partial_j v^i - \partial_i P - \partial_i \Pi - \partial_j\pi^{ij} + \frac{\eta}{\gamma \tau} v^k (\partial_k u^i + \partial_i u^k) + p_1 \pi^{0i} \f] |
/// | `minshv`\f$^*\f$ | \f[\pi^{0i}\f]                        | | |
/// \f$^*\f$ = Only computed if shear viscosity is enabled.
///
/// ## `compute_bulk_derivative`
///
/// The fourth lambda function computes the time derivative of the bulk pressure.
/// As usual, we also fill a number of cache variables as well. We list the
/// expressions below
/// | Variable          | Expression | Variable | Expression    |
/// |-------------------|-----------------------------------------------------------------|----------  |------------|
/// | `div_u`           | \f[\partial_i u^i  = \frac{1}{\gamma}u_i \frac{du^i}{d\tau} \f] | `bigtheta` | \f [\theta = \tau \partial_i u^i + \gamma \f] |
/// | `d_dt_specific_s` | \f[\frac{ds^*}{d\tau} -\frac{\Pi\Theta}{\sigma T} \f]           | `dBulk_dt` | \f[\frac{d\Pi}{d\tau} = -\frac{\zeta}{\tau_\Pi \sigma \Theta} - \frac{\tilde \Pi}{\gamma\tau_\Pi} \f] |
///
/// ## `compute_shear_derivative`
///
/// The fifth lambda function computes the time derivative of the shear tensor.
/// It is executed only if shear viscosity is enabled. The expressions being
/// computed are
/// | Variable          | Expression | Variable          | Expression |
/// |-------------------|------------|-------------------|------------|
/// | `partU`           | \f[ U^{ij} = \partial_i u^j + \partial_j u^i\f] | `ududt`| \f[ u^i \frac{du^j}{d\tau}\f] |
/// | `vduk`            | \f[ v^k \frac{du^k}{d\tau}\f]                   | `ulpi` | \f[ u^i\pi^{0j}\f] |
/// | `Ipi`             | \f[ I_\pi = -\frac{2}{3} \frac{\eta}{\tau_\pi} (u^i u^j+\delta^{ij}) + \frac{4}{3} \pi^{ij} \f] | | |
///
/// | Variable          | Expression |
/// |-------------------|------------|
/// | `dshv_dt`         | \f[ \frac{d\pi^{ij}}{d\tau} = -\frac{1}{\gamma \tau_\pi} (\pi^{ij} + \eta U^{ij}) - \frac{\eta}{\tau_\pi}\left(u^i \frac{du^j}{d\tau} + u^j \frac{du^i}{d\tau}\right) + I_\pi\left(\frac{1}{\tau}\frac{d\sigma}{d\tau}-\frac{1}{\tau}\right) - v^k \frac{du^k}{d\tau} (u^i \pi^{0j} + u^j \pi^{0i} + I_\pi/\gamma) + (u^i \pi^{jk} + u^j \pi^{ik})\frac{du^k}{d\tau} \f] |
///
/// We also add the contribution of the shear viscosity to the time derivative 
/// of specific entropy. This is given by
///
/// \f[ \frac{ds^*}{d\tau} += \tau (\pi^{00} v^k - \pi^{0k})\frac{du^k}{d\tau} 
/// -\left(\pi^{ij} + \frac{\pi^{00}}{\gamma^2} u^i u^j 
///      - \frac{\pi^{0i}u^j+\pi^{0j}u^i}{\gamma}\right)\partial_i u^j 
/// -\gamma\tau \pi^{33}\f]
///
/// @todo There are several quantities that are considered particle properties
/// but I believe are just cache variables used during the computation. These
/// should be moved to an AoSoA declared in the EoM class and filled here. This
/// will make the code more readable, easier to maintain and will allow for 
/// faster data movement between device and host.
/// @todo These equations of motion are specific for the transverse plane. We
/// need to review this to make sure that the equations are correct for the
/// longitudinal direction as well. A step in this direction is to make sure
/// that the equations are written in a covariant way and then implement them
/// in the code. Later, we can look for missing terms (Christoffel symbols) in 
/// the equations.
/// @todo We also need to include the expresions for the Beta_Bulk derivative.
/// @tparam D The number of spatial dimensions.
/// @param sysPtr A pointer to the object of class SystemState.
template<unsigned int D>
void EoM_default<D>::evaluate_time_derivatives( std::shared_ptr<SystemState<D>> sysPtr, std::shared_ptr<Settings> settingsPtr)
{
  // #ifdef DEBUG
  // ofstream outfile;
  // outfile.open("gradients.dat");
  // #endif
  double t = (sysPtr->t);
  double t2 = t*t;
  CREATE_VIEW(device_,sysPtr->cabana_particles);
  bool using_shear = settingsPtr->using_shear; 
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, sysPtr->cabana_particles.size());

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

    
    double dsigma_dt = -sigma * milne::tr(gradV,t2);
    double g2        = gamma*gamma;
    double gt        = gamma*t;
    double dwdsT1    = 1 - dwds/T;
    double bigPI     = Bulk*sigma/gt;
    double C         = w + bigPI;

    double eta_o_tau = using_shear ? setas/stauRelax : 0.0;

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

  Cabana::simd_parallel_for(simd_policy, fill_auxiliary_variables, "fill_auxiliary_variables");
  Kokkos::fence();

  auto fill_Btot = KOKKOS_LAMBDA(const int is, const int ia){
    double bsub = 0;
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

      piu    = milne::rowp1(0,shv)*u;
      piutot = piu + milne::transpose(piu);

      double pig  = shv(0,0)/g2;

      gradU        = gamma*gradV+g3*(v*(v*gradV));
      for (int i=0; i<D; i++)
      for (int j=0; j<D; j++)
        bsub = bsub + gradU(i,j) * ( pimin(i,j) + pig*uu(j,i) - piutot(i,j)/gamma );
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

  Cabana::simd_parallel_for(simd_policy, fill_Btot, "fill_Btot");
  Kokkos::fence();

  auto compute_velocity_derivative = KOKKOS_LAMBDA(const int is, const int ia)
  {
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
        milne::Matrix <double,D,D> partU = gradU+milne::transpose(gradU);
        milne::Vector<double,D> minshv = milne::rowp1(0, shv);
        F += pre*v*partU + p1*minshv;
      }

      milne::Matrix<double,D,D> MI = milne::inverse(M);
      milne::Vector<double,D> du_dt = MI*F;
      for(int idir=0; idir<D; ++idir) device_hydro_vector.access(is,ia,hydro_info::du_dt, idir) = du_dt(idir);

      #ifdef DEBUG
      double pos = device_position.access(is, ia, 0);
      if (pos < 7.52 && pos > 7.48) {
	std::cout << t << " " << pos << " " << \
        device_hydro_scalar.access(is, ia, hydro_info::gamma)*device_hydro_space_matrix.access(is, ia, hydro_info::gradV, 0, 0) \
        << " " << device_hydro_vector.access(is, ia, hydro_info::u, 0) << " " <<  gradP(0) << std::endl;
      }
      if (t >= 1.2){
        exit(1);
      }
      #endif
  };
  Cabana::simd_parallel_for(simd_policy, compute_velocity_derivative, "compute_velocity_derivative");
  Kokkos::fence();

  auto compute_bulk_derivative = KOKKOS_LAMBDA(const int is, const int ia)
  {

    double eta_o_tau = device_hydro_scalar.access(is, ia, hydro_info::eta_o_tau);
    double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
    double dsigma_dt = device_hydro_scalar.access(is, ia, hydro_info::dsigma_dt);
    double g2 = device_hydro_scalar.access(is, ia, hydro_info::gamma_squared);
    double T = device_thermo.access(is, ia, thermo_info::T);
    double bigPI = device_hydro_scalar.access(is, ia, hydro_info::bigPI);
    double zeta = device_hydro_scalar.access(is, ia, hydro_info::zeta);
    double Bulk = device_hydro_scalar.access(is, ia, hydro_info::Bulk);
    double tauRelax = device_hydro_scalar.access(is, ia, hydro_info::tauRelax);
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);

    milne::Vector<double,D> v, du_dt, u;
    for(int idir=0; idir<D; ++idir){
      v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
      du_dt(idir) = device_hydro_vector.access(is, ia, hydro_info::du_dt, idir);
      u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
    }
    milne::Matrix<double,D+1,D+1> shv;
    for(int idir=0; idir<D+1; ++idir)
    for(int jdir=0; jdir<D+1; ++jdir)
      shv(idir, jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);

    milne::Matrix<double,D,D> uu, pimin, piutot;
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir)
    {
      uu(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::uu, idir, jdir);
      pimin(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::pimin, idir, jdir);
      piutot(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::piutot, idir, jdir);
    }

    //===============
    // "coordinate" divergence
    double div_u = (1./gamma)*inner(u, du_dt) - ( gamma/ sigma ) * dsigma_dt;
    //===============
    // "covariant" divergence
    double bigtheta = div_u*t+gamma;

    // time derivative of ``specific entropy density per particle"
    double d_dt_specific_s = 1./sigma/T*(-bigPI*bigtheta);
	  //Bulk evolution equation
    double dBulk_dt = ( -zeta/sigma*bigtheta - Bulk/gamma )/tauRelax;//w/out hi.dBeta_dt

    device_hydro_scalar.access(is, ia, hydro_info::div_u) = div_u;
    device_hydro_scalar.access(is, ia, hydro_info::bigtheta) = bigtheta;
    device_hydro_scalar.access(is, ia, hydro_info::dBulk_dt) = dBulk_dt;
    device_d_dt_spec.access(is, ia, densities_info::s) = d_dt_specific_s;

	  //formulating simple setup for Beta_Bulk derivative
    ///TODO: Implement this
    //hi.finite_diff_cs2   =  (ti.cs2 - hi.prev_cs2)/0.05; // Asadek
	  //hi.finite_diff_T   =  (ti.T - hi.prev_T)/0.05; // Asadek
	  //hi.finite_diff_w   =  (ti.w - hi.prev_w)/0.05; // Asadek
	  //hi.dBeta_dt      = 0.5*((-hi.finite_diff_T/(ti.T*ti.T))*(1/ti.w)*(1/((1/3-ti.cs2)*(1/3-ti.cs2))))
	  //                 + 0.5*((-hi.finite_diff_w/(ti.w*ti.w))*(1/ti.T)*(1/((1/3-ti.cs2)*(1/3-ti.cs2))))
		//			   + 0.5*((4*ti.cs2*hi.finite_diff_cs2*(1/((1/3-ti.cs2)*(1/3-ti.cs2)*(1/3-ti.cs2))))*(1/ti.T)*(1/ti.w));//Asadek 

  };
  Cabana::simd_parallel_for(simd_policy, compute_bulk_derivative, "compute_bulk_derivative");
  Kokkos::fence();

  if (using_shear){
    auto compute_shear_derivative = KOKKOS_LAMBDA(const int is, const int ia)
    {
      double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
      double stauRelax = device_hydro_scalar.access(is, ia, hydro_info::stauRelax);
      double setas = device_hydro_scalar.access(is, ia, hydro_info::setas);
      double eta_o_tau = device_hydro_scalar.access(is, ia, hydro_info::eta_o_tau);
      double sigl = device_hydro_scalar.access(is, ia, hydro_info::sigl);
      double g2 = device_hydro_scalar.access(is, ia, hydro_info::gamma_squared);
      double shv33 = device_hydro_scalar.access(is, ia, hydro_info::shv33);
      double sigma = device_hydro_scalar.access(is, ia, hydro_info::sigma);
      double T = device_thermo.access(is, ia, thermo_info::T);

      milne::Matrix<double,D,D> pimin, gradU, uu, piutot;
      for(int idir=0; idir<D; ++idir)
      for(int jdir=0; jdir<D; ++jdir){
        pimin(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::pimin, idir, jdir);
        uu(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::uu, idir, jdir);
        gradU(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::gradU, idir, jdir);
        piutot(idir,jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::piutot, idir, jdir);
      }

      milne::Vector<double,D> v, du_dt, u;
      for(int idir=0; idir<D; ++idir){
        du_dt(idir) = device_hydro_vector.access(is, ia, hydro_info::du_dt, idir);
        u(idir) = device_hydro_vector.access(is, ia, hydro_info::u, idir);
        v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
      }

      milne::Matrix<double,D+1,D+1> shv;
      for(int idir=0; idir<D+1; ++idir)
      for(int jdir=0; jdir<D+1; ++jdir)
        shv(idir, jdir) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);


      milne::Matrix<double,D,D> partU = gradU+milne::transpose(gradU);
      milne::Matrix<double,D,D> ududt = u*du_dt;
      double gamt = 1.0/gamma/stauRelax;
      double vduk = milne::inner(v, du_dt);
      milne::Matrix <double,D,D> ulpi  = u*milne::colp1(0, shv);
      milne::Matrix <double,D,D> Ipi   = - 2.0*eta_o_tau/3.0 * uu + 4./3.*pimin;
      for(int idir=0; idir<D; ++idir) Ipi(idir, idir) -= 2.0*eta_o_tau/3.0;

      milne::Matrix <double,D,D> dpidtsub;

      for (int i=0; i<D; i++)
      for (int j=0; j<D; j++)
      for (int k=0; k<D; k++)
        dpidtsub(i,j) += ( u(i)*pimin(j,k) + u(j)*pimin(i,k) )*du_dt(k);


      milne::Matrix <double,D,D> dshv_dt = - gamt*( pimin + setas*partU ) 
                                 - eta_o_tau*( ududt + milne::transpose(ududt) )
                                 + dpidtsub + sigl*Ipi
                                 - vduk*( ulpi + milne::transpose(ulpi) + (1/gamma)*Ipi );

      milne::Matrix<double,D,D> sub = pimin + (shv(0,0)/g2)*uu
                                  - 1./gamma*piutot;
      milne::Vector<double,D> minshv = milne::rowp1(0, shv);
      double inside = t*( milne::inner( shv(0,0)*v-minshv, du_dt )
                                      - milne::con2(sub, gradU) - gamma*t*shv33 );

      // time derivative of ``specific entropy density per particle"
      double d_dt_specific_s = 1./sigma/T*inside;
      device_hydro_scalar.access(is, ia, hydro_info::inside) = inside;
      device_d_dt_spec.access(is, ia, densities_info::s) += d_dt_specific_s;
      for(int idir=0; idir<D; ++idir)
      for(int jdir=0; jdir<D; ++jdir){
        device_hydro_space_matrix.access(is, ia, hydro_info::dshv_dt, idir, jdir) = dshv_dt(idir, jdir);

      }

    };
    Cabana::simd_parallel_for(simd_policy, compute_shear_derivative, "compute_shear_derivative");
    Kokkos::fence();

  }
}
