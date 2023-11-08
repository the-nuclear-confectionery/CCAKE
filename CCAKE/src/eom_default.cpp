#include "eom_default.h"

namespace ccake{
  template class EoM_default<1>;
  template class EoM_default<2>;
  template class EoM_default<3>;

template<unsigned int D>
void EoM_default<D>::reset_pi_tensor(std::shared_ptr<SystemState<D>> sysPtr)
{
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
double EoM_default<D>::gamma_calc(double u[D], const double &time_squared)
{
    double dot_u = 0;
    for (unsigned int i=0; i<D-1; i++)
      dot_u+= u[i]*u[i];
    dot_u += time_squared*u[D-1]*u[D-1];
    return sqrt(1.0+dot_u);
}

/// @brief Calculates the gamma factor = u^0.
/// @details This is the default implementation of the gamma factor calculation.
/// It assumes that the last component of the velocity vector is the longitudinal
/// velocity, and that the metric is Milne.
/// @param u The velocity vector (space components only).
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return the value of gamma
template<> KOKKOS_FUNCTION
double EoM_default<2>::gamma_calc(double u[2], const double &time_squared)
{
    double dot_u = 0;
    for (unsigned int i=0; i<2; i++)
      dot_u += u[i]*u[i];
    return sqrt(1.0+dot_u);
}

/// @brief Transforms a scalar from the lab frame to the LRF.
/// @tparam D
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
      for (int i=0; i<=1; i++)
      for (int j=0; j<=1; j++)
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
        milne::Matrix <double,D,D> partU;
        for (int i=0; i<=1; i++)
        for (int j=0; j<=1; j++)
          partU(i,j) = gradU(i,j) + gradU(j,i);
        milne::Vector<double,D> minshv = milne::rowp1(0, shv);
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
      cout << h_matrix_view(c,hydro_info::piutot, i, j) << " ";
    cout << endl;
  }
  cout << "pimin......: " << endl;
  for(int j=0;  j<D; ++j){
    for(int i=0;  i<D; ++i)
      cout << h_matrix_view(c,hydro_info::pimin, i, j) << " ";
    cout << endl;
  }
  #endif
}
}