#include "sph_workstation.h"
//The cpp below needs to be included for the Kokkos lambda to be
//properly compiled in the GPU
#include "kernel.cpp"
#include "eom_default.cpp"
#include "transport_coefficients.cpp"

using namespace constants;
namespace tc = ccake::transport_coefficients;

namespace ccake{
//Template instantiations
template class SPHWorkstation<1,EoM_default>;
template class SPHWorkstation<2,EoM_default>;
template class SPHWorkstation<3,EoM_default>;

/// @brief Initialize the several components of the SPHWorkstation.
/// @details This function will initialize the several components of the
// SPHWorkstation. It will create and initialize the equation of motion
/// object, initialize the system state, create and initialize transport
/// coefficients, create and initialize the freeze out object, and create
/// and initialize the Runge-Kutta evolver.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::initialize()
{
  //----------------------------------------
  // set up system state
  systemPtr->initialize();

  //----------------------------------------
  // set up equation of state
  eos.set_SettingsPtr( settingsPtr );
  eos.init();

  //----------------------------------------
  // set up transport coefficients
  transp_coeff_params = tc::setup_parameters(settingsPtr);

  //----------------------------------------
  // set up freeze out (constant energy density)

  systemPtr->efcheck = eos.efreeze(settingsPtr->Freeze_Out_Temperature);
  systemPtr->sfcheck = eos.sfreeze(settingsPtr->Freeze_Out_Temperature);
  formatted_output::detail("freeze out energy density = "
                           + to_string(systemPtr->efcheck*hbarc_GeVfm) + " GeV/fm^3");
  formatted_output::detail("freeze out entropy density = "
                           + to_string(systemPtr->sfcheck) + " 1/fm^3");

  freeze_out.initialize( settingsPtr, systemPtr, systemPtr->efcheck );

  //----------------------------------------
  // set up RK evolver
  //evolver.initialize( settingsPtr, systemPtr );

}


/// @brief This function will ensure that the shear tensor has the correct properties.
/// @details This function will ensure that the shear tensor has the correct properties.
/// The desired properties are that it is traceless, symmetric and orthogonal to the
/// fluid velocity. From $u_\mu \pi^{\mu i}$ one can compute the components $\pi^{0i}$.
/// of the tensor. From $u_\mu \pi^{\mu 0}$ one can compute the component $\pi^{00}$.
/// Lastly, the component $pi^{33}$ can be computed from the tracelessness condition.
/// The remaining components are computed from the symmetry condition.
/// @tparam D
/// @param time_squared The square of the time step whwere the shear tensor is being computed.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::reset_pi_tensor(double time_squared)
{
    CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto reset_pi_tensor_lambda  = KOKKOS_LAMBDA(const int is, const int ia)
  {
    double u[D]; // Cache for fluid velocity
    double pi_space[D][D]; //Cache for pi^{ij}
    double pi_time_vector[D]; //Cache for the components pi^{0i}
    double pi00; //Cache for pi^{00}
    double pi_diag[D+1]; //Cache for the diagonal components of pi tensor

    //Retrieve information from the particle to local variables
    for( unsigned int idir=0; idir<D; idir++ ){
      u[idir] = device_hydro_vector.access(is, ia, ccake::hydro_info::u, idir);
      pi_diag[idir+1] = device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, idir+1, idir+1);
      for( unsigned int jdir=0; jdir<D; jdir++ )
        pi_space[idir][jdir] = device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, idir+1, jdir+1);
    }

    //Computes the gamma factor
    double u0 = TEOM<D>::gamma_calc(u,time_squared);

    //Computes pi^{0i} = - u_j pi^{ij}/gamma (obtained from the requirement u_\mu pi^{\mu i} = 0)
    for( unsigned int idir=0; idir<D; idir++ )
      pi_time_vector[idir] = TEOM<D>::dot(pi_space[idir],u,time_squared)/u0;

    //Set the shear tensor pi^{0 0} = - pi^{0j} u_j /gamma (from the requirement u_\mu pi^{\mu 0} = 0)
    //This is equivalent to pi^{i j} u^i u^j /gamma^2
    pi00 = TEOM<D>::dot(pi_time_vector,u,time_squared)/u0;
    pi_diag[0] = pi00;

    //Lastly, we need to ensure that the tensor is traceless. We do this by setting
    //the value of the last component pi^{33} = (pi^{00}-pi^{11}-pi^{22})/tau^22
    double pi33 = TEOM<D>::get_shvDD(pi_diag, time_squared);
    if( D != 2 )
      pi_diag[D] = pi33;


    //Set the values of the pi tensor
    device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, 0, 0) = pi_diag[0]; //pi^{00}
    for( unsigned int idir=1; idir<D+1; idir++ ){
      device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, 0, idir) = pi_time_vector[idir-1]; //pi^{0i}
      device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, idir, 0) = pi_time_vector[idir-1]; //pi^{i0}
      device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, idir, idir) = pi_diag[idir]; //pi^{ii}
      for (unsigned int jdir=idir+1; jdir<D+1; jdir++){
        device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, idir, jdir) = pi_space[idir-1][jdir-1]; //pi^{ij}
        device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, jdir, idir) = pi_space[idir-1][jdir-1]; //pi^{ji}
      }
    }
    //Set values of pi^{ij} and u^i u^j tensors
    for( unsigned int idir=0; idir<D; idir++ )
    for( unsigned int jdir=0; jdir<D; jdir++ )
    {
      device_hydro_space_matrix.access(is, ia, ccake::hydro_info::pimin, idir, jdir) = pi_space[idir][jdir]; //pi^{ij}
      device_hydro_space_matrix.access(is, ia, ccake::hydro_info::uu, jdir, idir) = pi_space[idir][jdir]; //pi^{ji}
    }
    //Update gamma factor
    device_hydro_scalar.access(is, ia, ccake::hydro_info::gamma) = u0;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::shv33) = pi33;
  };

  Cabana::simd_parallel_for(*(systemPtr->simd_policy), reset_pi_tensor_lambda, "reset_pi_tensor");
  Kokkos::fence();

};


////////////////////////////////////////////////////////////////////////////////
///\brief Set up entropy density and other thermodynamic variables
///\details This function sets up the entropy density using the initial energy
///         density and the initial charge densities. It also sets up all other
///         thermodynamic variables.
///
///         The norm_spec densities are also updated here.
///\todo This function is not yet implemented for realistic EoS.

template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::initialize_entropy_and_charge_densities()
{
	Stopwatch sw, swTotal;
	swTotal.Start();
	long long failCounter = 0;

  if ( VERBOSE > 5 )
  {
    cout << "----------------------------------------"
            "----------------------------------------" << "\n";
    cout << "----------------------------------------"
            "----------------------------------------" << endl;
  }

  systemPtr->n_particles = systemPtr->particles.size();

  //Compute thermal properties
  systemPtr->copy_device_to_host();
  ///\TODO: Trivially parallelizable via openmp. This should be done.
  for (auto & p : systemPtr->particles){
    p.input.s = locate_phase_diagram_point_eBSQ( p,
                  p.input.e, p.input.rhoB, p.input.rhoS, p.input.rhoQ );
  }
  systemPtr->copy_host_to_device();

  double t0 = settingsPtr->t0;
  // loop over all particles
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto init_particles = KOKKOS_LAMBDA(const int iparticle)
  {
    double e_input = device_input(iparticle, densities_info::e);
    double s_input = device_input(iparticle, densities_info::s);
    double rhoB_input = device_input(iparticle, densities_info::rhoB);
    double rhoS_input = device_input(iparticle, densities_info::rhoS);
    double rhoQ_input = device_input(iparticle, densities_info::rhoQ);

    double gamma = device_hydro_scalar(iparticle, ccake::hydro_info::gamma);

    device_norm_spec(iparticle, ccake::densities_info::e )    *= e_input*gamma*t0;    // constant after this
    device_norm_spec(iparticle, ccake::densities_info::s )    *= s_input*gamma*t0;       // constant after this
    device_norm_spec(iparticle, ccake::densities_info::rhoB ) *= rhoB_input*gamma*t0; // constant after this
    device_norm_spec(iparticle, ccake::densities_info::rhoS ) *= rhoS_input*gamma*t0; // constant after this
    device_norm_spec(iparticle, ccake::densities_info::rhoQ ) *= rhoQ_input*gamma*t0; // constant after this

    if (s_input < 0.0)
      device_freeze(iparticle) = 4;

    ///------------------------------------------------------------------------------------------
  };

  Kokkos::parallel_for("init_particles", systemPtr->n_particles, init_particles);

  //Allocate cache for evolver
  evolver.allocate_cache();


	swTotal.Stop();
  formatted_output::update("finished initializing particle densities in "
                              + to_string(swTotal.printTime()) + " s");

}

//==============================================================================
template<unsigned int D,  template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::initial_smoothing()
{
  Stopwatch sw;
  sw.Start();
  double t_squared = pow(settingsPtr->t0,2);

  //Creates the list of neighbours for each particle
  systemPtr->reset_neighbour_list();

  //Computes not independent components of the shear tensor
  reset_pi_tensor(t_squared);

  // smooth fields over particles
  smooth_all_particle_fields(t_squared);
  sw.Stop();
  formatted_output::update("Finished initial smoothing "
                            + to_string(sw.printTime()) + " s.");

}

//==============================================================================

///@brief Smooth all SPH fields
///@details This function updates the densities s, rhoB, rhoB, rhoS and sigma
///(auxiliary density) by performing the smoothing procedure. At its end, it
///takes the opportunity to update the remaining thermodynamic quantities
///(energy density, pressure, chemical potentials, temperature and speed of
/// sound squared). The routine that triggers the update of thermo quantities
/// will also update gamma and the non-relativistic velocity
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::smooth_all_particle_fields(double time_squared)
{

  CREATE_VIEW(device_, systemPtr->cabana_particles);
  double hT = settingsPtr->hT;

  //Reset smoothed fields
  auto reset_fields = KOKKOS_LAMBDA(const int is, const int ia) //First index for loop over struct, second for loop over array
  {
    device_smoothed.access(is, ia, ccake::densities_info::s) = 0.0;
    device_smoothed.access(is, ia, ccake::densities_info::rhoB) = 0.0;
    device_smoothed.access(is, ia, ccake::densities_info::rhoS) = 0.0;
    device_smoothed.access(is, ia, ccake::densities_info::rhoQ) = 0.0;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::sigma) = 0.0;
  };
  Cabana::simd_parallel_for( *(systemPtr->simd_policy), reset_fields, "reset_fields" );
  Kokkos::fence();

  auto smooth_fields = KOKKOS_CLASS_LAMBDA(const int iparticle, const int jparticle ){

    double r1[D] ,r2[D]; // cache for positions of particles 1 and 2
    for (int idir = 0; idir < D; ++idir){
      r1[idir] = device_position(iparticle,idir);
      r2[idir] = device_position(jparticle,idir);
    }
    double distance = SPHkernel<D>::distance(r1,r2);
    double kern = SPHkernel<D>::kernel(distance,hT);

    //Update sigma (reference density)
    Kokkos::atomic_add( &device_hydro_scalar(iparticle, ccake::hydro_info::sigma), device_norm_spec(jparticle, ccake::densities_info::s)*kern);
    ////Update entropy density
    Kokkos::atomic_add( &device_smoothed(iparticle, ccake::densities_info::s),
                         device_norm_spec(jparticle, ccake::densities_info::s)*device_specific_density(jparticle, ccake::densities_info::s)*kern);
    Kokkos::atomic_add( &device_smoothed(iparticle, ccake::densities_info::rhoB),
                          device_norm_spec(jparticle, ccake::densities_info::rhoB)*device_specific_density(jparticle, ccake::densities_info::rhoB)*kern);
    Kokkos::atomic_add( &device_smoothed(iparticle, ccake::densities_info::rhoS),
                          device_norm_spec(jparticle, ccake::densities_info::rhoS)*device_specific_density(jparticle, ccake::densities_info::rhoS)*kern);
    Kokkos::atomic_add( &device_smoothed(iparticle, ccake::densities_info::rhoQ),
                          device_norm_spec(jparticle, ccake::densities_info::rhoQ)*device_specific_density(jparticle, ccake::densities_info::rhoQ)*kern);
  };

  Cabana::neighbor_parallel_for( systemPtr->range_policy, smooth_fields,
                                 systemPtr->neighbour_list, Cabana::FirstNeighborsTag(),
                                 Cabana::TeamOpTag(), "smooth_fields_kernel");
  Kokkos::fence();

  update_all_particle_thermodynamics(time_squared);

}

//------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::locate_phase_diagram_point_sBSQ( Particle<D> & p,
                 double s_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{

  // default: use particle's current location as initial guess
  eos.tbqs( p.T(), p.muB(), p.muQ(), p.muS(), p.get_current_eos_name() );

  bool update_s_success = eos.update_s( s_In, rhoB_In, rhoS_In, rhoQ_In,
                                        p.print_this_particle );

  if ( update_s_success )
  {
    // check that enthalpy derivatives are not crazy
    thermodynamic_info p0 = p.thermo;

    eos.set_thermo( p.thermo );

    // do not permit cs2 to be negative if using C library
    // or if using the tanh_conformal EoS
    if (  ( eos.currently_using_static_C_library()
            && p.get_current_eos_name() == "table" )
          || p.get_current_eos_name() == "tanh_conformal" )
    {
      if ( p.thermo.cs2 < 0.0 )
      {
        cout << "WARNING: cs2 went negative in eos_type == "
            << p.get_current_eos_name() << " for particle " << p.ID << "\n"
            << "input thermo: " << s_In << "   "
            << rhoB_In << "   "
            << rhoS_In << "   "
            << rhoQ_In << endl
            << "check thermo: " << systemPtr->t << "   "
            << p.thermo.T << "   "
            << p.thermo.muB << "   "
            << p.thermo.muS << "   "
            << p.thermo.muQ << "   "
            << p.thermo.p << "   "
            << p.thermo.s << "   "
            << p.thermo.rhoB << "   "
            << p.thermo.rhoS << "   "
            << p.thermo.rhoQ << "   "
            << p.thermo.e << "   "
            << p.thermo.cs2 << "   "
            << p.thermo.eos_name << endl;
        p.thermo.cs2 = std::max( p.thermo.cs2, 0.0001 );
      }
    }

    if ( p0.eos_name == p.get_current_eos_name() )
    {
      auto abslogabs = [](double x){ return std::abs(std::log(std::abs(x))); };
      bool dwds_possibly_unstable = abslogabs((p.dwds()+1e-6)/(p0.dwds+1e-6)) > 1e6;
      bool dwdB_possibly_unstable = abslogabs((p.dwdB()+1e-6)/(p0.dwdB+1e-6)) > 1e6;
      bool dwdS_possibly_unstable = abslogabs((p.dwdS()+1e-6)/(p0.dwdS+1e-6)) > 1e6;
      bool dwdQ_possibly_unstable = abslogabs((p.dwdQ()+1e-6)/(p0.dwdQ+1e-6)) > 1e6;
      if ( dwds_possibly_unstable || dwdB_possibly_unstable
            || dwdS_possibly_unstable || dwdQ_possibly_unstable )
        std::cerr << "WARNING: thermodynamics may be unstable!\n"
                  << "\t --> particle: " << p.ID << "   " << systemPtr->t << "\n"
                  << "\t --> (before)  " << p0.T << "   " << p0.muB << "   "
                  << p0.muS << "   " << p0.muQ << "\n"
                  << "\t               " << p0.s << "   " << p0.rhoB << "   "
                  << p0.rhoS << "   " << p0.rhoQ << "\n"
                  << "\t               " << p0.dwds << "   " << p0.dwdB << "   "
                  << p0.dwdS << "   " << p0.dwdQ << "\n"
                  << "\t               " << p0.p << "   " << p0.e << "   "
                  << p0.cs2 << "   " << "\n"
                  << "\t --> (after)   " << p.T() << "   " << p.muB() << "   "
                  << p.muS() << "   " << p.muQ() << "\n"
                  << "\t               " << p.s() << "   " << p.rhoB() << "   "
                  << p.rhoS() << "   " << p.rhoQ() << "\n"
                  << "\t               "
                  << p.dwds() << "   " << p.dwdB() << "   "
                  << p.dwdS() << "   " << p.dwdQ() << "\n"
                  << "\t               " << p.p() << "   " << p.e() << "   "
                  << p.cs2() << "   " << std::endl;
    }


    if (p.thermo.cs2<0)
    {
    cout << "input thermo: " << s_In << "   "
          << rhoB_In << "   "
          << rhoS_In << "   "
          << rhoQ_In << endl
          << "check thermo: " << systemPtr->t << "   "
          << p.thermo.T << "   "
          << p.thermo.muB << "   "
          << p.thermo.muS << "   "
          << p.thermo.muQ << "   "
          << p.thermo.p << "   "
          << p.thermo.s << "   "
          << p.thermo.rhoB << "   "
          << p.thermo.rhoS << "   "
          << p.thermo.rhoQ << "   "
          << p.thermo.e << "   "
          << p.thermo.cs2 << "   "
          << p.thermo.eos_name << endl;
      cout << __LINE__ << ": cs2 was negative!" << endl;
      exit(8);
    }
  }

  return;
}



///////////////////////////////////////////////////////////////////////////////
// smoothing routines: first smoothing covers all hydrodyanmical fields


///////////////////////////////////////////////////////////////////////////////
//Second smoothing smoothes the gradients after constructing all the fields
//and derivatives using the equation of state
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::smooth_all_particle_gradients(double time_squared)
{
  double hT = settingsPtr->hT;
  CREATE_VIEW(device_, systemPtr->cabana_particles);

  //Reset gradients
  auto reset_gradients = KOKKOS_LAMBDA(const int iparticle){
    for (int idir=0; idir < D; ++idir)
    {
      device_hydro_vector(iparticle,ccake::hydro_info::gradP, idir) = 0.0;
      device_hydro_vector(iparticle,ccake::hydro_info::gradBulk, idir) = 0.0;
      device_hydro_vector(iparticle,ccake::hydro_info::gradshear, idir) = 0.0;
      device_hydro_vector(iparticle,ccake::hydro_info::divshear, idir) = 0.0;
      for (int jdir=0; jdir < D; ++jdir)
        device_hydro_space_matrix(iparticle,ccake::hydro_info::gradV, idir, jdir) = 0.0;
    }

    int btrack = device_btrack(iparticle);
    if ( btrack != -1 ) btrack = 0;
    device_btrack(iparticle) = btrack;
  };
  Kokkos::parallel_for("reset_gradients", systemPtr->n_particles, reset_gradients);
  Kokkos::fence();

  ///particle_a is the one that will be updated. Particle b is its neighbors.
  bool using_shear = settingsPtr->using_shear;
  auto smooth_gradients = KOKKOS_LAMBDA(const int particle_a, const int particle_b )
  {
    //Cache quantities locally
    double sigma_a        = device_hydro_scalar(particle_a, ccake::hydro_info::sigma);
    double sigma_b        = device_hydro_scalar(particle_b, ccake::hydro_info::sigma);
    double entropy_norm_b = device_norm_spec(particle_b, ccake::densities_info::s); //"Mass of particle b"
    double sigsqra        = 1.0/(sigma_a*sigma_a);
    double sigsqrb        = 1.0/(sigma_b*sigma_b);
    double pressure_a     = device_thermo(particle_a, ccake::thermo_info::p);
    double pressure_b     = device_thermo(particle_b, ccake::thermo_info::p);
    double energy_a       = device_thermo(particle_a, ccake::densities_info::e);
    double energy_b       = device_thermo(particle_b, ccake::densities_info::e);
    double bulk_a         = device_hydro_scalar(particle_a, ccake::hydro_info::Bulk);
    double bulk_b         = device_hydro_scalar(particle_b, ccake::hydro_info::Bulk);
    double gamma_a        = device_hydro_scalar(particle_a, ccake::hydro_info::gamma);
    double gamma_b        = device_hydro_scalar(particle_b, ccake::hydro_info::gamma);
    int pa_btrack         = device_btrack(particle_a);
    int pb_btrack         = device_btrack(particle_b);

    //Auxiliary local variables
    double rel_sep[D]; // cache for relative separation
    double pos_a[D]; // cache for position of particle a
    double pos_b[D]; // cache for position of particle b
    double vel_a[D], vel_b[D]; // cache for velocity of particles a and b
    double ua[D], ub[D]; // cache for fluid four-velocity of particles a and b
    double va[D], vb[D];
    double sigsigK[D];
    double vminia[D][D], vminib[D][D];
    //double pa_qmom[D], pb_qmom[D], qmom_difference[D]; //<Only for pressure weight

    //Initialize vectors and matrices for cache and auxiliary variables
    for (int idir = 0; idir < D; ++idir){
      pos_a[idir] = device_position(particle_a,idir);
      pos_b[idir] = device_position(particle_b,idir);
      rel_sep[idir] = pos_a[idir] - pos_b[idir];
      va[idir] = device_hydro_spacetime_matrix(particle_a, ccake::hydro_info::shv, 0, idir+1); //pi^{0i}
      vb[idir] = device_hydro_spacetime_matrix(particle_b, ccake::hydro_info::shv, 0, idir+1);
      vel_a[idir] = device_hydro_vector(particle_a, ccake::hydro_info::v, idir);
      vel_b[idir] = device_hydro_vector(particle_b, ccake::hydro_info::v, idir);
      ua[idir] = device_hydro_vector(particle_a, ccake::hydro_info::u, idir);
      ub[idir] = device_hydro_vector(particle_b, ccake::hydro_info::u, idir);
      //pa_qmom[idir] = ((energy_a+pressure_a)*gamma_a/sigma_a)*ua[idir]; //< Only for pressure weight
      //pb_qmom[idir] = ((energy_b+pressure_b)*gamma_b/sigma_b)*ub[idir]; //< Only for pressure weight
      //qmom_difference[idir] = pa_qmom[idir] - pb_qmom[idir]; //< Only for pressure weight
      for (int jdir=0; jdir<D; jdir++){
        vminia[idir][jdir] = device_hydro_spacetime_matrix(particle_a, ccake::hydro_info::shv, idir+1, jdir+1);
        vminib[idir][jdir] = device_hydro_spacetime_matrix(particle_b, ccake::hydro_info::shv, idir+1, jdir+1);
      }
    }

    //Compute kernel gradients
    double rel_sep_norm = SPHkernel<D>::distance(pos_a,pos_b);
    double gradK[D];
    SPHkernel<D>::gradKernel( rel_sep, rel_sep_norm, hT, gradK );

    //Get gradients of vectors
    for (int idir=0; idir<D; idir++){
      sigsigK[idir] = entropy_norm_b*sigma_a * gradK[idir];
      double gradP = (sigsqrb*pressure_b + sigsqra*pressure_a)*sigsigK[idir];
      double gradBulk = (bulk_b/sigma_b/gamma_b+bulk_a/sigma_a/gamma_a)/sqrt(time_squared)*sigsigK[idir]; ///wmatioli: Why not use the same approach as pressure?
      for (int jdir=0; jdir<D; jdir++)
      {
        double gradV = (entropy_norm_b/sigma_a)*( vel_b[idir] - vel_a[idir] )*gradK[jdir]; ///wmatioli: Related to dsigma^*/dx^0?
        Kokkos::atomic_add( &device_hydro_space_matrix(particle_a, ccake::hydro_info::gradV, idir, jdir), gradV);
      }
      Kokkos::atomic_add( &device_hydro_vector(particle_a, ccake::hydro_info::gradP, idir), gradP);
      Kokkos::atomic_add( &device_hydro_vector(particle_a, ccake::hydro_info::gradBulk, idir), gradBulk);
    }

    ///Counts nearest neighbors btrack
    ///TODO: Cabana has methods to do this. We should use them instead of doing this.
    double relative_distance_by_h = rel_sep_norm / hT;
    if ( ( relative_distance_by_h <= 2.0 ) && ( particle_a != particle_b ) && ( pa_btrack != -1 ) )
      Kokkos::atomic_add( &device_btrack(particle_a), 1 ); // effectively counts nearest neighbors

    //===============
    // add shear terms
    if ( using_shear )
    {
      double vel_dot_grad = TEOM<D>::dot(vel_a,sigsigK,time_squared);
      double aux_a[D];
      double aux_b[D];
      for(int idir=0; idir<D; ++idir){
        aux_a[idir]=0; aux_b[idir]=0;
      }
      ///TODO: This is right for (2+1)D. Need to check for (3+1)D and (1+1)D
      for(int idir=0; idir<D; ++idir){
        for(int jdir=0; jdir<D; ++jdir){
          aux_a[idir] += vminia[idir][jdir]*sigsigK[jdir];
          aux_b[idir] += vminib[idir][jdir]*sigsigK[jdir];
        }
      }
      for(int idir=0; idir<D; ++idir){
        Kokkos::atomic_add( &device_hydro_vector(particle_a, ccake::hydro_info::gradshear, idir),
                             vel_dot_grad*(sigsqrb*vb[idir]+sigsqra*va[idir]) );
        Kokkos::atomic_add( &device_hydro_vector(particle_a, ccake::hydro_info::divshear, idir),
                            sigsqrb*aux_b[idir]+sigsqra*aux_a[idir] );
      }
    }
    //Computes gradP weights - Seems relevant only if pressure is NaN
    //double alpha_q    = 1.0;
    //double v_signal_q = sqrt(1.0/3.0);
    //double innerp = EoMPtr<D>::dot( rel_sep, qmom_difference, time_squared );
    //double innerr = EoMPtr<D>::dot( rel_sep, rel_sep, time_squared );
    //innerp = 2.0*alpha_q*v_signal_q
    //      / ( sigma_a/gamma_a + sigma_b/gamma_b )
    //      / sqrt(innerr) * innerp;
//
    //if ( innerp > 0.0 || iparticle == jparticle ) innerp = 0.0;
    //double pressureWeight = entropy_norm_b*sigma_a
    //                        * ( pressure_b / (sigma_b*sigma_b)
    //                          + pressure_a / (sigma_a*sigma_b) - innerp );

  };
  Cabana::neighbor_parallel_for( systemPtr->range_policy, smooth_gradients, systemPtr->neighbour_list, Cabana::FirstNeighborsTag(),
                                   Cabana::TeamOpTag(), "smooth_gradients");
  Kokkos::fence();

  double hc = constants::hbarc_MeVfm;
  auto freeze_out_step = KOKKOS_LAMBDA(const int iparticle){
    int btrack = device_btrack(iparticle);
    double temperature = device_thermo(iparticle, ccake::thermo_info::T)*hc;
    if (( btrack == 1 ) && ( temperature >= 150.0 )) {
      //fo.frz2[pa.ID].t=tin;   // If particle a has only one nearest neighbor (and T>150) set penultimate timestep;
    } ///< TODO: Implement freeze out corrections. Actually I think this is not the best place to implement it
    else if ( (btrack == 0) && (temperature >= 150.0)
            //&& pa.Freeze < 4
    ){
      // otherwise, if a has no nearest neighbors
        // but has T>150 and isn't frozen out,
        // just print this warning message

        //cout << boolalpha << "Missed " << a << " "
         //<< tin << " " << pa.T()*hc << " "
         //<< rdis << " " << systemPtr->do_freeze_out << endl;
    }
  };

  ///TODO: When implemented freeze-out, uncomment this
  //Kokkos::parallel_for("freeze_out_step", systemPtr->n_particles, freeze_out_step);
  //Kokkos::fence();
  std::cout << std::flush;
  return;
}

////////////////////////////////////////////////////////////////////////////////
///\brief Loop over all particles and prepare them for initialization
///
///\details This function loops over all particles and prepares them for
///         initialization. This includes computing gamma factor, defining
///         normalization of norm_spec density and setup of reference density
///         for each particle.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::process_initial_conditions()
{
  formatted_output::report("Processing initial conditions");

  //============================================================================
  // TOTAL NUMBER OF PARTICLES FIXED AFTER THIS POINT
  systemPtr->n_particles = systemPtr->particles.size();

  // fill out initial particle information
  //int TMP_particle_count = 0;
	for (auto & p : systemPtr->particles)
  {

    // set area element for each SPH particle
    double dA = (settingsPtr->stepx)*(settingsPtr->stepy);

    // Set the rest of particle elements using area element
		//p.u(0)          = 0.0;  // flow is set in Particle constructor!!!
		//p.u(1)          = 0.0;  // flow is set in Particle constructor!!!
    p.hydro.gamma     = p.gamcalc();
    p.hydro.gamma_squared        = p.hydro.gamma*p.hydro.gamma;
    p.hydro.gamma_cube        = p.hydro.gamma_squared*p.hydro.gamma;


    // normalize all specific densities to 1
		p.specific.s      = 1.0;
		p.specific.rhoB   = 1.0;
		p.specific.rhoS   = 1.0;
		p.specific.rhoQ   = 1.0;

    // normalization of each density include transverse area element dA
		p.norm_spec.s     = dA;
		p.norm_spec.rhoB  = dA;
		p.norm_spec.rhoS  = dA;
		p.norm_spec.rhoQ  = dA;

		p.hydro.Bulk      = 0.0;

		// make educated initial guess here for this particle's (T, mu_i) coordinates
		// (improve this in the future)
		p.thermo.T        = 800.0/hbarc_MeVfm;	// rootfinder seems to work better going downhill than "uphill"

    // use max T of default EoS method instead
    //p.thermo.T        = (eos.get_default_eos()->tbqs_maxima)[0];

		p.thermo.muB      = 0.0/hbarc_MeVfm;
		p.thermo.muS      = 0.0/hbarc_MeVfm;
		p.thermo.muQ      = 0.0/hbarc_MeVfm;
		p.thermo.eos_name = "default";  // uses whatever the default EoS is

    p.efcheck = systemPtr->efcheck;
		if ( p.input.e > systemPtr->efcheck )	// impose freeze-out check for e, not s
			p.Freeze = 0;
		else
		{
			p.Freeze = 4;
			systemPtr->number_part++;
		}
  }

  formatted_output::detail("particles frozen out: "
                           + to_string(systemPtr->number_part) );
  formatted_output::detail("particles not frozen out: "
                           + to_string(systemPtr->n_particles
                                        - systemPtr->number_part) );


  // make sure freeze-out vectors, etc. are correct size
  freeze_out.resize_vectors( systemPtr->n_particles );


  //============================================================================
  // with particles vector now fully initialized, specify or initialize any
  // remaining quantities which depend on this

  // assign particles IDs
  for ( int i = 0; i < systemPtr->n_particles; i++ )
  {
    systemPtr->particles[i].ID       = i;
    systemPtr->particles[i].hydro.ID = i;
  }

  // set particles to print
  for ( int & p : settingsPtr->particles_to_print )
  {
    if ( p >= systemPtr->n_particles )
    {
      std::cerr << "Particle with ID #" << p << " does not exist, so you cannot"
                   " print it!  Please fix this and re-run." << std::endl;
      abort();
    }
    systemPtr->particles[p].print_this_particle  = true;
    systemPtr->particles[p].hydro.print_particle = true;
  }
}


//==============================================================================
// initialize bulk Pi
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::set_bulk_Pi()
{

  if ( settingsPtr->initializing_with_full_Tmunu ){
    double t0 = settingsPtr->t0;
    CREATE_VIEW( device_, systemPtr->cabana_particles);

    auto set_bulk = KOKKOS_LAMBDA( const int is, const int ia )
    {
      double varsigma = device_hydro_scalar.access(is, ia, ccake::hydro_info::varsigma);
      double p = device_thermo.access(is, ia, ccake::thermo_info::p);
      double u0 = device_hydro_scalar.access(is, ia, ccake::hydro_info::gamma);
      double sigma = device_hydro_scalar.access(is, ia, ccake::hydro_info::sigma);

      device_hydro_scalar.access(is, ia, ccake::hydro_info::Bulk) = (varsigma - p)
                      * u0 * t0 / sigma;
    };
    Cabana::simd_parallel_for( *(systemPtr->simd_policy), set_bulk, "set_bulk_Pi_kernel");
    systemPtr->copy_device_to_host();
    //Kokkos::parallel_for( "set_bulk_Pi", systemPtr->n_particles, set_bulk );
    Kokkos::fence();
  }
  return;
}

/*
//==============================================================================
template<unsigned int D>
void SPHWorkstation<D>::freeze_out_particles()
{
  //---------------------------------------
  // perform freeze out checks
  int n_freezing_out = 0;
  for ( auto & p : systemPtr->particles )
    freeze_out.check_freeze_out_status( p, systemPtr->t, n_freezing_out,
                                systemPtr->n() );

  //---------------------------------------
  // update global quantities accordingly
  systemPtr->number_part += n_freezing_out;
  systemPtr->list.resize(n_freezing_out);

  //---------------------------------------
  // update freeze out status/lists
  int m = 0;
  for ( auto & p : systemPtr->particles )
    if ( p.Freeze == 3 )
    {
      systemPtr->list[m++] = p.ID;
      p.Freeze         = 4;
    }

  //---------------------------------------
  // finalize frozen out particles
  freeze_out.bsqsvfreezeout( n_freezing_out );
}


*/

////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::get_time_derivatives()
{
  Stopwatch sw;
  sw.Start();

  double t = systemPtr->t;
  double t2 = t*t;


  // reset nearest neighbors
  systemPtr->reset_neighbour_list();
  // reset pi tensor to be consistent
  // with all essential symmetries
  reset_pi_tensor(t2);
  // smooth all particle fields - s, rhoB, rhoQ and rhoS and sigma
  smooth_all_particle_fields(t2);
  // update viscosities for all particles
  update_all_particle_viscosities();
  //Computes gradients to obtain dsigma/dt
  smooth_all_particle_gradients(t2);
  //calculate time derivatives needed for equations of motion
  TEOM<D>::evaluate_time_derivatives( systemPtr->cabana_particles, t );

  sw.Stop();
  formatted_output::update("Finished computing time derivatives in "
                            + to_string(sw.printTime()) + " s.");
/*
  // identify and handle particles which have frozen out
  if ( systemPtr->do_freeze_out )
    freeze_out_particles();

  // check/update conserved quantities
  systemPtr->conservation_energy();
  */
  return;
}

template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::update_all_particle_viscosities()
{
  CREATE_VIEW(device_, systemPtr->cabana_particles);

  auto set_transport_coefficients = KOKKOS_CLASS_LAMBDA(const int is, const int ia ){
    double thermo[thermo_info::NUM_THERMO_INFO];
    for (int ithermo=0; ithermo < thermo_info::NUM_THERMO_INFO; ++ithermo)
      thermo[ithermo] = device_thermo.access(is, ia, ithermo);
    device_hydro_scalar.access(is, ia, ccake::hydro_info::setas)
      = tc::eta(thermo, this->transp_coeff_params );
    device_hydro_scalar.access(is, ia, ccake::hydro_info::stauRelax)
      = tc::tau_pi(thermo,this->transp_coeff_params);
    device_hydro_scalar.access(is, ia, ccake::hydro_info::zeta)
      = tc::zeta(thermo, this->transp_coeff_params);
    device_hydro_scalar.access(is, ia, ccake::hydro_info::tauRelax)
      = tc::tau_Pi(thermo, this->transp_coeff_params);
  };
  Cabana::simd_parallel_for(*(systemPtr->simd_policy),set_transport_coefficients, "set_transport_coefficients");
  Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::update_all_particle_thermodynamics(double time_squared)
{
  Stopwatch sw;
  sw.Start();
  double t = systemPtr->t;

  systemPtr->copy_device_to_host();
  ///TODO: This seems trivially parallelizable with openMP. It should be implemented.
  for ( auto & p : systemPtr->particles ){
    double u[D];
    for (int idir = 0; idir < D; ++idir)
      u[idir] = p.hydro.u(idir);
    p.hydro.gamma     = TEOM<D>::gamma_calc(u,time_squared);
    p.hydro.v         = (1.0/p.hydro.gamma)*p.hydro.u;
    double s_LRF      = TEOM<D>::get_LRF(p.smoothed.s, p.hydro.gamma, t);
    double rhoB_LRF   = TEOM<D>::get_LRF(p.smoothed.rhoB, p.hydro.gamma, t);
    double rhoS_LRF   = TEOM<D>::get_LRF(p.smoothed.rhoS, p.hydro.gamma, t);
    double rhoQ_LRF   = TEOM<D>::get_LRF(p.smoothed.rhoQ, p.hydro.gamma, t);
	  locate_phase_diagram_point_sBSQ( p, s_LRF, rhoB_LRF , rhoS_LRF, rhoQ_LRF );
  }
  systemPtr->copy_host_to_device();

  sw.Stop();
  formatted_output::update("got particle thermodynamics in "
                            + to_string(sw.printTime()) + " s.");
}


////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> class TEOM>
double SPHWorkstation<D, TEOM>::locate_phase_diagram_point_eBSQ( Particle<D> & p,
                 double e_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
//cout << "Finding thermodynamics of particle #" << p.ID << endl;

  // default: use particle's current location as initial guess
  // (pass in corresponding EoS as well!)
  eos.tbqs( p.T(), p.muB(), p.muQ(), p.muS(), p.get_current_eos_name() );

  bool solution_found = false;
  double sVal = eos.s_out( e_In, rhoB_In, rhoS_In, rhoQ_In, solution_found,
                           p.print_this_particle );

  if ( solution_found )
  {
    // check that enthalpy derivatives are not crazy
    thermodynamic_info p0 = p.thermo;

    eos.set_thermo( p.thermo );


    ///*
    // do not permit cs2 to be negative if using C library
    // or if using the tanh_conformal EoS
    if (  ( eos.currently_using_static_C_library()
            && p.get_current_eos_name() == "table" )
          || p.get_current_eos_name() == "tanh_conformal" )
    {
      if ( p.thermo.cs2 < 0.0 )
      {
        cout << "WARNING: cs2 went negative in eos_type == "
            << p.get_current_eos_name() << " for particle " << p.ID << "\n"
            << "input thermo: " << e_In << "   "
            << rhoB_In << "   "
            << rhoS_In << "   "
            << rhoQ_In << endl
            << "check thermo: " << systemPtr->t << "   "
            << p.thermo.T << "   "
            << p.thermo.muB << "   "
            << p.thermo.muS << "   "
            << p.thermo.muQ << "   "
            << p.thermo.p << "   "
            << p.thermo.s << "   "
            << p.thermo.rhoB << "   "
            << p.thermo.rhoS << "   "
            << p.thermo.rhoQ << "   "
            << p.thermo.e << "   "
            << p.thermo.cs2 << "   "
            << p.thermo.eos_name << endl;
        p.thermo.cs2 = std::max( p.thermo.cs2, 0.0001 );
      }
    }

    if (p.cs2()<0)
    {
    cout << "input thermo: " << e_In << "   "
          << rhoB_In << "   "
          << rhoS_In << "   "
          << rhoQ_In << endl
          << "check thermo: " << systemPtr->t << "   "
          << p.thermo.T << "   "
          << p.thermo.muB << "   "
          << p.thermo.muS << "   "
          << p.thermo.muQ << "   "
          << p.thermo.p << "   "
          << p.thermo.s << "   "
          << p.thermo.rhoB << "   "
          << p.thermo.rhoS << "   "
          << p.thermo.rhoQ << "   "
          << p.thermo.e << "   "
          << p.thermo.cs2 << "   "
          << p.thermo.eos_name << endl;
      cout << __LINE__ << ": cs2 was negative!" << endl;
      exit(8);
    }

  }

  return sVal;
}

template<unsigned int D, template<unsigned int> class TEOM>
double SPHWorkstation<D, TEOM>::locate_phase_diagram_point_eBSQ(Particle<D> & p, double e_In)
                 { return locate_phase_diagram_point_eBSQ( p, e_In, 0.0, 0.0, 0.0 ); }

////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::locate_phase_diagram_point_sBSQ(Particle<D> & p, double s_In) // previously update_s
               { locate_phase_diagram_point_sBSQ(p, s_In, 0.0, 0.0, 0.0 ); }

/*



///////////////////////////////////////////////////////////////////////////////////
template<unsigned int D>
double SPHWorkstation<D>::gradPressure_weight(const int a, const int b)
{
  auto & pa = systemPtr->particles[a];
  auto & pb = systemPtr->particles[b];

  Vector<double,2> pa_qmom = ( (pa.e()+pa.p())*pa.hydro.gamma/pa.hydro.sigma )*pa.hydro.u;
  Vector<double,2> pb_qmom = ( (pb.e()+pb.p())*pb.hydro.gamma/pb.hydro.sigma )*pb.hydro.u;


  double alpha_q    = 1.0;
  double v_signal_q = sqrt(1.0/3.0);

  double innerp = inner( pa.r - pb.r, pa_qmom - pb_qmom );
  double innerr = inner( pa.r - pb.r, pa.r    - pb.r    );
  innerp = 2.0*alpha_q*v_signal_q
          / ( pa.hydro.sigma/pa.hydro.gamma + pb.hydro.sigma/pb.hydro.gamma )
          / sqrt(innerr) * innerp;

  if ( innerp > 0.0 || a == b ) innerp = 0.0;

  return pb.norm_spec.s * pa.hydro.sigma
        * ( pb.p() / (pb.hydro.sigma*pb.hydro.sigma)
          + pa.p() / (pa.hydro.sigma*pb.hydro.sigma) - innerp );
}


*/
//==============================================================================
// currently add a particle to every grid point which doesn't have one yet
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::add_buffer(double default_e)
{
  if ( settingsPtr->IC_type != "ICCING" )
  {
    std::cerr << "WARNING: ADDING A BUFFER IS CURRENTLY ONLY DEFINED FOR ICCING"
                 " INITIAL CONDITIONS.  NO BUFFER ADDED." << std::endl;
    return;
  }

  constexpr double TINY = 1e-10;

  double xmin  = settingsPtr->xmin;
  double ymin  = settingsPtr->ymin;
  double stepx = settingsPtr->stepx;
  double stepy = settingsPtr->stepy;


  // find maximum distance from origin first
  double max_r = -1.0;
  if ( settingsPtr->circular_buffer )
  {
    for ( Particle<D> & p : systemPtr->particles )
      max_r = std::max( max_r, Norm(p.r) );
    max_r *= 1.0 + settingsPtr->padding_thickness;

    xmin = -stepx * ceil( max_r / stepx );
    ymin = -stepy * ceil( max_r / stepy );

    if ( max_r < 0.0 )
    {
      cout << "max_r = " << max_r << " cannot be less than zero!" << endl;
      abort();
    }
  }


  const int nx = 1 - 2*int(round(xmin/stepx));  // xmin is negative
  const int ny = 1 - 2*int(round(ymin/stepy));  // ymin is negative
  bool particle_exists[nx][ny];
  for ( int ix = 0; ix < nx; ix++ )
  for ( int iy = 0; iy < ny; iy++ )
    particle_exists[ix][iy] = false;


  // specify which particles are already in grid
  int particle_count = 0;
  for ( auto & p : systemPtr->particles )
  {
    int ix = int(round((p.r(0)-xmin)/stepx));
    int iy = int(round((p.r(1)-ymin)/stepy));
    particle_exists[ix][iy] = true;
    particle_count++;
  }

  cout << "nx = " << nx << endl;
  cout << "ny = " << ny << endl;
  cout << "particle_count = " << particle_count << endl;

  int checked_particles = 0;
  int added_particles = 0;

  // loop over all points to consider,
  // initialize those not yet in grid
  // (if in padded circle, for circular_buffer == true)
  for ( int ix0 = 0; ix0 < nx; ix0++ )
  for ( int iy0 = 0; iy0 < ny; iy0++ )
  {
    checked_particles++;

    // don't initialize particles that already exist!
    if ( particle_exists[ix0][iy0] ) continue;

    double x0 = xmin + ix0*stepx;
    double y0 = ymin + iy0*stepy;

    // if circular_buffer == true,
    // don't initialize particles outside padded circle!
    if ( settingsPtr->circular_buffer && sqrt(x0*x0+y0*y0) > max_r ) continue;

    Particle<D> p;

    p.r(0)       = x0;
    p.r(1)       = y0;
    p.input.e    = default_e;
    p.input.rhoB = 0.0;
    p.input.rhoS = 0.0;
    p.input.rhoQ = 0.0;
    p.hydro.u(0) = 0.0;
    p.hydro.u(1) = 0.0;

    systemPtr->particles.push_back( p );
    added_particles++;

  }

  cout << "checked_particles = " << checked_particles << endl;
  cout << "added_particles = " << added_particles << endl;

  return;
}

//============================================================================
// decide whether to continue evolving
template<unsigned int D, template<unsigned int> class TEOM>
bool SPHWorkstation<D, TEOM>::continue_evolution()
{
  std::cout << "t = " << systemPtr->t << std::endl;
  return ( systemPtr->t < settingsPtr->tend )
          && ( systemPtr->number_part < systemPtr->n() );
}

//============================================================================
// Advance one timestep
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::advance_timestep( double dt, int rk_order )
{
  Stopwatch sw;
  sw.Start();
  // turn on freeze-out flag initially
  //systemPtr->do_freeze_out = true;
  ///TODO: do_freeze_out should be a flag in the settings file
  // use evolver to actually do RK evolution
  // (pass workstation's own time derivatives function as lambda)

  //Bulk of code evaluation is done below
  evolver.execute_timestep( dt, rk_order,
                            [this]{ this->get_time_derivatives(); } );

  // set number of particles which have frozen out
  //systemPtr->number_part = systemPtr->get_frozen_out_count();
  // keep track of how many timesteps have elapsed
  systemPtr->number_of_elapsed_timesteps++;

  sw.Stop();
  formatted_output::report("finished timestep in "
                            + to_string(sw.printTime()) + " s");
  if ( settingsPtr->max_number_of_timesteps >= 0
        && systemPtr->number_of_elapsed_timesteps
            > settingsPtr->max_number_of_timesteps )
    exit(-1);
  return;
}
};