#include "sph_workstation.h"
//The cpp below needs to be included for the Kokkos lambda to be
//properly compiled in the GPU
#include "kernel.cpp"
#include "eom_default.cpp"
#include "eom_cartesian.cpp"
#include "transport_coefficients.cpp"

#define ONLINE_INVERTER ///todo: This should be in the config file. For now, comment, uncomment it to choose eos method to be used
using namespace constants;
namespace tc = ccake::transport_coefficients;

namespace ccake{
//Template instantiations
template class SPHWorkstation<1,EoM_default>;
template class SPHWorkstation<2,EoM_default>;
template class SPHWorkstation<3,EoM_default>;
template class SPHWorkstation<1,EoM_cartesian>;
template class SPHWorkstation<2,EoM_cartesian>;
template class SPHWorkstation<3,EoM_cartesian>;

/// @brief Initialize the several components of the SPHWorkstation.
/// @details SPHWorstation needs to load and initiualize a few objects before 
/// @details SPHWorstation needs to load and initiualize a few objects before
/// it can work. Those are
/// - The equation of state
/// - The system state (containing the SPH particles)
/// - The transport coefficients
/// - The freeze out energy density
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
  if(!settingsPtr->online_inverter_enabled)
    eos_interpolatorPtr = std::make_shared<EoS_Interpolator>(
          settingsPtr->preinverted_eos_path);
  //----------------------------------------
  // set up transport coefficients
  transp_coeff_params = tc::setup_parameters(settingsPtr);

  //----------------------------------------
  // set up freeze out  at constant energy density
  formatted_output::report("Setting up freeze out");
  formatted_output::detail("Freeze out temperature = "
                             + to_string(settingsPtr->Freeze_Out_Temperature) + " MeV");
  systemPtr->efcheck = eos.efreeze(settingsPtr->Freeze_Out_Temperature/hbarc_MeVfm); //Factor 1000 to convert to MeV from GeV
  formatted_output::detail("freeze out energy density = "
                           + to_string(systemPtr->efcheck*hbarc_GeVfm) + " GeV/fm^3");

}


/// @brief This function will ensure that the shear tensor has the correct
/// properties.
/// @details This function will ensure that the shear tensor has the correct
/// properties. The desired properties are that it is traceless, symmetric and 
/// orthogonal to the fluid velocity. That is, we want to ensure that 
/// \f$\pi^{\mu \nu} = \pi^{\nu \mu}\f$, \f$u_\mu \pi^{\mu \nu} = 0\f$ and
/// \f$\pi^\mu_\mu = 0\f$. Because this depends on the metric tensor being used,
/// a simple call to the reset_pi_tensor function is not enough. This function
/// should be implemented in the templates TEOM equation of motion class.
/// @todo We do not need the time_squared parameter here. We should remove it.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation. Must
/// have a reset_pi_tensor function.
/// @param time_squared The square of the time step where the shear tensor is
/// being computed.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::reset_pi_tensor(double time_squared)
{
  TEOM<D>::reset_pi_tensor(systemPtr);
};

// @brief calculate gamma and velocties
/// @details This function calculates the Lorentz contraction factor \f$ \gamma \f$ 
/// and the fluid velocity \f$ v \f$ from the fluid four velocity \f$ u \f$.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::calculate_gamma_and_velocities()
{
 TEOM<D>::update_velocity(systemPtr);
};


/// @brief Set up entropy density and other thermodynamic variables
/// @details This function sets up the entropy density using the initial energy
/// density and the initial charge densities via the equation of state. It also
/// sets up all other thermodynamic variables. The sph_mass densities are also
/// updated here.
///
/// sph_mass densities are are defined as
/// \f[
/// \rho_{\text{sph_mass}} = 1
/// \f]
/// @todo At the moment, the equation of state functions do not support running
/// in device. Because of this, a copy of the particle data is made to the host
/// and the equation of state is run on the host. Then the data is copied back
/// to the device. This is can be a performance bottleneck, specially for large
/// simulations. The equation of state functions should be modified to run on
/// the device.
/// @todo The formula for the sph_mass densities is not general. An
/// implementation inside the TEOM class should be used instead.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::initialize_entropy_and_charge_densities()
{
	Stopwatch sw, swTotal;
	swTotal.Start();

  //Compute thermal properties
  systemPtr->copy_device_to_host();
  if(settingsPtr->input_as_entropy==true){
    for (auto & p : systemPtr->particles){
      locate_phase_diagram_point_sBSQ( p,
                    p.input.s, p.input.rhoB, p.input.rhoS, p.input.rhoQ );
      p.input.e = p.thermo.e;
      p.efcheck = systemPtr->efcheck;
      if ( p.input.e > systemPtr->efcheck )	p.Freeze = 0;
		  else
		    {
		    	p.Freeze = 4;
		    	systemPtr->number_part_fo++;
		    }
    }
  }
  else{
    for (auto & p : systemPtr->particles){
      p.input.s = locate_phase_diagram_point_eBSQ( p,
                    p.input.e, p.input.rhoB, p.input.rhoS, p.input.rhoQ );
      p.efcheck = systemPtr->efcheck;
      if ( p.input.e > systemPtr->efcheck )	p.Freeze = 0;
		  else
		    {
		    	p.Freeze = 4;
		    	systemPtr->number_part_fo++;
		    }
    }
  }
  
  formatted_output::detail("particles frozen out: "
                           + to_string(systemPtr->number_part_fo) );
  formatted_output::detail("particles not frozen out: "
                           + to_string(systemPtr->n_particles
                                        - systemPtr->number_part_fo) );



  systemPtr->copy_host_to_device();

  double t0 = settingsPtr->t0;
  // loop over all particles
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  //calculate initial volumes
  calculate_intial_sigma(t0);

  auto init_particles = KOKKOS_LAMBDA(const int iparticle)
  {
    double s_input = device_input(iparticle, densities_info::s);
    double rhoB_input = device_input(iparticle, densities_info::rhoB);
    double rhoS_input = device_input(iparticle, densities_info::rhoS);
    double rhoQ_input = device_input(iparticle, densities_info::rhoQ);
    double gamma = device_hydro_scalar(iparticle, hydro_info::gamma);
    double s_sph_mass = device_sph_mass(iparticle, densities_info::s);
    double rhoB_sph_mass = device_sph_mass(iparticle, densities_info::rhoB);
    double rhoS_sph_mass = device_sph_mass(iparticle, densities_info::rhoS);
    double rhoQ_sph_mass = device_sph_mass(iparticle, densities_info::rhoQ);
    double sigma = device_hydro_scalar(iparticle, hydro_info::sigma);


    device_sph_mass(iparticle, ccake::densities_info::s )    = 1;    // constant after this
    device_sph_mass(iparticle, ccake::densities_info::rhoB ) = 1; // constant after this
    device_sph_mass(iparticle, ccake::densities_info::rhoS ) = 1; // constant after this
    device_sph_mass(iparticle, ccake::densities_info::rhoQ ) = 1; // constant after this
    device_extensive(iparticle, ccake::densities_info::s )    = s_input/sigma;
    device_extensive(iparticle, ccake::densities_info::rhoB ) = rhoB_input/sigma;
    device_extensive(iparticle, ccake::densities_info::rhoS ) = rhoS_input/sigma;
    device_extensive(iparticle, ccake::densities_info::rhoQ ) = rhoQ_input/sigma;
    if (s_input < 0.0)
      device_freeze(iparticle) = 4;

  };
  Kokkos::parallel_for("init_particles", systemPtr->n_particles, init_particles);
  Kokkos::fence();
  //Allocate cache for evolver
  evolver.allocate_cache();



	swTotal.Stop();
  formatted_output::update("finished initializing particle densities in "
                              + to_string(swTotal.printTime()) + " s");

}

/// @brief Smooth all SPH fields
/// @details This function updates the densities \f$s\f$, \f$\rho_B\f$, 
/// \f$rho_Q\f$, \f$\rho_S\f$ and \f$\sigma_lab\f$ Iauxiliary density) by performing
/// the smoothing procedure. At its end, it
/// takes the opportunity to update the remaining thermodynamic quantities
/// (energy density, pressure, chemical potentials, temperature and speed of
/// sound squared). The routine that triggers the update of thermo quantities
/// will also update gamma and the non-relativistic velocity. At the end, it 
/// also computes the freeze out status of the particles. In this function, we
/// also take the opportunity to enforce the constraints for the shear viscous
/// tensor.
/// @see reset_pi_tensor
/// @see smooth_all_particle_fields
/// @see update_all_particle_thermodynamics
/// @see FreezeOut::check_freeze_out_status
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D,  template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::initial_smoothing()
{
  Stopwatch sw;
  sw.Start();
  double t_squared = pow(settingsPtr->t0,2);
  // calculate gamma and velocities
  calculate_gamma_and_velocities();
  /*#ifdef DEBUG
  std::cout << "initial_smoothing() --> reset_p_tensor()" << std::endl;
  for (auto & p : systemPtr->particles){
    if(p.r(0) > -0.05 && p.r(0) < 0.05){
      std::cout << "eta=  " << p.r(0) << "   ";
      std::cout << "thermo.e=   " << p.thermo.e*hbarc_GeVfm << "   input.e=   " << p.input.e*hbarc_GeVfm << \
              "   thermo.s=   " << p.thermo.s << "   input.s=   " << p.input.s << "   " << std::endl;
    }
  }
  #endif*/
  // smooth fields over particles
  smooth_all_particle_fields(t_squared);
  //calculate extensive/extensive shear tensor
  calculate_extensive_shv();
  /* #ifdef DEBUG
  std::cout << "initial_smoothing() --> smooth_all_particle_fields()" << std::endl;
  for (auto & p : systemPtr->particles){
    if(p.r(0) > -0.05 && p.r(0) < 0.05){
      std::cout << "eta=  " << p.r(0) << "   ";
      std::cout << "thermo.e=   " << p.thermo.e*hbarc_GeVfm << "   input.e=   " << p.input.e*hbarc_GeVfm << \
              "   thermo.s=   " << p.thermo.s << "   input.s=   " << p.input.s << "   " << std::endl;
    }
  }
  #endif */
  //Computes not independent components of the shear tensor
  reset_pi_tensor(t_squared);
  // Update particle thermodynamic properties
  update_all_particle_thermodynamics();
  //std::cout << "check initial conservation" << std::endl;
  systemPtr->conservation_energy(true,settingsPtr->t0);
  //std::cout << "check initial conservation done" << std::endl;

  /* #ifdef DEBUG
  std::cout << "initial_smoothing() --> update_all_particle_thermodynamics()" << std::endl;
  for (auto & p : systemPtr->particles){
    if(p.r(0) > -0.05 && p.r(0) < 0.05){
      std::cout << "eta=  " << p.r(0) << "   ";
      std::cout << "thermo.e=   " << p.thermo.e*hbarc_GeVfm << "   input.e=   " << p.input.e*hbarc_GeVfm << \
              "   thermo.s=   " << p.thermo.s << "   input.s=   " << p.input.s << "   " << std::endl;
    }
  }
  exit(1);
  #endif */

  //Performs the initial freeze-out
  int count1=0;
  freezePtr->check_freeze_out_status(count1);
  sw.Stop();
  formatted_output::update("Finished initial smoothing "
                            + to_string(sw.printTime()) + " s.");

}

///@brief Smooth all SPH fields
///@details This function updates the densities \f$s\f$, \f$\rho_B\f$, 
/// \f$\rho_Q\f$, \f$\rho_S\f$ and \f$\sigma_lab\f$ (auxiliary density) by 
/// performing the smoothing procedure. At its end, it takes the opportunity to
/// update the remaining thermodynamic quantities (energy density, pressure, 
/// chemical potentials, temperature and speed of  sound squared). The routine 
/// that triggers the update of thermo quantities will also update gamma and 
///@details This function updates the densities \f$s\f$, \f$\rho_B\f$,
/// \f$\rho_Q\f$, \f$\rho_S\f$ and \f$\sigma\f$ (auxiliary density) by
/// performing the smoothing procedure. At its end, it takes the opportunity to
/// update the remaining thermodynamic quantities (energy density, pressure,
/// chemical potentials, temperature and speed of  sound squared). The routine
/// that triggers the update of thermo quantities will also update gamma and
/// the non-relativistic velocity
///
/// The computation of the smoothed fields is done using the following formula:
/// \f[
/// \rho_{\text{smoothed},\,i} = \sum_{j = \text{nearest neighbors}}
/// \rho_{\text{sph_mass},\,i} \rho_{\text{spec},\,j} W_{ij}
/// \f]
/// where \f$\rho_{\text{sph_mass},\,i}\f$ is the normalized extensive density
/// of particle \f$i\f$, \f$\rho_{\text{spec},\,j}\f$ is the extensive density
/// of particle \f$j\f$ and \f$W_{ij}\f$ is the kernel function evaluated at
/// the distance between particles \f$i\f$ and \f$j\f$. The list of nearest
/// neighbors over which the sum is performed includes the particle \f$i\f$
/// itself.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
/// @param time_squared The square of the current time step.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::smooth_all_particle_fields(double time_squared)
{
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
  double hT = settingsPtr->hT;
  double t = systemPtr->t;

  //creating outfile to for kernal function output
  /*ofstream outfile;
  std::fname = "kernel_" + std::to_string(hT) + ".dat";
  outfile.open(fname);*/

  //Reset smoothed fields. Initializes using the contribution of the particle to the density evaluated
  //on top of itself because the loop over the neighbour particles does not include the particle itself.
  double kern0 = SPHkernel<D>::kernel(0,hT); //The value of the Kernel evaluated on top of the particle
  auto reset_fields = KOKKOS_LAMBDA(const int is, const int ia) //First index for loop over struct, second for loop over array
  {
    device_smoothed.access(is, ia, densities_info::s)     = device_sph_mass.access(is, ia, densities_info::s)*   device_extensive.access(is, ia, densities_info::s)*kern0;
    device_smoothed.access(is, ia, densities_info::rhoB)  = device_sph_mass.access(is, ia, densities_info::rhoB)*device_extensive.access(is, ia, densities_info::rhoB)*kern0;
    device_smoothed.access(is, ia, densities_info::rhoQ)  = device_sph_mass.access(is, ia, densities_info::rhoQ)*device_extensive.access(is, ia, densities_info::rhoQ)*kern0;
    device_smoothed.access(is, ia, densities_info::rhoS)  = device_sph_mass.access(is, ia, densities_info::rhoS)*device_extensive.access(is, ia, densities_info::rhoS)*kern0;
    device_hydro_scalar.access(is, ia, hydro_info::sigma_lab) = device_sph_mass.access(is, ia, densities_info::s)*kern0;
  };
  Cabana::simd_parallel_for( simd_policy, reset_fields, "reset_fields" );
  Kokkos::fence();

  auto smooth_fields = KOKKOS_LAMBDA(const int iparticle, const int jparticle ){

    double r1[D] ,r2[D]; // cache for positions of particles 1 and 2
    for (int idir = 0; idir < D; ++idir){
      r1[idir] = device_position(iparticle,idir);
      r2[idir] = device_position(jparticle,idir);
    }
    double distance = SPHkernel<D>::distance(r1,r2);
    double kern = SPHkernel<D>::kernel(distance,hT);

    //outputting kernel function to the outfile "kernel_{hT}.dat"
    /*outfile << kern << distance << endl;
    outfile.close();*/

    //Update sigma_lab (reference density)
    Kokkos::atomic_add( &device_hydro_scalar(iparticle, ccake::hydro_info::sigma_lab), device_sph_mass(jparticle, ccake::densities_info::s)*kern);
    ////Update entropy density
    Kokkos::atomic_add( &device_smoothed(iparticle, ccake::densities_info::s),
                         device_sph_mass(jparticle, ccake::densities_info::s)*device_extensive(jparticle, ccake::densities_info::s)*kern);
    Kokkos::atomic_add( &device_smoothed(iparticle, ccake::densities_info::rhoB),
                          device_sph_mass(jparticle, ccake::densities_info::rhoB)*device_extensive(jparticle, ccake::densities_info::rhoB)*kern);
    Kokkos::atomic_add( &device_smoothed(iparticle, ccake::densities_info::rhoS),
                          device_sph_mass(jparticle, ccake::densities_info::rhoS)*device_extensive(jparticle, ccake::densities_info::rhoS)*kern);
    Kokkos::atomic_add( &device_smoothed(iparticle, ccake::densities_info::rhoQ),
                          device_sph_mass(jparticle, ccake::densities_info::rhoQ)*device_extensive(jparticle, ccake::densities_info::rhoQ)*kern);
  };

  auto range_policy = Kokkos::RangePolicy<ExecutionSpace>(0, systemPtr->cabana_particles.size());
  Cabana::neighbor_parallel_for( range_policy, smooth_fields,
                                 systemPtr->neighbour_list, Cabana::FirstNeighborsTag(),
                                 Cabana::TeamOpTag(), "smooth_fields_kernel");

  Kokkos::fence();
  auto update_input_thermodynamics = KOKKOS_LAMBDA(const int is, const int ia){
    double s = device_smoothed.access(is, ia, densities_info::s);
    double rhoB = device_smoothed.access(is, ia, densities_info::rhoB);
    double rhoS = device_smoothed.access(is, ia, densities_info::rhoS);
    double rhoQ = device_smoothed.access(is, ia, densities_info::rhoQ);
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double sigma_lab = device_hydro_scalar.access(is, ia, hydro_info::sigma_lab);

    double sigma  = TEOM<D>::get_LRF(sigma_lab, gamma, t);
    s = TEOM<D>::get_LRF(s, gamma, t);
    rhoB = TEOM<D>::get_LRF(rhoB, gamma, t);
    rhoS = TEOM<D>::get_LRF(rhoS, gamma, t);
    rhoQ = TEOM<D>::get_LRF(rhoQ, gamma, t);
    device_hydro_scalar.access(is, ia, hydro_info::sigma) = sigma;
    device_thermo.access(is, ia, thermo_info::s) = s;
    device_thermo.access(is, ia, thermo_info::rhoB) = rhoB;
    device_thermo.access(is, ia, thermo_info::rhoS) = rhoS;
    device_thermo.access(is, ia, thermo_info::rhoQ) = rhoQ;
  };
  Cabana::simd_parallel_for( simd_policy, update_input_thermodynamics, "update_input_thermodynamics" );
  Kokkos::fence();

}



template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::calculate_intial_sigma(double t)
{
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
  double hT = settingsPtr->hT;
  //creating outfile to for kernal function output
  /*ofstream outfile;
  std::fname = "kernel_" + std::to_string(hT) + ".dat";
  outfile.open(fname);*/

  //Reset smoothed fields. Initializes using the contribution of the particle to the density evaluated
  //on top of itself because the loop over the neighbour particles does not include the particle itself.
  double kern0 = SPHkernel<D>::kernel(0,hT); //The value of the Kernel evaluated on top of the particle
  auto reset_fields = KOKKOS_LAMBDA(const int is, const int ia) //First index for loop over struct, second for loop over array
  {
    device_hydro_scalar.access(is, ia, hydro_info::sigma_lab) = device_sph_mass.access(is, ia, densities_info::s)*kern0;
  };
  Cabana::simd_parallel_for( simd_policy, reset_fields, "reset_fields" );
  Kokkos::fence();

  auto smooth_fields = KOKKOS_LAMBDA(const int iparticle, const int jparticle ){

    double r1[D] ,r2[D]; // cache for positions of particles 1 and 2
    for (int idir = 0; idir < D; ++idir){
      r1[idir] = device_position(iparticle,idir);
      r2[idir] = device_position(jparticle,idir);
    }
    double distance = SPHkernel<D>::distance(r1,r2);
    double kern = SPHkernel<D>::kernel(distance,hT);

    //outputting kernel function to the outfile "kernel_{hT}.dat"
    /*outfile << kern << distance << endl;
    outfile.close();*/

    //Update sigma_lab (reference density)
    Kokkos::atomic_add( &device_hydro_scalar(iparticle, ccake::hydro_info::sigma_lab), device_sph_mass(jparticle, ccake::densities_info::s)*kern);
  };

  auto range_policy = Kokkos::RangePolicy<ExecutionSpace>(0, systemPtr->cabana_particles.size());
  Cabana::neighbor_parallel_for( range_policy, smooth_fields,
                                 systemPtr->neighbour_list, Cabana::FirstNeighborsTag(),
                                 Cabana::TeamOpTag(), "smooth_fields_kernel");

  Kokkos::fence();
  auto update_sigma = KOKKOS_LAMBDA(const int is, const int ia){
    double gamma = device_hydro_scalar.access(is, ia, hydro_info::gamma);
    double sigma_lab = device_hydro_scalar.access(is, ia, hydro_info::sigma_lab);

    double sigma  = TEOM<D>::get_LRF(sigma_lab, gamma, t);
    device_hydro_scalar.access(is, ia, hydro_info::sigma) = sigma;

  };
  Cabana::simd_parallel_for( simd_policy, update_sigma, "update_input_thermodynamics" );
  Kokkos::fence();

}


/// @brief Shell function to update thermodynamic particle properties.
/// @details This will updated Particle<D>::thermo for the input partiucle,
/// assuming the smoothed fields \f$s\f$, \f$\rho_B\f$, \f$\rho_S\f$, 
/// \f$\rho_Q\f$ are known.
///
/// @todo I am not an expert in the EoS part of the code. Someone, please
/// comnplement this documentation.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
/// @param p The particle for which the thermodynamic properties are to be
/// updated.
/// @param s_In The entropy density of the particle.
/// @param rhoB_In The baryon density of the particle.
/// @param rhoS_In The strangeness density of the particle.
/// @param rhoQ_In The electric charge density of the particle.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::locate_phase_diagram_point_sBSQ( Particle<D> & p,
                 double s_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{

  // default: use particle's current location as initial guess
  #ifdef DEBUG_SLOW
  std::cout << "In locate_phase_diagram_point_sBSQ" << std::endl;
  for (auto & p : systemPtr->particles){
    if(p.r(0) > -0.01 && p.r(0) < 0.01){
      std::cout << "eta=  " << p.r(0) << "   ";
      std::cout << "  temp=  " << p.thermo.T << "  thermo.e=   " << p.thermo.e*hbarc_GeVfm << "   input.e=   " << p.input.e*hbarc_GeVfm << \
              "   thermo.s=   " << p.thermo.s << "   input.s=   " << p.input.s << "   " << std::endl;
    }
  }
  #endif

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
  #ifdef DEBUG_SLOW
  std::cout << "End of locate_phase_diagram_point_sBSQ" << std::endl;
  for (auto & p : systemPtr->particles){
    if(p.r(0) > -0.02 && p.r(0) < 0.02){
      std::cout << "eta=  " << p.r(0) << "   ";
      std::cout << "  temp=  " << p.thermo.T << "  thermo.e=   " << p.thermo.e*hbarc_GeVfm << "   input.e=   " << p.input.e*hbarc_GeVfm << \
              "   thermo.s=   " << p.thermo.s << "   input.s=   " << p.input.s << "   " << \
              "   eos:   " << p.thermo.eos_name << std::endl;
    }
  }
  #endif

  return;
}

/// @brief Computes all the gradients necessary for the hydrodynamic evolution.
/// @details This function computes the gradients which are necessary for
/// evaluating the time derivatives of the hydrodynamic fields. The gradients
/// computed are `gradP`, `gradE`, `gradBulk`, `gradV`, `gradshear` and
/// `divshear`. The expression for these are, respectivelly:
/// \f[\begin{align*}
/// \partial_i P_a & = \sum_{
///                    \substack{ b=\text{Nearest} \\ \text{neighbors}}
///                   }
/// \sigma_lab_a^* s_{\text{sph_mass},\,b} \left( \frac{P_b}{{\sigma_lab_b^*}^2} +
/// \frac{P_a}{{\sigma_lab_a^*}^2} \right) \partial_{i,\,a} W(\vec r_a - \vec r_b,\,h) \\
/// \partial_i E_a & = \sum_{
///                    \substack{ b=\text{Nearest} \\ \text{neighbors}}
///                   }
/// \sigma_lab_a^* s_{\text{sph_mass},\,b} \left( \frac{E_b}{{\sigma_lab_b^*}^2} +
/// \frac{E_a}{{\sigma_lab_a^*}^2} \right)
/// \partial_{i,\,a} W(\vec r_a - \vec r_b,\,h) \\
/// \partial_i \Pi_a & = \sum_{
///                    \substack{ b=\text{Nearest} \\ \text{neighbors}}
///                   }
/// \frac{\sigma_lab_a^* s_{\text{sph_mass},\,b}}{t}
/// \left( \frac{\tilde\Pi_a}{\gamma_a {\sigma_lab_a^*}} +
///        \frac{\tilde\Pi_b}{\gamma {\sigma_lab_b^*}} \right)
/// \partial_{i,\,a} W(\vec r_a - \vec r_b,\,h) \\
/// \partial_j v^i_a & = \sum_{
///                    \substack{ b=\text{Nearest} \\ \text{neighbors}}
///                   }
/// \frac{s_{\text{sph_mass},\,b}}{\sigma_lab_a^*}(v_b^i-v_a^i)
/// \partial_{j,\,a} W(\vec r_a - \vec r_b,\,h) \\
/// v^j\partial_j\pi^{0i}_a & = \sum_{
///                    \substack{ b=\text{Nearest} \\ \text{neighbors}}
///                   }
/// s_{\text{sph_mass},\,b} \sigma_lab_a^*
/// \left[\frac{\pi^{0i}_b}{{\sigma_lab^*_b}^2}
///       +\frac{\pi^{0i}_a}{{\sigma_lab^*_a}^2} \right]
/// v^j_a \partial_{j,\,a} W(\vec{r}_a-\vec{r}_b,h) \\
/// \partial_j\pi^{ij}_a & = \sum_{
///                    \substack{ b=\text{Nearest} \\ \text{neighbors}}
///                   }
/// s_{\text{sph_mass},\,b} \sigma_lab_a^* \left[
/// \frac{\pi^{ij}_a}{{\sigma_lab^*_a}^2}+\frac{\pi^{ij}_b}{{\sigma_lab^*_b}^2}\right]
/// \partial_{j,\,a} W(\vec{r}_a-\vec{r}_b,h)
/// \end{align*}\f]
/// @todo Notice the denominator of `gradBulk` could use the replacement
/// \f$\sigma = \sigma^* \gamma t\f$ and
/// we would have a more general expression. Time should not be a parameter of
/// the function.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
/// @param time_squared The square of the current time step.
/// Cartesian specializations at end of file 
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::smooth_all_particle_gradients(double time_squared)
{
  double hT = settingsPtr->hT;
  double t = systemPtr->t;
  bool using_shear = settingsPtr->using_shear;
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  CREATE_VIEW(device_, systemPtr->cabana_particles);
  //Reset gradients - No need to take into account self interaction, since gradient of Kernel for r = 0 is zero
  auto reset_gradients = KOKKOS_LAMBDA(const int is, const int ia){
    for (int idir=0; idir < D; ++idir)
    {
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradP, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradBulk, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradshear, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::divshear, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradE,idir)  = 0.0;
      for (int jdir=0; jdir < D; ++jdir)
        device_hydro_space_matrix.access(is, ia,ccake::hydro_info::gradV, idir, jdir) = 0.0;
    }
  };
  Cabana::simd_parallel_for(simd_policy,reset_gradients, "reset_gradients");
  Kokkos::fence();

  ///particle_a is the one that will be updated. Particle b is its neighbors.
  auto smooth_gradients = KOKKOS_LAMBDA(const int particle_a, const int particle_b )
  {
    //Compute kernel gradient
    double rel_sep[D], pos_a[D], pos_b[D]; // cache for relative separation
    for (int idir = 0; idir < D; ++idir){
      pos_a[idir] = device_position(particle_a,idir);
      pos_b[idir] = device_position(particle_b,idir);
      rel_sep[idir] = pos_a[idir] - pos_b[idir];
    }
    double rel_sep_norm = SPHkernel<D>::distance(pos_a,pos_b);
    double gradK_aux[D];
    ccake::SPHkernel<D>::gradKernel( rel_sep, rel_sep_norm, hT, gradK_aux );

    double pressure_a     = device_thermo(particle_a, thermo_info::p);
    double pressure_b     = device_thermo(particle_b, thermo_info::p);
    double energy_a       = device_thermo(particle_a, thermo_info::e);
    double energy_b       = device_thermo(particle_b, thermo_info::e);
    double sigma_lab_a        = device_hydro_scalar(particle_a, hydro_info::sigma_lab);
    double sigma_lab_b        = device_hydro_scalar(particle_b, hydro_info::sigma_lab);
    double sigsqra        = 1.0/(sigma_lab_a*sigma_lab_a);
    double sigsqrb        = 1.0/(sigma_lab_b*sigma_lab_b);
    double sph_mass_s_b  = device_sph_mass(particle_b, densities_info::s);
    double bulk_a         = device_hydro_scalar(particle_a, hydro_info::bulk);
    double bulk_b         = device_hydro_scalar(particle_b, hydro_info::bulk);
    double gamma_a        = device_hydro_scalar(particle_a, hydro_info::gamma);
    double gamma_b        = device_hydro_scalar(particle_b, hydro_info::gamma);

    milne::Vector<double,D> gradK, sigsigK, gradP, gradE, v_a, v_b, gradBulk, gradshear, divshear;
    milne::Matrix<double, D, D> gradV;
    for(int idir=0; idir<D; ++idir){
      gradK(idir) = gradK_aux[idir];
      v_a(idir) = device_hydro_vector(particle_a, hydro_info::v, idir);
      v_b(idir) = device_hydro_vector(particle_b, hydro_info::v, idir);
    };

    sigsigK = sph_mass_s_b * sigma_lab_a * gradK;
    gradP = (sigsqrb*pressure_b + sigsqra*pressure_a ) * sigsigK;
    gradE = (sph_mass_s_b/sigma_lab_a)*( energy_b-energy_a )*gradK;
    //gradE = (sigsqrb*energy_b + sigsqra*energy_a ) * sigsigK;
    gradBulk = (bulk_b/sigma_lab_b/sigma_lab_b + bulk_a/sigma_lab_a/sigma_lab_a)*sigsigK;
    //gradV    = (sph_mass_s_b/sigma_lab_a)*( v_b -  v_a )*gradK;
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir)
      gradV(idir,jdir) = (sph_mass_s_b/sigma_lab_a)*( v_b(jdir) -  v_a(jdir) )*gradK(idir);


    if (using_shear){
      if (D > 1) {
        milne::Matrix<double,D+1,D+1> shv_a, shv_b;
        for(int idir=0; idir<D+1;++idir)
        for(int jdir=0; jdir<D+1;++jdir){
          shv_a(idir,jdir) = device_hydro_spacetime_matrix(particle_a, hydro_info::shv, idir, jdir);
          shv_b(idir,jdir) = device_hydro_spacetime_matrix(particle_b, hydro_info::shv, idir, jdir);
        }
        milne::Vector<double,D> shv_0i_a  = milne::rowp1(0, shv_a);
        milne::Vector<double,D> shv_0i_b  = milne::rowp1(0, shv_b);

        milne::Matrix<double,D,D> vminia, vminib;
        milne::mini(vminia, shv_a);
        milne::mini(vminib, shv_b);
        gradshear = milne::contract(sigsigK, v_a)*( sigsqrb*shv_0i_b + sigsqra*shv_0i_a );
        divshear  = sigsqrb*sigsigK*milne::transpose(vminib)
                  + sigsqra*sigsigK*milne::transpose(vminia);
      }
      else{
        double shv_0i_a  = device_hydro_spacetime_matrix(particle_a, hydro_info::shv, 0, 3);
        double shv_0i_b  = device_hydro_spacetime_matrix(particle_b, hydro_info::shv, 0, 3);
        milne::Matrix<double,D,D> vminia, vminib;
        vminia(0,0) = device_hydro_spacetime_matrix(particle_a, hydro_info::shv, 3, 3);
        vminib(0,0) = device_hydro_spacetime_matrix(particle_b, hydro_info::shv, 3, 3);
        gradshear = milne::contract(sigsigK, v_a)*( sigsqrb*shv_0i_b + sigsqra*shv_0i_a );
        ///\todo do we need this transpose here?
        divshear  = sigsqrb*sigsigK*vminib
                  + sigsqra*sigsigK*vminia;        
      }

    }

    //Accumulate
    for(int idir=0; idir<D; ++idir) {
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradP, idir), gradP(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradE, idir), gradE(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradBulk, idir), gradBulk(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradshear, idir), gradshear(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::divshear, idir), divshear(idir));
      for(int jdir=0; jdir<D; ++jdir)
        Kokkos::atomic_add(&device_hydro_space_matrix(particle_a, hydro_info::gradV, idir, jdir), gradV(idir, jdir));
    }
  };
  auto range_policy = Kokkos::RangePolicy<ExecutionSpace>(0, systemPtr->cabana_particles.size());
  Cabana::neighbor_parallel_for(range_policy, smooth_gradients,
                                 systemPtr->neighbour_list, Cabana::FirstNeighborsTag(),
                                 Cabana::TeamOpTag(), "smooth_gradients_kernel");
  Kokkos::fence();
  return;
}

/// @brief Loop over all particles and initialize them
/// @details This function loops over all particles and initializes quantities
/// that do not depend on the initial conditions. These are
/// - The lorentz factor \f$\gamma\f$.
/// - The square of the lorentz factor \f$\gamma^2\f$.
/// - The cube of the lorentz factor \f$\gamma^3\f$.
/// - The extensive densities \f$\rho_{\text{spec}}\f$, which are set to 1.
/// - The normalized extensive densities \f$\rho_{\text{sph_mass}}\f$, which are
/// set to the volume assigned to the particle.
/// - The bulk pressure, which is set to 0.
/// - Initial guesses for the temperature and chemical potentials 
/// - Initial guesses for the temperature and chemical potentials
/// (800, 0, 0 ,0) MeV.
/// - The freeze out status, which is set to 4 if the energy density is below
/// @note The volume assigned to each particle is assumed to be the grid step in
/// the initial conditions, as specified in the initial condition header.
/// @note Freeze out 4 means the particle is frozen out. 0 means it is not.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::process_initial_conditions()
{
  formatted_output::report("Processing initial conditions");

  //============================================================================
  // TOTAL NUMBER OF PARTICLES FIXED AFTER THIS POINT
  systemPtr->n_particles = systemPtr->particles.size();
  double t = settingsPtr->t0;
  double t2 = t*t;

  double dA = 1;
  switch (D)
  {
  case 1:
    dA = t*settingsPtr->stepEta;
    break;
  case 2:
    dA = settingsPtr->stepx*settingsPtr->stepy;
    break;
  case 3:
    dA = t*settingsPtr->stepx*settingsPtr->stepy*settingsPtr->stepEta;
    break;
  default:
    std::cerr << "Invalid dimensionality!" << std::endl;
    abort();
    break;
  }

  // fill out initial particle information
  //int TMP_particle_count = 0;
	for (auto & p : systemPtr->particles)
  {

    // Set the rest of particle elements using area element
    double u[D];
    for (int i=0; i<D;++i) u[i] = p.hydro.u(i);
    p.hydro.gamma     = TEOM<D>::gamma_calc(u,t2) ;



    // normalize all extensive densities to 1
		p.extensive.s      = 1.0;
		p.extensive.rhoB   = 1.0;
		p.extensive.rhoS   = 1.0;
		p.extensive.rhoQ   = 1.0;

    // normalization of each density include transverse area element dA
		p.sph_mass.s     = 1.0;
		p.sph_mass.rhoB  = 1.0;
		p.sph_mass.rhoS  = 1.0;
		p.sph_mass.rhoQ  = 1.0;

		p.hydro.bulk      = 0.0;

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
		//if ( p.input.e > systemPtr->efcheck )	// impose freeze-out check for e, not s
		//	p.Freeze = 0;
		//else
		//{
		//	p.Freeze = 4;
		//	systemPtr->number_part_fo++;
		//}
  }

  //formatted_output::detail("particles frozen out: "
  //                         + to_string(systemPtr->number_part_fo) );
  //formatted_output::detail("particles not frozen out: "
  //                         + to_string(systemPtr->n_particles
  //                                      - systemPtr->number_part_fo) );

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


/// @brief Sets the bulk pressure \f$\Pi\f$ to carry the trace of the
/// energy-momentum tensor.
/// @todo The math here is not clear to me, and the documentation is not
/// sufficient. Someone should complement this documentation.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::set_bulk_Pi()
{
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
  if ( settingsPtr->initializing_with_full_Tmunu ){
    double t0 = settingsPtr->t0;
    CREATE_VIEW( device_, systemPtr->cabana_particles);

    auto set_bulk = KOKKOS_LAMBDA( const int is, const int ia )
    {
      double tmunu_trace = device_hydro_scalar.access(is, ia, ccake::hydro_info::tmunu_trace);
      double p = device_thermo.access(is, ia, ccake::thermo_info::p);
      double gamma = device_hydro_scalar.access(is, ia, ccake::hydro_info::gamma);
      double sigma_lab = device_hydro_scalar.access(is, ia, ccake::hydro_info::sigma_lab);
      double sigma  = TEOM<D>::get_LRF(sigma_lab, gamma, t0);
      //device_hydro_scalar.access(is, ia, ccake::hydro_info::extensive_bulk) = (tmunu_trace - p)*sigma;
      device_hydro_scalar.access(is, ia, ccake::hydro_info::extensive_bulk) = 0.0;

    };
    Cabana::simd_parallel_for( simd_policy, set_bulk, "set_bulk_Pi_kernel");
    Kokkos::fence();
  }
  return;
}

/// @brief Calculate the extensive shear used for evolution of the shear tensor.
/// @details This function calculates the extensive shear tensor, which is used
/// for the evolution of the shear tensor. The extensive shear tensor is given by
/// \f[\tilde{pi}^{\mu\nu} = \frac{\pi^\mu\nu}{\sigma} = \gamma_0 tau_0 \frac{\pi^\mu\nu}{\sigma*}\f]
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::calculate_extensive_shv()
{
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
  CREATE_VIEW( device_, systemPtr->cabana_particles);
  double t0 = settingsPtr->t0;
  auto calculate_shear = KOKKOS_LAMBDA( const int is, const int ia )
  {
    double sigma_lab = device_hydro_scalar.access(is, ia, ccake::hydro_info::sigma_lab);
    double u0 = device_hydro_scalar.access(is, ia, ccake::hydro_info::gamma);
    for (int idir=0; idir<2; ++idir)
    for (int jdir=idir; jdir<3; ++jdir)
    {
      double shv = device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, idir+1, jdir+1);
      //shv = 0.;
      device_hydro_shear_aux_vector.access(is, ia, ccake::hydro_info::extensive_shv, idir, jdir) = shv * u0 * t0 / sigma_lab;
      //std::cout << "shv = " << shv << "  extensive_shv = " << device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::extensive_shv, idir, jdir) << std::endl;
    }
  };
  Cabana::simd_parallel_for( simd_policy, calculate_shear, "calculate_extensive_shear_kernel");
  Kokkos::fence();
}


/// @brief Shell function for performing the freeze out procedure.
/// @details This function performs the freeze out procedure, which is
/// implemented in the FreezeOut class. The first step is to check which is to
/// update which particles are to be frozen out by looking at the freeze out
/// energy density, their current energy density and their crrent status. Once
/// the particles to be frozen out are identified, the computation of the freeze
/// out hypersurface is performed.
/// @see FreezeOut::check_freeze_out_status
/// @see FreezeOut::bsqsvfreezeout
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::freeze_out_particles()
{
  //---------------------------------------
  // perform freeze out checks
  int n_freezing_out = 0;
  freezePtr->check_freeze_out_status( n_freezing_out );

  //---------------------------------------
  // finalize frozen out particles
  freezePtr->bsqsvfreezeout( n_freezing_out );
}

/// @brief Shell function to obtain the time derivatives of each quantity.
/// @details This function is a shell function to obtain the time derivatives of
/// each hydrodynamic quantity that needs to be evolved. The steps are the 
/// folloing:
/// - Reset the nearest neighbors list, since the particles moved in the last
/// time step.
/// - Reset the shear viscous tensor to be consistent with all symmetries.
/// - Smooth all particle extensive fields, such as the entropy density etc
/// - Update all particle thermodynamic properties (temperature, pressure etc)
/// - Update transport coefficients, such as the shear and bulk viscosities and
/// the relaxation times.
/// - Compute the gradients necessary for the hydrodynamic evolution.
/// - Invoke the equation of motion class to evaluate the time derivatives.
/// @note The TEOM class, which is a template parameter of this class, should
/// contain a method called `evaluate_time_derivatives` which takes a pointer to
/// the SystemState as an argument.
/// @see SystemState::reset_neighbour_list
/// @see SPHWorkstation::reset_pi_tensor
/// @see SPHWorkstation::smooth_all_particle_fields
/// @see SPHWorkstation::update_all_particle_thermodynamics
/// @see SPHWorkstation::update_all_particle_viscosities
/// @see SPHWorkstation::smooth_all_particle_gradients
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::get_time_derivatives()
{
  Stopwatch sw;
  sw.Start();

  double t = systemPtr->t;
  double t2 = t*t;

  #ifdef DEBUG_SLOW
  systemPtr->copy_device_to_host();
  std::ofstream file;
  file.open("probe.dbg", std::ios::app);
  if (file.is_open()) {
    auto & p =  systemPtr->particles[0];
    file << p.ID << " " << std::setprecision(5) << systemPtr->t << " " << std::setprecision(15) << p.extensive.s << " " << p.r << std::endl;
    file << "\t" << std::setprecision(15) << p.hydro.inside << " " << p.hydro.du_dt << " " << p.hydro.shv(0,0)
    << endl;
    file.close();
  }

  // reset nearest neighbors
  bool fail = false;
  for(auto & p: systemPtr->particles){
    if (p.extensive.s < 0){
        cout << "Negative entropy in particle " << p.ID << " Before reset_neighbour_list " << endl;
        fail = true;
        //cout << p << endl;
    }
  }
  if (fail) exit(8);
  #endif

  // reset nearest neighbors
  systemPtr->reset_neighbour_list();
  // calcuate gamma and velocities
  calculate_gamma_and_velocities();
  // smooth all particle fields - s, rhoB, rhoQ and rhoS and sigma
  smooth_all_particle_fields(t2);
    // Update particle thermodynamic properties
  update_all_particle_thermodynamics();
  // reset pi tensor to be consistent
  // with all essential symmetries
  reset_pi_tensor(t2);
  //add source terms to the energy momentum tensor
  //add_source();
    // update viscosities for all particles
  update_all_particle_viscosities();
    //Computes gradients to obtain dsigma_lab/dt
  smooth_all_particle_gradients(t2);
    //calculate time derivatives needed for equations of motion
  TEOM<D>::evaluate_time_derivatives( systemPtr, settingsPtr );

  // check for causality
  if (settingsPtr->check_causality)
  {
    TEOM<D>::check_causality(systemPtr);
  }

  #ifdef DEBUG_SLOW
  systemPtr->copy_device_to_host();
  file.open("probe2.dbg", std::ios::app);
  if (file.is_open()) {
    auto & p =  systemPtr->particles[0];
    file << p.ID << " " << std::setprecision(5) << systemPtr->t << " " << std::setprecision(15) << p.r ; //Particle info
    file << " " << p.smoothed.s << p.hydro.sigma_lab << " " << p.btrack; //Smoothing info
    file << " " << p.hydro.gradP << " " << p.hydro.gradE << " " << p.hydro.gradBulk << " " << p.hydro.gradV << " " << p.hydro.gradshear << " " << p.hydro.divshear; //gradient info
    file << endl;
    file.close();
  }
  #endif

  sw.Stop();
  formatted_output::update("Finished computing time derivatives in "
                            + to_string(sw.printTime()) + " s.");
  // check/update conserved quantities
  if (settingsPtr->print_conservation_status )
    systemPtr->conservation_energy(false,t);
  return;
}

/// @brief Updates the transport coefficients for all particles.
/// @see transport_coefficients::eta
/// @see transport_coefficients::tau_pi
/// @see transport_coefficients::zeta
/// @see transport_coefficients::tau_Pi
/// @see SPHWorkstation::transp_coeff_params
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::update_all_particle_viscosities()
{
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  auto set_transport_coefficients = KOKKOS_CLASS_LAMBDA(const int is, const int ia ){
    double thermo[thermo_info::NUM_THERMO_INFO];
    for (int ithermo=0; ithermo < thermo_info::NUM_THERMO_INFO; ++ithermo)
      thermo[ithermo] = device_thermo.access(is, ia, ithermo);
    device_hydro_scalar.access(is, ia, ccake::hydro_info::eta_pi)
      = tc::eta(thermo, this->transp_coeff_params );
    device_hydro_scalar.access(is, ia, ccake::hydro_info::tau_pi)
      = tc::tau_pi(thermo,this->transp_coeff_params);
    device_hydro_scalar.access(is, ia, ccake::hydro_info::zeta_Pi)
      = tc::zeta(thermo, this->transp_coeff_params);
    device_hydro_scalar.access(is, ia, ccake::hydro_info::tau_Pi)
      = tc::tau_Pi(thermo, this->transp_coeff_params);
    ///@todo: Read the coefficients instead of hardcoding them
    device_hydro_scalar.access(is, ia, ccake::hydro_info::delta_PiPi)
      = tc::tau_Pi(thermo, this->transp_coeff_params);
    device_hydro_scalar.access(is, ia, ccake::hydro_info::delta_pipi)
      = 4.*tc::tau_pi(thermo, this->transp_coeff_params)/3.;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::lambda_piPi)
      = 0.;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::lambda_Pipi)
      = 0.;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::phi1)
      = 0.;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::phi3)
      = 0.;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::phi6)
      = 0.;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::phi7)
      = 0.;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::tau_pipi)
      = 0.;
  };
  Cabana::simd_parallel_for(simd_policy,set_transport_coefficients, "set_transport_coefficients");
  Kokkos::fence();
  }

/// @brief Update the thermodynamic properties of all particles.
/// @details This function uses the Equation of State to update all the
/// thermodynamic properties of the particles.
/// @todo At the moment, the equation of state functions do not support running
/// in device. Because of this, a copy of the particle data is made to the host
/// and the equation of state is run on the host. Then the data is copied back
/// to the device. This is can be a performance bottleneck, specially for large
/// simulations. The equation of state functions should be modified to run on
/// the device. An interpolator of pre-inverted EoS tables is implemented in
/// experimental status.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::update_all_particle_thermodynamics()
{
  Stopwatch sw;
  sw.Start();
  double t = systemPtr->t;
  double t2 = t*t;
  #ifdef DEBUG_SLOW
  std::cout << "update_all_thermodynamics() before copying back and forth" << std::endl;
  for (auto & p : systemPtr->particles){
    if(p.r(0) > -0.02 && p.r(0) < 0.02){
      std::cout << "eta=  " << p.r(0) << "   ";
      std::cout << "  temp=  " << p.thermo.T << "  thermo.e=   " << p.thermo.e*hbarc_GeVfm << "   input.e=   " << p.input.e*hbarc_GeVfm << \
              "   thermo.s=   " << p.thermo.s << "   input.s=   " << p.input.s << "   " << std::endl;
    }
  }
  #endif
  if(!settingsPtr->online_inverter_enabled){
    eos_interpolatorPtr->fill_thermodynamics(systemPtr->cabana_particles, t);
    #ifdef DEBUG_SLOW
    systemPtr->copy_device_to_host();
    #endif
  }
  else{
  systemPtr->copy_device_to_host();
  for ( auto & p : systemPtr->particles ){
    double s_LRF      = p.thermo.s;
    double rhoB_LRF   = p.thermo.rhoB;
    double rhoS_LRF   = p.thermo.rhoS;
    double rhoQ_LRF   = p.thermo.rhoQ;
    double pos = p.r(0);
    double e_thermo = p.thermo.e;
    double e_input = p.input.e;
    double s_input = p.input.s;
    double T = p.thermo.T;
    try
    {
	    locate_phase_diagram_point_sBSQ( p, s_LRF, rhoB_LRF , rhoS_LRF, rhoQ_LRF );
    }
    catch(const std::exception& e)
    {
      std::cerr << e.what() << '\n';
      cout << "---------------- t = " << systemPtr->t << "fm/c" << endl
           << p;
      systemPtr->print_neighbors(p.ID);
      Kokkos::finalize();
      exit(404);
    }  
  }
  systemPtr->copy_host_to_device();
  }


  sw.Stop();
  formatted_output::detail("Got particle thermodynamics in "
                            + to_string(sw.printTime()) + " s.");
}


/// @brief Shell function to update thermodynamic particle properties.
/// @details This will updated Particle<D>::thermo for the input particle,
/// assuming the smoothed fields \f$\varepsilon\f$, \f$\rho_B\f$, \f$\rho_S\f$, 
/// \f$\rho_Q\f$ are known.
///
/// @todo I am not an expert in the EoS part of the code. Someone, please
/// comnplement this documentation.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
/// @param p The particle for which the thermodynamic properties are to be
/// updated.
/// @param e_In The energy density of the particle.
/// @param rhoB_In The baryon density of the particle.
/// @param rhoS_In The strangeness density of the particle.
/// @param rhoQ_In The electric charge density of the particle.
/// @return The entropy density of the particle.
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

/// @brief Shell function to update thermodynamic particle properties.
/// @details This will updated Particle<D>::thermo for the input particle,
/// assuming the smoothed field \f$\varepsilon\f$is known. The charge densities
/// are assumed to be zero.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
/// @return The entropy density of the particle.
template<unsigned int D, template<unsigned int> class TEOM>
double SPHWorkstation<D, TEOM>::locate_phase_diagram_point_eBSQ(Particle<D> & p, double e_In)
                 { return locate_phase_diagram_point_eBSQ( p, e_In, 0.0, 0.0, 0.0 ); }

/// @brief Shell function to update thermodynamic particle properties.
/// @details This will updated Particle<D>::thermo for the input particle,
/// assuming the smoothed field \f$s\f$ is known. The charge densities
/// are assumed to be zero.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::locate_phase_diagram_point_sBSQ(Particle<D> & p, double s_In) // previously update_s
               { locate_phase_diagram_point_sBSQ(p, s_In, 0.0, 0.0, 0.0 ); }


/// @brief Add a buffer of particles to the simulation.
/// @details This function adds a buffer of particles to the simulation.
/// At the moment, this function is only implemented for the ICCING initial
/// conditions.
///
/// The buffer particles are added with a small energy density, and are placed
/// in the grid points where particles are absent.
/// @todo Details of this function still needs to be documented.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::add_buffer(double default_e)
{
  if ( settingsPtr->IC_type != "ICCING" )
  {
    std::cerr << "WARNING: ADDING A BUFFER IS CURRENTLY ONLY DEFINED FOR ICCING"
                 " INITIAL CONDITIONS.  NO BUFFER ADDED." << std::endl;
    return;
  }

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

/// @brief Decides if the simulation should continue evolving.
/// @todo This function should be aware if we have freeze out turned on or off.
/// @details This function decides if the simulation should continue evolving.
/// The criteria to stop the simulation are:
/// - If all particles have frozen out, if the freeze out flag is set.
/// - If the simulation goes above the maximum time set in the settings.
/// @tparam D Dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
/// @return `true` if the simulation should continue, `false` otherwise.
template<unsigned int D, template<unsigned int> class TEOM>
bool SPHWorkstation<D, TEOM>::continue_evolution()
{
  std::cout << "t = " << systemPtr->t << std::endl;

  bool keep_going =  true;
  if(settingsPtr->particlization_enabled )
    keep_going = (systemPtr->number_part_fo != systemPtr->n_particles); // all particles have frozen out. Break evolution.

  keep_going = keep_going && (systemPtr->t <= settingsPtr->max_tau); // time is up. Break evolution.
  return keep_going;
}

/// @brief Resets advance the simulation one timestep.
/// @details This function advances the simulation one timestep. For each
/// time step, the following operations are performed:
/// - Call the evolver to execute the timestep. The get_time_derivatives
/// function will be passed to the evolver and is its responsability to call it
/// according to the method used.
/// - Perform the freeze out procedure if the freeze out flag is set.
/// - A safeguard is placed to stop the simulation if the maximum number of
/// timesteps is reached.
/// @todo The criteria to stop the simulation seems to be scattered in the code.
/// It should be centralized in a single place. I am of the opinion that it
/// should be in the BSQHydro::run function.
/// @param dt The size of the timestep.
/// @param rk_order The order of the Runge-Kutta method to be used.
/// @tparam D The dimensionality of the simulation.
/// @tparam TEOM The equation of motion class to be used in the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::advance_timestep( double dt, int rk_order )
{
  Stopwatch sw;
  sw.Start();

  //Bulk of code evaluation is done below
  evolver.execute_timestep( dt, rk_order,
                            [this]{ this->get_time_derivatives(); } );

  // Perform freeze out
  if ( settingsPtr->particlization_enabled ) freeze_out_particles();

  // set number of particles which have frozen out
  //systemPtr->number_part = systemPtr->get_frozen_out_count();
  // keep track of how many timesteps have elapsed
  systemPtr->number_of_elapsed_timesteps++;

  sw.Stop();
  formatted_output::detail("finished timestep in "
                            + to_string(sw.printTime()) + " s");
  return;
}


///@brief Function to regulate negative extensive entropy.
///@details This function is responsible for regulating hydrodynamic evolution,
// ensuring that some constraints are always met. Those are
/// - The viscous shear tensor is traceless and orthogonal to the velocity.
/// - The viscous bulk pressure and shear viscous tensor are regulated to
/// prevent the Reynolds number from being too high. A user-specifiable
/// threshold is used to decide when to regulate the viscous terms. The
/// dampening of the viscous terms is done by multiplying them by
/// \f$\frac{\text{Reynolds threshold}}{1+\text{Reynolds threshold}}\f$.
/// - We avoid negative extensive entropy. We only do so if the particle is
/// below the freeze-out energy density. If this criteria is met, we freeze the
/// particle out and set its entropy to zero.
///@tparam D The number of spatial dimensions.
///@param reset_pi A functional that ensures that the viscous shear tensor is
/// traceless and orthogonal to the velocity.
//template <unsigned int D, template<unsigned int> class TEOM>
//void SPHWorkstation<D, TEOM>::regulator(){
//
//  TEOM<D>::reset_pi_tensor(systemPtr);
//
//  if (!settingsPtr->regulate_dissipative_terms) return;
//
//  CREATE_VIEW(device_, systemPtr->cabana_particles);
//  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
//
//  double eFO = systemPtr->efcheck;
//  double t = systemPtr->t;
//  double t4 = t*t*t*t;
//  double t2 = t*t;
//  double Reynolds_threshold = settingsPtr->regulator_threshold;
//  // double dampening = max(1., Reynolds_threshold);
//
//  //Simplified management of freeze-out status if we disabled particlization.
//  if ( !settingsPtr->particlization_enabled ){
//    //Do simplified check of freeze-out if we disabled particlization
//    auto simple_freeze = KOKKOS_LAMBDA(const int is, const int ia){
//      double epsilon = device_thermo.access(is, ia, ccake::thermo_info::e);
//      double freeze = device_freeze.access(is, ia);
//      if (epsilon < eFO && freeze == 0){
//        device_freeze.access(is, ia) = 1;
//      } else if (epsilon > eFO && freeze != 4){
//        device_freeze.access(is, ia) = 0;
//      } else if (epsilon < eFO && freeze == 1){
//        device_freeze.access(is, ia) = 4;
//      }
//    };
//    Cabana::simd_parallel_for(simd_policy, simple_freeze, "simple_freeze");
//    Kokkos::fence();
//  }
//
//  auto regulate_reynolds = KOKKOS_LAMBDA(const int is, const int ia){
//    //Loading variables
//    double p = device_thermo.access(is, ia, ccake::thermo_info::p);
//    double Bulk = device_hydro_scalar.access(is, ia, hydro_info::extensive_bulk);
//    double pi[D+1][D+1];
//    for (int idir=0; idir<D+1; ++idir)
//    for (int jdir=0; jdir<D+1; ++jdir)
//      pi[idir][jdir] = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);
//    double shv33 = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv,3,3);
//
//    double R_shear = pi[0][0]*pi[0][0];
//
//    // //Space diagonal components contributions
//    for (int idir=0; idir<D; ++idir)
//      R_shear += pi[idir+1][idir+1]*pi[idir+1][idir+1];
//    if (D == 3){
//      R_shear += pi[D][D]*pi[D][D]*(t4-1);
//    } else if (D == 2 || D == 1){
//      R_shear += shv33*shv33*t4;
//    }
//
//    //Off-diagonal time components contributions
//    for (int idir=0; idir<D; ++idir)
//      R_shear -= 2.*pi[0][idir+1]*pi[0][idir+1];
//    if (D == 3)
//      R_shear -= 2.*pi[0][3]*pi[0][3]*(t2-1);
//
//    // //Off-diagonal space components contributions
//    for (int idir=0; idir<D; ++idir)
//    for (int jdir=idir+1; jdir<D; ++jdir)
//      R_shear += 2.*pi[idir+1][jdir+1]*pi[idir+1][jdir+1];
//    if (D == 3)
//      for (int idir=0; idir<D; ++idir)
//        R_shear += 2.*pi[idir+1][3]*pi[idir+1][3]*(t2-1);
//
//    R_shear = Kokkos::sqrtf(R_shear)/p;
//    double R_bulk = Kokkos::fabs(Bulk)/p;
//
//    //Regulate if necessary
//    double dampening = 0;
//    if (R_shear > Reynolds_threshold){
//      dampening = Reynolds_threshold/R_shear;
//      if (dampening > 1.){
//        cout << "Dampening does not dampens" << endl;
//        exit(EXIT_FAILURE);
//      }
//      for (int idir=0; idir<D+1; ++idir)
//      for (int jdir=0; jdir<D+1; ++jdir)
//        device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir) =
//          pi[idir][jdir]*dampening;
//      device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv,3,3) = shv33*dampening;
//    }
//
//    if (R_bulk > Reynolds_threshold){
//      dampening = Reynolds_threshold/R_bulk;
//      device_hydro_scalar.access(is, ia, hydro_info::extensive_bulk) = Bulk*dampening;
//    }
//
//  };
//  if(settingsPtr->regulate_dissipative_terms){
//    Cabana::simd_parallel_for(simd_policy, regulate_reynolds, "regulate_reynolds");
//    Kokkos::fence();
//  }
//
//}


template<>
void SPHWorkstation<3, EoM_cartesian>::smooth_all_particle_gradients(double time_squared)
{
  double hT = settingsPtr->hT;
  double t = systemPtr->t;
  bool using_shear = settingsPtr->using_shear;
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  CREATE_VIEW(device_, systemPtr->cabana_particles);
  //Reset gradients - No need to take into account self interaction, since gradient of Kernel for r = 0 is zero
  auto reset_gradients = KOKKOS_LAMBDA(const int is, const int ia){
    for (int idir=0; idir < 3; ++idir)
    {
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradP, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradBulk, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradshear, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::divshear, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradE,idir)  = 0.0;
      for (int jdir=0; jdir < 3; ++jdir)
        device_hydro_space_matrix.access(is, ia,ccake::hydro_info::gradV, idir, jdir) = 0.0;
    }
  };
  Cabana::simd_parallel_for(simd_policy,reset_gradients, "reset_gradients");
  Kokkos::fence();

  ///particle_a is the one that will be updated. Particle b is its neighbors.
  auto smooth_gradients = KOKKOS_LAMBDA(const int particle_a, const int particle_b )
  {
    //Compute kernel gradient
    double rel_sep[3], pos_a[3], pos_b[3]; // cache for relative separation
    for (int idir = 0; idir < 3; ++idir){
      pos_a[idir] = device_position(particle_a,idir);
      pos_b[idir] = device_position(particle_b,idir);
      rel_sep[idir] = pos_a[idir] - pos_b[idir];
    }
    double rel_sep_norm = SPHkernel<3>::distance(pos_a,pos_b);
    double gradK_aux[3];
    ccake::SPHkernel<3>::gradKernel( rel_sep, rel_sep_norm, hT, gradK_aux );

    double pressure_a     = device_thermo(particle_a, thermo_info::p);
    double pressure_b     = device_thermo(particle_b, thermo_info::p);
    double energy_a       = device_thermo(particle_a, thermo_info::e);
    double energy_b       = device_thermo(particle_b, thermo_info::e);
    double sigma_lab_a        = device_hydro_scalar(particle_a, hydro_info::sigma_lab);
    double sigma_lab_b        = device_hydro_scalar(particle_b, hydro_info::sigma_lab);
    double sigsqra        = 1.0/(sigma_lab_a*sigma_lab_a);
    double sigsqrb        = 1.0/(sigma_lab_b*sigma_lab_b);
    double sph_mass_s_b  = device_sph_mass(particle_b, densities_info::s);
    double bulk_a         = device_hydro_scalar(particle_a, hydro_info::bulk);
    double bulk_b         = device_hydro_scalar(particle_b, hydro_info::bulk);
    double gamma_a        = device_hydro_scalar(particle_a, hydro_info::gamma);
    double gamma_b        = device_hydro_scalar(particle_b, hydro_info::gamma);

    cartesian::Vector<double,3> gradK, sigsigK, gradP, gradE, v_a, v_b, gradBulk, gradshear, divshear;
    cartesian::Matrix<double, 3, 3> gradV;
    for(int idir=0; idir<3; ++idir){
      gradK(idir) = gradK_aux[idir];
      v_a(idir) = device_hydro_vector(particle_a, hydro_info::v, idir);
      v_b(idir) = device_hydro_vector(particle_b, hydro_info::v, idir);
    };

    sigsigK = sph_mass_s_b * sigma_lab_a * gradK;
    gradP = (sigsqrb*pressure_b + sigsqra*pressure_a ) * sigsigK;
    gradE = (sph_mass_s_b/sigma_lab_a)*( energy_b-energy_a )*gradK;
    //gradE = (sigsqrb*energy_b + sigsqra*energy_a ) * sigsigK;
    gradBulk = (bulk_b/sigma_lab_b/sigma_lab_b + bulk_a/sigma_lab_a/sigma_lab_a)*sigsigK;
    //gradV    = (sph_mass_s_b/sigma_lab_a)*( v_b -  v_a )*gradK;
    for(int idir=0; idir<3; ++idir)
    for(int jdir=0; jdir<3; ++jdir)
      gradV(idir,jdir) = (sph_mass_s_b/sigma_lab_a)*( v_b(jdir) -  v_a(jdir) )*gradK(idir);


    if (using_shear){
        cartesian::Matrix<double,4,4> shv_a, shv_b;
        for(int idir=0; idir<4;++idir)
        for(int jdir=0; jdir<4;++jdir){
          shv_a(idir,jdir) = device_hydro_spacetime_matrix(particle_a, hydro_info::shv, idir, jdir);
          shv_b(idir,jdir) = device_hydro_spacetime_matrix(particle_b, hydro_info::shv, idir, jdir);
        }
        cartesian::Vector<double,3> shv_0i_a  = cartesian::rowp1(0, shv_a);
        cartesian::Vector<double,3> shv_0i_b  = cartesian::rowp1(0, shv_b);

        cartesian::Matrix<double,3,3> vminia, vminib;
        cartesian::mini(vminia, shv_a);
        cartesian::mini(vminib, shv_b);
        gradshear = cartesian::contract(sigsigK, v_a)*( sigsqrb*shv_0i_b + sigsqra*shv_0i_a );
        divshear  = sigsqrb*sigsigK*cartesian::transpose(vminib)
                  + sigsqra*sigsigK*cartesian::transpose(vminia);
    }

    //Accumulate
    for(int idir=0; idir<3; ++idir) {
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradP, idir), gradP(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradE, idir), gradE(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradBulk, idir), gradBulk(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradshear, idir), gradshear(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::divshear, idir), divshear(idir));
      for(int jdir=0; jdir<3; ++jdir)
        Kokkos::atomic_add(&device_hydro_space_matrix(particle_a, hydro_info::gradV, idir, jdir), gradV(idir, jdir));
    }
  };
  auto range_policy = Kokkos::RangePolicy<ExecutionSpace>(0, systemPtr->cabana_particles.size());
  Cabana::neighbor_parallel_for(range_policy, smooth_gradients,
                                 systemPtr->neighbour_list, Cabana::FirstNeighborsTag(),
                                 Cabana::TeamOpTag(), "smooth_gradients_kernel");
  Kokkos::fence();
  return;
}



template<>
void SPHWorkstation<2, EoM_cartesian>::smooth_all_particle_gradients(double time_squared)
{
  double hT = settingsPtr->hT;
  double t = systemPtr->t;
  bool using_shear = settingsPtr->using_shear;
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  CREATE_VIEW(device_, systemPtr->cabana_particles);
  //Reset gradients - No need to take into account self interaction, since gradient of Kernel for r = 0 is zero
  auto reset_gradients = KOKKOS_LAMBDA(const int is, const int ia){
    for (int idir=0; idir < 2; ++idir)
    {
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradP, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradBulk, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradshear, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::divshear, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradE,idir)  = 0.0;
      for (int jdir=0; jdir < 2; ++jdir)
        device_hydro_space_matrix.access(is, ia,ccake::hydro_info::gradV, idir, jdir) = 0.0;
    }
  };
  Cabana::simd_parallel_for(simd_policy,reset_gradients, "reset_gradients");
  Kokkos::fence();

  ///particle_a is the one that will be updated. Particle b is its neighbors.
  auto smooth_gradients = KOKKOS_LAMBDA(const int particle_a, const int particle_b )
  {
    //Compute kernel gradient
    double rel_sep[2], pos_a[2], pos_b[2]; // cache for relative separation
    for (int idir = 0; idir < 2; ++idir){
      pos_a[idir] = device_position(particle_a,idir);
      pos_b[idir] = device_position(particle_b,idir);
      rel_sep[idir] = pos_a[idir] - pos_b[idir];
    }
    double rel_sep_norm = SPHkernel<2>::distance(pos_a,pos_b);
    double gradK_aux[2];
    ccake::SPHkernel<2>::gradKernel( rel_sep, rel_sep_norm, hT, gradK_aux );

    double pressure_a     = device_thermo(particle_a, thermo_info::p);
    double pressure_b     = device_thermo(particle_b, thermo_info::p);
    double energy_a       = device_thermo(particle_a, thermo_info::e);
    double energy_b       = device_thermo(particle_b, thermo_info::e);
    double sigma_lab_a        = device_hydro_scalar(particle_a, hydro_info::sigma_lab);
    double sigma_lab_b        = device_hydro_scalar(particle_b, hydro_info::sigma_lab);
    double sigsqra        = 1.0/(sigma_lab_a*sigma_lab_a);
    double sigsqrb        = 1.0/(sigma_lab_b*sigma_lab_b);
    double sph_mass_s_b  = device_sph_mass(particle_b, densities_info::s);
    double bulk_a         = device_hydro_scalar(particle_a, hydro_info::bulk);
    double bulk_b         = device_hydro_scalar(particle_b, hydro_info::bulk);
    double gamma_a        = device_hydro_scalar(particle_a, hydro_info::gamma);
    double gamma_b        = device_hydro_scalar(particle_b, hydro_info::gamma);

    cartesian::Vector<double,2> gradK, sigsigK, gradP, gradE, v_a, v_b, gradBulk, gradshear, divshear;
    cartesian::Matrix<double, 2, 2> gradV;
    for(int idir=0; idir<2; ++idir){
      gradK(idir) = gradK_aux[idir];
      v_a(idir) = device_hydro_vector(particle_a, hydro_info::v, idir);
      v_b(idir) = device_hydro_vector(particle_b, hydro_info::v, idir);
    };

    sigsigK = sph_mass_s_b * sigma_lab_a * gradK;
    gradP = (sigsqrb*pressure_b + sigsqra*pressure_a ) * sigsigK;
    gradE = (sph_mass_s_b/sigma_lab_a)*( energy_b-energy_a )*gradK;
    //gradE = (sigsqrb*energy_b + sigsqra*energy_a ) * sigsigK;
    gradBulk = (bulk_b/sigma_lab_b/sigma_lab_b + bulk_a/sigma_lab_a/sigma_lab_a)*sigsigK;
    //gradV    = (sph_mass_s_b/sigma_lab_a)*( v_b -  v_a )*gradK;
    for(int idir=0; idir<2; ++idir)
    for(int jdir=0; jdir<2; ++jdir)
      gradV(idir,jdir) = (sph_mass_s_b/sigma_lab_a)*( v_b(jdir) -  v_a(jdir) )*gradK(idir);


    if (using_shear){
        cartesian::Matrix<double,3,3> shv_a, shv_b;
        for(int idir=0; idir<3;++idir)
        for(int jdir=0; jdir<3;++jdir){
          shv_a(idir,jdir) = device_hydro_spacetime_matrix(particle_a, hydro_info::shv, idir, jdir);
          shv_b(idir,jdir) = device_hydro_spacetime_matrix(particle_b, hydro_info::shv, idir, jdir);
        }
        cartesian::Vector<double,2> shv_0i_a  = cartesian::rowp1(0, shv_a);
        cartesian::Vector<double,2> shv_0i_b  = cartesian::rowp1(0, shv_b);

        cartesian::Matrix<double,2,2> vminia, vminib;
        cartesian::mini(vminia, shv_a);
        cartesian::mini(vminib, shv_b);
        gradshear = cartesian::contract(sigsigK, v_a)*( sigsqrb*shv_0i_b + sigsqra*shv_0i_a );
        divshear  = sigsqrb*sigsigK*cartesian::transpose(vminib)
                  + sigsqra*sigsigK*cartesian::transpose(vminia);
    }

    //Accumulate
    for(int idir=0; idir<2; ++idir) {
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradP, idir), gradP(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradE, idir), gradE(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradBulk, idir), gradBulk(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradshear, idir), gradshear(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::divshear, idir), divshear(idir));
      for(int jdir=0; jdir<2; ++jdir)
        Kokkos::atomic_add(&device_hydro_space_matrix(particle_a, hydro_info::gradV, idir, jdir), gradV(idir, jdir));
    }
  };
  auto range_policy = Kokkos::RangePolicy<ExecutionSpace>(0, systemPtr->cabana_particles.size());
  Cabana::neighbor_parallel_for(range_policy, smooth_gradients,
                                 systemPtr->neighbour_list, Cabana::FirstNeighborsTag(),
                                 Cabana::TeamOpTag(), "smooth_gradients_kernel");
  Kokkos::fence();
  return;
}



template<>
void SPHWorkstation<1, EoM_cartesian>::smooth_all_particle_gradients(double time_squared)
{
  double hT = settingsPtr->hT;
  double t = systemPtr->t;
  bool using_shear = settingsPtr->using_shear;
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  CREATE_VIEW(device_, systemPtr->cabana_particles);
  //Reset gradients - No need to take into account self interaction, since gradient of Kernel for r = 0 is zero
  auto reset_gradients = KOKKOS_LAMBDA(const int is, const int ia){
    for (int idir=0; idir < 1; ++idir)
    {
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradP, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradBulk, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradshear, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::divshear, idir) = 0.0;
      device_hydro_vector.access(is, ia, ccake::hydro_info::gradE,idir)  = 0.0;
      for (int jdir=0; jdir < 1; ++jdir)
        device_hydro_space_matrix.access(is, ia,ccake::hydro_info::gradV, idir, jdir) = 0.0;
    }
  };
  Cabana::simd_parallel_for(simd_policy,reset_gradients, "reset_gradients");
  Kokkos::fence();

  ///particle_a is the one that will be updated. Particle b is its neighbors.
  auto smooth_gradients = KOKKOS_LAMBDA(const int particle_a, const int particle_b )
  {
    //Compute kernel gradient
    double rel_sep[1], pos_a[1], pos_b[1]; // cache for relative separation
    for (int idir = 0; idir < 1; ++idir){
      pos_a[idir] = device_position(particle_a,idir);
      pos_b[idir] = device_position(particle_b,idir);
      rel_sep[idir] = pos_a[idir] - pos_b[idir];
    }
    double rel_sep_norm = SPHkernel<1>::distance(pos_a,pos_b);
    double gradK_aux[1];
    ccake::SPHkernel<1>::gradKernel( rel_sep, rel_sep_norm, hT, gradK_aux );

    double pressure_a     = device_thermo(particle_a, thermo_info::p);
    double pressure_b     = device_thermo(particle_b, thermo_info::p);
    double energy_a       = device_thermo(particle_a, thermo_info::e);
    double energy_b       = device_thermo(particle_b, thermo_info::e);
    double sigma_lab_a        = device_hydro_scalar(particle_a, hydro_info::sigma_lab);
    double sigma_lab_b        = device_hydro_scalar(particle_b, hydro_info::sigma_lab);
    double sigsqra        = 1.0/(sigma_lab_a*sigma_lab_a);
    double sigsqrb        = 1.0/(sigma_lab_b*sigma_lab_b);
    double sph_mass_s_b  = device_sph_mass(particle_b, densities_info::s);
    double bulk_a         = device_hydro_scalar(particle_a, hydro_info::bulk);
    double bulk_b         = device_hydro_scalar(particle_b, hydro_info::bulk);
    double gamma_a        = device_hydro_scalar(particle_a, hydro_info::gamma);
    double gamma_b        = device_hydro_scalar(particle_b, hydro_info::gamma);

    cartesian::Vector<double,1> gradK, sigsigK, gradP, gradE, v_a, v_b, gradBulk, gradshear, divshear;
    cartesian::Matrix<double, 1, 1> gradV;
    for(int idir=0; idir<1; ++idir){
      gradK(idir) = gradK_aux[idir];
      v_a(idir) = device_hydro_vector(particle_a, hydro_info::v, idir);
      v_b(idir) = device_hydro_vector(particle_b, hydro_info::v, idir);
    };

    sigsigK = sph_mass_s_b * sigma_lab_a * gradK;
    gradP = (sigsqrb*pressure_b + sigsqra*pressure_a ) * sigsigK;
    gradE = (sph_mass_s_b/sigma_lab_a)*( energy_b-energy_a )*gradK;
    //gradE = (sigsqrb*energy_b + sigsqra*energy_a ) * sigsigK;
    gradBulk = (bulk_b/sigma_lab_b/sigma_lab_b + bulk_a/sigma_lab_a/sigma_lab_a)*sigsigK;
    //gradV    = (sph_mass_s_b/sigma_lab_a)*( v_b -  v_a )*gradK;
    for(int idir=0; idir<1; ++idir)
    for(int jdir=0; jdir<1; ++jdir)
      gradV(idir,jdir) = (sph_mass_s_b/sigma_lab_a)*( v_b(jdir) -  v_a(jdir) )*gradK(idir);


    if (using_shear){
        double shv_0i_a  = device_hydro_spacetime_matrix(particle_a, hydro_info::shv, 0, 3);
        double shv_0i_b  = device_hydro_spacetime_matrix(particle_b, hydro_info::shv, 0, 3);
        cartesian::Matrix<double,1,1> vminia, vminib;
        vminia(0,0) = device_hydro_spacetime_matrix(particle_a, hydro_info::shv, 3, 3);
        vminib(0,0) = device_hydro_spacetime_matrix(particle_b, hydro_info::shv, 3, 3);
        gradshear = cartesian::contract(sigsigK, v_a)*( sigsqrb*shv_0i_b + sigsqra*shv_0i_a );
        ///\todo do we need this transpose here?
        divshear  = sigsqrb*sigsigK*vminib
                  + sigsqra*sigsigK*vminia;  
    }

    //Accumulate
    for(int idir=0; idir<1; ++idir) {
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradP, idir), gradP(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradE, idir), gradE(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradBulk, idir), gradBulk(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::gradshear, idir), gradshear(idir));
      Kokkos::atomic_add(&device_hydro_vector(particle_a, hydro_info::divshear, idir), divshear(idir));
      for(int jdir=0; jdir<1; ++jdir)
        Kokkos::atomic_add(&device_hydro_space_matrix(particle_a, hydro_info::gradV, idir, jdir), gradV(idir, jdir));
    }
  };
  auto range_policy = Kokkos::RangePolicy<ExecutionSpace>(0, systemPtr->cabana_particles.size());
  Cabana::neighbor_parallel_for(range_policy, smooth_gradients,
                                 systemPtr->neighbour_list, Cabana::FirstNeighborsTag(),
                                 Cabana::TeamOpTag(), "smooth_gradients_kernel");
  Kokkos::fence();
  return;
}

}; //namespace ccake
