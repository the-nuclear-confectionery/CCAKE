#include "sph_workstation.h"

using namespace constants;
using namespace ccake;

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
  // set up equation of motion object
  EoMPtr = std::make_shared<TEOM<D>>(settingsPtr);

  //----------------------------------------
  // set up system state
  systemPtr->initialize();

  //----------------------------------------
  // set up equation of state
  eos.set_SettingsPtr( settingsPtr );
  eos.init();

  //----------------------------------------
  // set up transport coefficients
  transport_coefficients.set_SettingsPtr( settingsPtr );
  transport_coefficients.initialize( settingsPtr->etaMode, settingsPtr->shearRelaxMode,
                 settingsPtr->zetaMode, settingsPtr->bulkRelaxMode );

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
  evolver.initialize( settingsPtr, systemPtr );

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
  auto reset_pi_tensor_lambda  = KOKKOS_LAMBDA(const unsigned int iparticle)
  {
    double u[D]; // Cache for fluid velocity
    double pi_space[D][D]; //Cache for pi^{ij}
    double pi_time_vector[D]; //Cache for the components pi^{0i}
    double pi00; //Cache for pi^{00}
    double pi_diag[D+1]; //Cache for the diagonal components of pi tensor

    //Retrieve information from the particle to local variables
    for( unsigned int idir=0; idir<D; idir++ ){
      u[idir] = device_hydro_vector(iparticle, ccake::hydro_info::u, idir);
      pi_diag[idir+1] = device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, idir+1, idir+1);
      for( unsigned int jdir=0; jdir<D; jdir++ )
        pi_space[idir][jdir] = device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, idir+1, jdir+1);
    }
    
    //Computes the gamma factor
    double u0 = EoMPtr->gamma_calc(u,time_squared);
    
    //Computes pi^{0i} = u_j pi^{ij}/gamma (obtained from the requirement u_\mu pi^{\mu i} = 0)
    for( unsigned int idir=0; idir<D; idir++ )
      pi_time_vector[idir] = EoMPtr->dot(pi_space[idir],u,time_squared)/u0;

    //Set the shear tensor pi^{0 0} = pi^{0j} u_j /gamma (from the requirement u_\mu pi^{\mu 0} = 0)
    //This is equivalent to pi^{i j} u^i u^j /gamma^2
    pi00 = EoMPtr->dot(pi_time_vector,u,time_squared)/u0;
    pi_diag[0] = pi00;

    //Lastly, we need to ensure that the tensor is traceless. We do this by setting
    //the value of the last component pi^{33} = pi^{00}-pi^{11}-pi^{22}
    double pi33 = EoMPtr->get_shvDD(pi_diag, time_squared);
    if( D != 2 )
      pi_diag[D] = pi33;

    
    //Set the values of the pi tensor
    device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, 0, 0) = pi_diag[0]; //pi^{00}
    device_hydro_scalar(iparticle, ccake::hydro_info::shv33) = pi33;  //pi^{33}
    for( unsigned int idir=1; idir<D+1; idir++ ){
      device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, 0, idir) = pi_time_vector[idir-1]; //pi^{0i}
      device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, idir, 0) = pi_time_vector[idir-1]; //pi^{i0}
      device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, idir, idir) = pi_diag[idir]; //pi^{ii}
      for (unsigned int jdir=idir+1; jdir<D+1; jdir++){
        device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, idir, jdir) = pi_space[idir-1][jdir-1]; //pi^{ij}
        device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, jdir, idir) = pi_space[idir-1][jdir-1]; //pi^{ji}
      }
    }
    //Set values of pi^{ij} and u^i u^j tensors
    for( unsigned int idir=0; idir<D; idir++ )
    for( unsigned int jdir=0; jdir<D; jdir++ )
    {
      device_hydro_space_matrix(iparticle, ccake::hydro_info::pimin, idir, jdir) = pi_space[idir][jdir]; //pi^{ij}
      device_hydro_space_matrix(iparticle, ccake::hydro_info::uu, jdir, idir) = pi_space[idir][jdir]; //pi^{ji}
    }
    //Update gamma factor
    device_hydro_scalar(iparticle, ccake::hydro_info::gamma) = u0;
  };


  Kokkos::parallel_for("reset_pi_tensor", systemPtr->n_particles, reset_pi_tensor_lambda);

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

  #ifndef DEBUG
  #error "STOP RIGHT THERE!!! THIS CODE IS USING CONFORMAL EOS! It is not ready for production yet!"
  #error "This is hard coded and was done to facilitate parallelism implementation."
  #error "If you want to use it in production with realistic EOS, you need to implement calls to the EoS properly."
  #endif

  double t0 = settingsPtr->t0;
  // loop over all particles
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto init_particles = KOKKOS_LAMBDA(const int iparticle)
  {
    double e_input = device_input(iparticle, densities_info::e);
    double rhoB_input = device_input(iparticle, densities_info::rhoB);
    double rhoS_input = device_input(iparticle, densities_info::rhoS);
    double rhoQ_input = device_input(iparticle, densities_info::rhoQ);

    ///------------------------------------------------------------------------------------------
    ///TODO: implement calls to EOS here. Using conformal EOS for now.
    // We are forcibly using a conformal EOS for now. A proper EOS will be implemented later.
    // Example on what this should look like:
    // ```
    // double s_input = locate_phase_diagram_point_eBSQ(e_input, rhoB_input, rhoS_input, rhoQ_input);
    // ```
    // Now, update the particle's thermodynamic quantities

    double static const C = 0.2281360133; // See Aguiar et al. (2001) for details
    double sVal = std::pow(e_input/(3*C),3./4.);
    double TVal = 4*C*std::pow(sVal,1./3.);
    double wVal = 4*e_input/3.;
    double pVal = e_input/3.;
    double cs2Val = 1./3.;
    double dwdsVal = 4*C*TVal;
    double AVal = wVal-sVal*dwdsVal;
    
    device_input(iparticle, ccake::densities_info::s ) =   sVal;
    device_thermo(iparticle, ccake::densities_info::s ) =   sVal;
    device_thermo(iparticle, ccake::densities_info::e ) =   e_input;
    device_thermo(iparticle, ccake::thermo_info::T ) =   TVal;
    device_thermo(iparticle, ccake::thermo_info::muB ) =   0;
    device_thermo(iparticle, ccake::thermo_info::muS ) =   0;
    device_thermo(iparticle, ccake::thermo_info::muQ ) =   0;
    device_thermo(iparticle, ccake::thermo_info::s ) =   sVal;
    device_thermo(iparticle, ccake::thermo_info::p ) =   pVal;
    device_thermo(iparticle, ccake::thermo_info::cs2 )  = cs2Val;
    device_thermo(iparticle, ccake::thermo_info::w ) =   wVal;
    device_thermo(iparticle, ccake::thermo_info::A ) =   AVal;
    device_thermo(iparticle, ccake::thermo_info::dwds ) = dwdsVal;
    device_thermo(iparticle, ccake::thermo_info::dwdB ) =   0;
    device_thermo(iparticle, ccake::thermo_info::dwdS ) =   0;
    device_thermo(iparticle, ccake::thermo_info::dwdQ ) =   0;

    double gamma = device_hydro_scalar(iparticle, ccake::hydro_info::gamma);

    device_norm_spec(iparticle, ccake::densities_info::e )    *= e_input*gamma*t0;    // constant after this
    device_norm_spec(iparticle, ccake::densities_info::s )    *= sVal*gamma*t0;       // constant after this
    device_norm_spec(iparticle, ccake::densities_info::rhoB ) *= rhoB_input*gamma*t0; // constant after this
    device_norm_spec(iparticle, ccake::densities_info::rhoS ) *= rhoS_input*gamma*t0; // constant after this
    device_norm_spec(iparticle, ccake::densities_info::rhoQ ) *= rhoQ_input*gamma*t0; // constant after this

    if (sVal < 0.0)
      device_freeze(iparticle) = 4;
    
    ///------------------------------------------------------------------------------------------
  };

  Kokkos::parallel_for("init_particles", systemPtr->n_particles, init_particles);


	swTotal.Stop();
  formatted_output::update("finished initializing particle densities in "
                              + to_string(swTotal.printTime()) + " s");

}

//==============================================================================
template<unsigned int D,  template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::initial_smoothing()
{
  double t_squared = pow(settingsPtr->t0,2);

  //Creates the list of neighbours for each particle
  systemPtr->reset_neighbour_list();

  //Computes not independent components of the shear tensor
  reset_pi_tensor(t_squared);

  // smooth fields over particles
  smooth_all_particle_fields(t_squared);
  
}

//==============================================================================

template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D, TEOM>::smooth_all_particle_fields(double time_squared)
{
  Stopwatch sw;
  sw.Start();

  CREATE_VIEW(device_, systemPtr->cabana_particles);

  //Reset smoothed fields
  auto reset_fields = KOKKOS_CLASS_LAMBDA(const int iparticle)
  {
    device_smoothed(iparticle, ccake::densities_info::s) = 0.0;
    device_smoothed(iparticle, ccake::densities_info::rhoB) = 0.0;
    device_smoothed(iparticle, ccake::densities_info::rhoS) = 0.0;
    device_smoothed(iparticle, ccake::densities_info::rhoQ) = 0.0;
    device_hydro_scalar(iparticle, ccake::hydro_info::sigma) = 0.0;
    
  };
  Kokkos::parallel_for("reset_fields", systemPtr->n_particles, reset_fields);
  Kokkos::fence();

  auto smooth_fields = KOKKOS_LAMBDA(const int iparticle, const int jparticle ){
    
    double r1[D] ,r2[D]; // cache for positions of particles 1 and 2
    for (int idir = 0; idir < D; ++idir){
      r1[idir] = device_position(iparticle,idir);
      r2[idir] = device_position(jparticle,idir);
    }
    double distance = EoMPtr->get_distance(r1,r2,time_squared);
    double kern = kernel::kernel(distance,settingsPtr->hT);

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
     
  Cabana::neighbor_parallel_for( systemPtr->range_policy, smooth_fields, systemPtr->neighbour_list, Cabana::FirstNeighborsTag(),
                                   Cabana::TeamOpTag(), "smooth fields");
  Kokkos::fence();

  double t0 = sqrt(time_squared);
  auto update_densities_lab = KOKKOS_LAMBDA(int iparticle){
    double u0 = device_hydro_scalar(iparticle, ccake::hydro_info::gamma);
    double smoothed_s_lab    = device_hydro_scalar(iparticle, ccake::hydro_info::sigma)/u0/t0;
    double smoothed_rhoB_lab = device_smoothed(iparticle, ccake::densities_info::rhoB)/u0/t0;
    double smoothed_rhoS_lab = device_smoothed(iparticle, ccake::densities_info::rhoS)/u0/t0;
    double smoothed_rhoQ_lab = device_smoothed(iparticle, ccake::densities_info::rhoQ)/u0/t0;

    // UNCOMMENT THIS AND DOCUMENT OUTPUT AS REFERENCE
    ///\todo: For now, we can't use the eos class on the device, so we are going to use conformal eos
    //locate_phase_diagram_point_sBSQ( p, smoothed_s_lab, smoothed_rhoB_lab,
    //                                  smoothed_rhoS_lab, smoothed_rhoQ_lab );

    double static const C = 0.2281360133; // See Aguiar et al. (2001) for details
    double eVal = 3*C*std::pow(smoothed_s_lab,4./3.);
    double sVal = smoothed_s_lab;
    double TVal = 4*C*std::pow(sVal,1./3.);
    double wVal = 4*eVal/3.;
    double pVal = eVal/3.;
    double cs2Val = 1./3.;
    double dwdsVal = 4*C*TVal;
    double AVal = wVal-sVal*dwdsVal;
    
    device_thermo(iparticle, ccake::densities_info::s ) =   sVal;
    device_thermo(iparticle, ccake::densities_info::e ) =   eVal;
    device_thermo(iparticle, ccake::thermo_info::T ) =   TVal;
    device_thermo(iparticle, ccake::thermo_info::muB ) =   0;
    device_thermo(iparticle, ccake::thermo_info::muS ) =   0;
    device_thermo(iparticle, ccake::thermo_info::muQ ) =   0;
    device_thermo(iparticle, ccake::thermo_info::s ) =   sVal;
    device_thermo(iparticle, ccake::thermo_info::p ) =   pVal;
    device_thermo(iparticle, ccake::thermo_info::cs2 )  = cs2Val;
    device_thermo(iparticle, ccake::thermo_info::w ) =   wVal;
    device_thermo(iparticle, ccake::thermo_info::A ) =   AVal;
    device_thermo(iparticle, ccake::thermo_info::dwds ) = dwdsVal;
    device_thermo(iparticle, ccake::thermo_info::dwdB ) =   0;
    device_thermo(iparticle, ccake::thermo_info::dwdS ) =   0;
    device_thermo(iparticle, ccake::thermo_info::dwdQ ) =   0;

    ///\todo: Freeze out will need to be done on device and this is a beast of its own
	  //fo.check_freeze_out_status(p, settingsPtr->t0, count1, systemPtr->n_particles); 

  };

  Kokkos::parallel_for("update densities lab", systemPtr->n_particles, update_densities_lab);

  sw.Stop();
  formatted_output::update("finished smoothing particle fields in "
                            + to_string(sw.printTime()) + " s.");


                  
}

//------------------------------------------------------------------------------




///////////////////////////////////////////////////////////////////////////////
// smoothing routines: first smoothing covers all hydrodyanmical fields
/*

///////////////////////////////////////////////////////////////////////////////
//Second smoothing smoothes the gradients after constructing all the fields 
//and derivatives using the equation of state
template<unsigned int D>
void SPHWorkstation<D>::smooth_gradients( Particle & pa, double tin )
{
  int a = pa.ID;

  auto & pah = pa.hydro;

  pah.gradP     = 0.0;
  pah.gradBulk  = 0.0;
  pah.gradV     = 0.0;
  pah.gradshear = 0.0;
  pah.divshear  = 0.0;

  if ( pa.btrack != -1 ) pa.btrack = 0;

  double rdis = 0;

  auto & a_neighbors = systemPtr->linklist.all_neighbors[a];

  for ( int b : a_neighbors )
  {
    auto & pb                = systemPtr->particles[b];
    auto & pbh               = pb.hydro;

    Vector<double,2> rel_sep = pa.r - pb.r;
    double rel_sep_norm      = Norm( rel_sep );
    Vector<double,2> gradK   = kernel::gradKernel( rel_sep, rel_sep_norm, settingsPtr->hT );
    Vector<double,2> va      = rowp1(0, pah.shv);
    Vector<double,2> vb      = rowp1(0, pbh.shv);
    Matrix<double,2,2> vminia, vminib;
    mini(vminia, pah.shv);
    mini(vminib, pbh.shv);

    double sigsqra           = 1.0/(pah.sigma*pah.sigma);
    double sigsqrb           = 1.0/(pbh.sigma*pbh.sigma);
    Vector<double,2> sigsigK = pb.norm_spec.s * pah.sigma * gradK;

    pah.gradP                += ( sigsqrb*pb.p() + sigsqra*pa.p() ) * sigsigK;

    //===============
    // print status
    if ( VERBOSE > 2 && pa.print_this_particle )
      std::cout << "CHECK grads: " << tin << "   "
                << pah.gradP << "   " << a << "   " << b << "   "
                << sigsqra << "   " << sigsqrb
                << "   " << pa.p() << "   " << pb.p()
                << "   " << pa.get_current_eos_name()
                << "   " << pb.get_current_eos_name()
                << "   " << gradK << "   " << sigsigK
                << "   " << pah.sigma << "\n";

    double relative_distance_by_h = rel_sep_norm / settingsPtr->hT;
    if ( ( relative_distance_by_h <= 2.0 ) && ( a != b ) )
    {
      if ( pa.btrack != -1 ) pa.btrack++;                   // effectively counts nearest neighbors
      if ( pa.btrack ==  1 ) rdis = relative_distance_by_h;
    }

    pah.gradBulk             += ( pbh.Bulk/pbh.sigma/pbh.gamma
                                + pah.Bulk/pah.sigma/pah.gamma)/tin*sigsigK;
    pah.gradV                += (pb.norm_spec.s/pah.sigma)*( pbh.v -  pah.v )*gradK;

    //===============
    // print status
    if ( VERBOSE > 2 && pa.print_this_particle )
        std::cout << "CHECK gradV: " << tin << "   " << a << "   " << b << "   "
                  << pb.norm_spec.s/pah.sigma << "   " << pbh.v -  pah.v
                  << "   " << gradK << "   " << pah.gradV << "\n";

    //===============
    // add shear terms
    if ( settingsPtr->using_shear )
    {
      pah.gradshear            += inner(sigsigK, pah.v)*( sigsqrb*vb + sigsqra*va );
      pah.divshear             += sigsqrb*sigsigK*transpose(vminib)
                                  + sigsqra*sigsigK*transpose(vminia);
    }

    //===============
    // check for nan pressure gradients
    if ( isnan( pah.gradP(0) ) )
    {
      cout << "gradP stopped working" << endl;
      cout << systemPtr->t <<" "  << pah.gradP << " " << a << " " << b << endl;
      cout << pb.norm_spec.s << " " << pah.sigma << " " << pb.p() << endl;
      cout << systemPtr->linklist.Size << " " << pb.s() << " " << pa.s() << endl;

      cout << pa.r << endl;
      cout << pb.r << endl;
      cout << kernel::kernel( pa.r - pb.r, settingsPtr->hT ) << endl;
    }
    else if ( isnan( pah.gradP(1) ) )
      cout << "1 " << gradPressure_weight(a, b)
           << " " << a << " " << b << endl;
  }

  const double hc = constants::hbarc_MeVfm;

  if ( ( pa.btrack == 1 )                               // if particle a has only
            && ( ( pa.T()*hc ) >= 150 ) )               // one nearest neighbor (and T>150),
    freeze_out.frz2[pa.ID].t=tin;                               // set penultimate timestep;
  else if ( ( pa.btrack == 0 )                          // otherwise, if a has no nearest neighbors
            && ( ( pa.T()*hc ) >= 150 )                 // but has T>150 and isn't frozen out,
            && ( pa.Freeze < 4 ) )                      // just print this warning message
    cout << boolalpha << "Missed " << a << " "
         << tin << " " << pa.T()*hc << " "
         << rdis << " " << systemPtr->do_freeze_out << endl;

  return;
}
*/

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
    p.hydro.g2        = p.hydro.gamma*p.hydro.gamma;
    p.hydro.g3        = p.hydro.g2*p.hydro.gamma;


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

    auto set_bulk = KOKKOS_LAMBDA( const int iparticle )
    {
      double varsigma = device_hydro_scalar(iparticle, ccake::hydro_info::varsigma);
      double p = device_thermo(iparticle, ccake::thermo_info::p);
      double u0 = device_hydro_scalar(iparticle, ccake::hydro_info::gamma);
      double sigma = device_hydro_scalar(iparticle, ccake::hydro_info::sigma);

      device_hydro_scalar(iparticle, ccake::hydro_info::Bulk) = (varsigma - p)
                      * u0 * t0 / sigma;
    };

    Kokkos::parallel_for( "set_bulk_Pi", systemPtr->n_particles, set_bulk );
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




////////////////////////////////////////////////////////////////////////////////
template<unsigned int D>
void SPHWorkstation<D>::get_time_derivatives()
{
  // reset nearest neighbors
  reset_linklist();

  // reset pi tensor to be consistent
  // with all essential symmetries
  reset_pi_tensor();

  // smooth all particle fields
  smooth_all_particle_fields();

  // update gamma, velocity, and thermodynamics for all particles
  update_all_particle_thermodynamics();

  // update viscosities for all particles
  update_all_particle_viscosities();

  //Computes gradients to obtain dsigma/dt
  smooth_all_particle_gradients();

  //calculate time derivatives needed for equations of motion
  evaluate_all_particle_time_derivatives();

  // identify and handle particles which have frozen out
  if ( systemPtr->do_freeze_out )
    freeze_out_particles();

  // check/update conserved quantities
  systemPtr->conservation_energy();

  return;
}
*/
////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> class TEOM>
double SPHWorkstation<D,TEOM>::locate_phase_diagram_point_eBSQ( 
                 double e_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{

  #ifndef DEBUG
  #error "STOP RIGHT THERE!!! THIS CODE IS USING CONFORMAL EOS! It is not ready for production yet!"
  #error "This is hard coded and was done to facilitate parallelism implementation."
  #error "If you want to use it in production with realistic EOS, you need to implement calls to the EoS properly."
  #endif

  double static const C = 3*0.2281360133; // See Aguiar et al. (2001) for details
  double sVal = std::pow(e_In/C,3./4.);
  

  return sVal;
}

////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> class TEOM>
double SPHWorkstation<D,TEOM>::locate_phase_diagram_point_eBSQ(double e_In)
                 { return locate_phase_diagram_point_eBSQ( e_In, 0.0, 0.0, 0.0 ); }

////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::locate_phase_diagram_point_sBSQ( Particle<D> & p,
                 double s_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
//  cout << "Rootfinder for p.ID = " << p.ID << endl;

  // default: use particle's current location as initial guess
  eos.tbqs( p.T(), p.muB(), p.muQ(), p.muS(), p.get_current_eos_name() );

  bool update_s_success = eos.update_s( s_In, rhoB_In, rhoS_In, rhoQ_In,
                                        p.print_this_particle );

  if ( update_s_success )
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
    //*/

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

////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> class TEOM>
void SPHWorkstation<D,TEOM>::locate_phase_diagram_point_sBSQ(Particle<D> & p, double s_In) // previously update_s
               { locate_phase_diagram_point_sBSQ(p, s_In, 0.0, 0.0, 0.0 ); }

/*

////////////////////////////////////////////////////////////////////////////////
//  Computes gamma and velocity
template<unsigned int D>
void SPHWorkstation<D>::set_phase_diagram_point(Particle & p)
{
  p.hydro.gamma     = p.gamcalc();
  p.hydro.v         = (1.0/p.hydro.gamma)*p.hydro.u;
  double s_lab      = p.smoothed.s/p.hydro.gamma/systemPtr->t;
  double rhoB_lab   = p.smoothed.rhoB/p.hydro.gamma/systemPtr->t;
  double rhoS_lab   = p.smoothed.rhoS/p.hydro.gamma/systemPtr->t;
  double rhoQ_lab   = p.smoothed.rhoQ/p.hydro.gamma/systemPtr->t;
	locate_phase_diagram_point_sBSQ( p, s_lab, rhoB_lab, rhoS_lab, rhoQ_lab );
}




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






////////////////////////////////////////////////////////////////////////////////
template<unsigned int D>
void SPHWorkstation<D>::set_transport_coefficients( Particle & p )
{
  transport_coefficients.setTherm( p.thermo );
  p.hydro.setas     = transport_coefficients.eta();
  p.hydro.stauRelax = transport_coefficients.tau_pi();
  p.hydro.zeta      = transport_coefficients.zeta();
  p.hydro.tauRelax  = transport_coefficients.tau_Pi();
//  if ( p.ID == 0 )
//  {
//    cout << "check thermo: "
//          << p.hydro.setas << "   "
//          << p.hydro.stauRelax << "   "
//          << p.hydro.zeta << "   "
//          << p.hydro.tauRelax << "   "
//          << p.thermo.cs2 << endl;
//  }
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