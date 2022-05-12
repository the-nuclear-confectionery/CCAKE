#include <algorithm>
#include <memory>

#include "constants.h"
#include "sph_workstation.h"
#include "Stopwatch.h"

using namespace constants;


void SPHWorkstation::set_SystemStatePtr( SystemState * systemPtr_in )
{
  systemPtr = systemPtr_in;
}

void SPHWorkstation::set_SettingsPtr( Settings * settingsPtr_in )
{
  settingsPtr = settingsPtr_in;
}

void SPHWorkstation::reset_pi_tensor()
{
  for ( auto & p : systemPtr->particles )
    p.reset_pi_tensor(systemPtr->t*systemPtr->t);
}


////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::initialize_entropy_and_charge_densities() // formerly updateIC
{
	Stopwatch sw, swTotal;
	swTotal.Start();
	long long failCounter = 0;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;

    systemPtr->n_particles = systemPtr->particles.size();
    cout << "systemPtr->n_particles = " << systemPtr->n_particles << "\n";

    for (int i=0; i<systemPtr->n_particles; i++)
    {
	cout << "----------------------------------------"
			"----------------------------------------" << "\n";


      auto & p = systemPtr->particles[i];


		
		{
			sw.Start();
			cout << setprecision(12) << "Doing this particle: "
					<< p.r(0) << "   " << p.r(1) << "\n";

      // solve for the entropy density
			p.input.s = locate_phase_diagram_point_eBSQ( p,
                    p.input.e, p.input.rhoB, p.input.rhoS, p.input.rhoQ );

			sw.Stop();
			string successString = (p.input.s < 0.0) ?
									"unsuccessfully" : "successfully";
      cout << "Print particle info:\n";
			cout << "    SPH particle " << i << ", locate_phase_diagram_point_eBSQ: completed "
					<< successString << " in " << sw.printTime() << "s." << "\n";

		}

		///////////////////////////////////////////////////////////////////////////
		// for now, if we failed to find a real entropy density for this
		// point, just freeze it out, set its entropy to the freeze-out value,
		// and continue without setting anything else
		if (p.input.s < 0.0)
		{
			// freeze this particle out!
      cout << "This shouldn't have happened!" << endl;
      exit(8);
			p.Freeze = 5;
			//p.Freeze = 4;
			//systemPtr->number_part++;
			////////////////////////////////////////////////////////
		}
		else
		{
			std::cout << "\t --> Solution info: "
                << p.r(0) << "   " << p.r(1) << "   "
                << p.get_current_eos_name() << "\n";
			std::cout << "\t\t - phase diagram point (TBQS): "
                << p.T()*hbarc_MeVfm << "   "
                << p.muB()*hbarc_MeVfm << "   "
                << p.muQ()*hbarc_MeVfm << "   "
                << p.muS()*hbarc_MeVfm << "\n";
			std::cout << "\t\t - input densities (eBSQ):    "
                << p.input.e*hbarc_MeVfm << "   "
                << p.input.rhoB << "   "
                << p.input.rhoS << "   "
                << p.input.rhoQ << "\n";
			cout << "\t\t - output densities (eBSQ):    "
                << p.e()*hbarc_MeVfm << "   "
                << p.rhoB() << "   "
                << p.rhoS() << "   "
                << p.rhoQ() << "\n";
			cout << "\t\t - freeze-out status:   " << p.Freeze << "\n";
		}

    p.hydro.gamma = p.gamcalc();

    p.sigmaweight *= p.input.s*p.hydro.gamma*settingsPtr->t0;	  // sigmaweight is constant after this
    p.rhoB_weight *= p.hydro.gamma*settingsPtr->t0; // rhoB_weight is constant after this
    p.rhoS_weight *= p.hydro.gamma*settingsPtr->t0; // rhoS_weight is constant after this
    p.rhoQ_weight *= p.hydro.gamma*settingsPtr->t0; // rhoQ_weight is constant after this

		p.B *= p.hydro.gamma*settingsPtr->t0;	// B does not evolve in ideal case
		p.S *= p.hydro.gamma*settingsPtr->t0;	// S does not evolve in ideal case
		p.Q *= p.hydro.gamma*settingsPtr->t0;	// Q does not evolve in ideal case

	cout << "----------------------------------------"
			"----------------------------------------" << "\n";

	if (false)
	{
		cout << "Exiting prematurely from " << __PRETTY_FUNCTION__
			<< "::" << __LINE__ << "!" << endl;
		cerr << "Exiting prematurely from " << __PRETTY_FUNCTION__
			<< "::" << __LINE__ << "!" << endl;
		exit(8);
	}

    }
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;

	swTotal.Stop();
	cout << "Finished function call to " << __FUNCTION__ << "(...) in "
			<< swTotal.printTime() << " s." << endl;

	if (false)
	{
		cout << "Exiting prematurely from " << __PRETTY_FUNCTION__
			<< "::" << __LINE__ << "!" << endl;
		cerr << "Exiting prematurely from " << __PRETTY_FUNCTION__
			<< "::" << __LINE__ << "!" << endl;
		exit(8);
	}


	if (failCounter > 0) exit(-1);

}


//==============================================================================
void SPHWorkstation::initial_smoothing()
{
  // reset linklist to update nearest neighbors
  reset_linklist();

  // reset linklist to update nearest neighbors
  reset_pi_tensor();

  // smooth fields over particles
  smooth_all_particle_fields();

	int count1=0;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;

	for ( auto & p : systemPtr->particles )
	{
    // must reset smoothed charge densities also
    double smoothed_s_lab    = p.hydro.sigma/p.hydro.gamma/settingsPtr->t0;
		double smoothed_rhoB_lab = p.rhoB_sub/p.hydro.gamma/settingsPtr->t0;
		double smoothed_rhoS_lab = p.rhoS_sub/p.hydro.gamma/settingsPtr->t0;
		double smoothed_rhoQ_lab = p.rhoQ_sub/p.hydro.gamma/settingsPtr->t0;

    // UNCOMMENT THIS AND DOCUMENT OUTPUT AS REFERENCE
//    locate_phase_diagram_point_sBSQ( p, smoothed_s_lab, smoothed_rhoB_lab,
//                                      smoothed_rhoS_lab, smoothed_rhoQ_lab );

		fo.frzcheck(p, settingsPtr->t0, count1, systemPtr->n_particles);
	}

	return;

}
///////////////////////////////////////////////////////////////////////////////
// smoothing routines: first smoothing covers all hydrodyanmical fields
void SPHWorkstation::smooth_fields(Particle & pa)
{
  int a = pa.ID;

  pa.hydro.sigma     = 0.0;
  pa.smoothed.s             = 0.0;
  pa.rhoB_sub        = 0.0;
  pa.rhoS_sub        = 0.0;
  pa.rhoQ_sub        = 0.0;

  auto & a_neighbors = systemPtr->linklist.all_neighbors[a];

  for ( int b : a_neighbors )
  {
    auto & pb       = systemPtr->particles[b];

    double kern     = kernel::kernel( pa.r - pb.r, settingsPtr->h );
    pa.hydro.sigma += pb.sigmaweight*kern;
    pa.smoothed.s         += pb.sigmaweight*pb.eta_sigma*kern;
    pa.rhoB_sub    += pb.B*kern;
    pa.rhoS_sub    += pb.S*kern;
    pa.rhoQ_sub    += pb.Q*kern;

    //===============
    // print status
    if ( ( VERBOSE > 2 && pa.print_this_particle )
          || pa.smoothed.s < 0 || isnan( pa.smoothed.s ) )
      std::cout << __FUNCTION__ << "(SPH particle == " << a << "): "
                << systemPtr->t << "   "
                << b << "   " << pa.r
                << "   " << pa.hydro.sigma
                << "   " << pa.smoothed.s
                << "   " << pb.r
                << "   " << pb.sigmaweight
                << "   " << pb.eta_sigma
                << "   " << pb.input.rhoB
                << "   " << pa.rhoB_sub
                << "   " << pb.input.rhoS
                << "   " << pa.rhoS_sub
                << "   " << pb.input.rhoQ
                << "   " << pa.rhoQ_sub
                << "   " << kern << "\n";

  }

  // check if particle has gone nan or negative entropy
  if ( (pa.smoothed.s<0) || isnan(pa.smoothed.s) )
  {
    std::cout << pa.ID <<  " has invalid entropy "
              << pa.T()*hbarc << " " << pa.smoothed.s << endl;
    pa.smoothed.s = TOLERANCE;
  }


  return;
}


///////////////////////////////////////////////////////////////////////////////
//Second smoothing smoothes the gradients after constructing all the fields 
//and derivatives using the equation of state
void SPHWorkstation::smooth_gradients( Particle & pa, double tin )
{
  int a = pa.ID;

  pa.hydro.gradP     = 0.0;
  pa.hydro.gradBulk  = 0.0;
  pa.hydro.gradV     = 0.0;
  pa.hydro.gradshear = 0.0;
  pa.hydro.divshear  = 0.0;

  if ( pa.btrack != -1 ) pa.btrack = 0;

  double rdis = 0;

  auto & a_neighbors = systemPtr->linklist.all_neighbors[a];

  for ( int b : a_neighbors )
  {
    auto & pb                = systemPtr->particles[b];

    Vector<double,2> rel_sep = pa.r - pb.r;
    double rel_sep_norm      = Norm( rel_sep );
    Vector<double,2> gradK   = kernel::gradKernel( rel_sep, rel_sep_norm, settingsPtr->h );
    Vector<double,2> va      = rowp1(0, pa.hydro.shv);
    Vector<double,2> vb      = rowp1(0, pb.hydro.shv);
    Matrix<double,2,2> vminia, vminib;
    mini(vminia, pa.hydro.shv);
    mini(vminib, pb.hydro.shv);

    double sigsqra           = 1.0/(pa.hydro.sigma*pa.hydro.sigma);
    double sigsqrb           = 1.0/(pb.hydro.sigma*pb.hydro.sigma);
    Vector<double,2> sigsigK = pb.sigmaweight * pa.hydro.sigma * gradK;

    pa.hydro.gradP                += ( sigsqrb*pb.p() + sigsqra*pa.p() ) * sigsigK;

    //===============
    // print status
    if ( VERBOSE > 2 && pa.print_this_particle )
      std::cout << "CHECK grads: " << tin << "   "
                << pa.hydro.gradP << "   " << a << "   " << b << "   "
                << sigsqra << "   " << sigsqrb
                << "   " << pa.p() << "   " << pb.p()
                << "   " << pa.get_current_eos_name()
                << "   " << pb.get_current_eos_name()
                << "   " << gradK << "   " << sigsigK
                << "   " << pa.hydro.sigma << "\n";

    double relative_distance_by_h = rel_sep_norm / settingsPtr->h;
    if ( ( relative_distance_by_h <= 2.0 ) && ( a != b ) )
    {
      if ( pa.btrack != -1 ) pa.btrack++;
      if ( pa.btrack ==  1 ) rdis = relative_distance_by_h;
    }

    pa.hydro.gradBulk             += ( pb.hydro.Bulk/pb.hydro.sigma/pb.hydro.gamma
                                  + pa.hydro.Bulk/pa.hydro.sigma/pa.hydro.gamma)/tin*sigsigK;
    pa.hydro.gradV                += (pb.sigmaweight/pa.hydro.sigma)*( pb.hydro.v -  pa.hydro.v )*gradK;

    //===============
    // print status
    if ( VERBOSE > 2 && pa.print_this_particle )
        std::cout << "CHECK gradV: " << tin << "   " << a << "   " << b << "   "
                  << pb.sigmaweight/pa.hydro.sigma << "   " << pb.hydro.v -  pa.hydro.v
                  << "   " << gradK << "   " << pa.hydro.gradV << "\n";

    //===============
    // add shear terms
    if ( settingsPtr->using_shear )
    {
      pa.hydro.gradshear            += inner(sigsigK, pa.hydro.v)*( sigsqrb*vb + sigsqra*va );
      pa.hydro.divshear             += sigsqrb*sigsigK*transpose(vminib)
                                  + sigsqra*sigsigK*transpose(vminia);
    }

    //===============
    // check for nan pressure gradients
    if ( isnan( pa.hydro.gradP(0) ) )
    {
      cout << "gradP stopped working" << endl;
      cout << systemPtr->t <<" "  << pa.hydro.gradP << " " << a << " " << b << endl;
      cout << pb.sigmaweight << " " << pa.hydro.sigma << " " << pb.p() << endl;
      cout << systemPtr->linklist.Size << " " << pb.s() << " " << pa.s() << endl;

      cout << pa.r << endl;
      cout << pb.r << endl;
      cout << kernel::kernel( pa.r - pb.r, settingsPtr->h ) << endl;
    }
    else if ( isnan( pa.hydro.gradP(1) ) )
      cout << "1 " << gradPressure_weight(a, b)
           << " " << a << " " << b << endl;
  }

  const double hc = constants::hbarc_MeVfm;

  if ( ( pa.btrack == 1 )
            && ( ( pa.T()*hc ) >= 150 ) )
    pa.frz2.t=tin;
  else if ( ( pa.btrack == 0 )
            && ( ( pa.T()*hc ) >= 150 )
            && ( pa.Freeze < 4 ) )
    cout << "Missed " << a << " " << tin << "  "
         << pa.T()*hc << " "
         << rdis << " " << systemPtr->cfon <<  endl;

  return;
}




////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::process_initial_conditions()
{
  //============================================================================
  // IMPOSE ENERGY/CHARGE CUTOFFS TO REGULATE EVENT (NO CUTOFFS FOR GUBSER)
  if ( settingsPtr->IC_type != "Gubser"
        && settingsPtr->IC_type != "Gubser_with_shear")
  {

    std::cout << "Length of particles at line " << __LINE__
              << " is " << systemPtr->particles.size() << std::endl;

    //==========================================================================
    // impose the energy cut-off before the initial time step of hydro
    // the original cut should be 0.15
    // try this; NOTE THAT THE 0.00301 IS HARDCODED IN FOR NOW
    systemPtr->particles.erase( std::remove_if(
      systemPtr->particles.begin(),
      systemPtr->particles.end(),
      [](Particle const & p) { return p.input.e <= 0.00301 / hbarc_GeVfm; } ),
      systemPtr->particles.end() );



    std::cout << "Length of particles at line " << __LINE__
              << " is " << systemPtr->particles.size() << std::endl;



    //==========================================================================
    // cut out particles whose energy density is too small for charge densities
    // remove particles with no possible solutions
    systemPtr->particles.erase( std::remove_if(
      systemPtr->particles.begin(),
      systemPtr->particles.end(),
      [this](Particle const & p)  // apply lambda to all particles;
        {                         // check if eBSQ combo has real solution
          return !(this->eos.eBSQ_has_solution_in_conformal_diagonal(
                    p.input.e, p.input.rhoB, p.input.rhoS, p.input.rhoQ ) );
        } ),
      systemPtr->particles.end() );




    std::cout << "Length of particles at line " << __LINE__
              << " is " << systemPtr->particles.size() << std::endl;

  }



  cout << "After e-cutoff and freeze-out: size = " << systemPtr->particles.size() << endl;

//if (1) exit(8);


  // fill out initial particle information
  //int TMP_particle_count = 0;
	for (auto & p : systemPtr->particles)
  {
    // set area element for each SPH particle
    double dA = (settingsPtr->stepx)*(settingsPtr->stepy);

    // Set the rest of particle elements using area element
		//p.u(0)          = 0.0;  // flow must be set in Particle constructor!!!
		//p.u(1)          = 0.0;  // flow must be set in Particle constructor!!!
		p.eta_sigma       = 1.0;
		p.sigmaweight     = dA;
		p.rhoB_weight     = dA;
		p.rhoS_weight     = dA;
		p.rhoQ_weight     = dA;
		p.hydro.Bulk      = 0.0;
		p.B               = p.input.rhoB*dA;
		p.S               = p.input.rhoS*dA;
		p.Q               = p.input.rhoQ*dA;

		// make educated initial guess here for this particle's (T, mu_i) coordinates
		// (improve this in the future)
		p.thermo.T        = 1000.0/hbarc_MeVfm;	// rootfinder seems to work better going downhill than "uphill"
		p.thermo.muB      = 0.0/hbarc_MeVfm;
		p.thermo.muS      = 0.0/hbarc_MeVfm;
		p.thermo.muQ      = 0.0/hbarc_MeVfm;
		p.thermo.eos_name = "default";  // uses whatever the default EoS is

		if ( p.input.e > systemPtr->efcheck )	// impose freeze-out check for e, not s
			p.Freeze=0;
		else
		{
			p.Freeze=4;
			systemPtr->number_part++;
		}
  }

  cout << "After freezeout (redundant): size = "
      << systemPtr->particles.size() - systemPtr->number_part << endl;
  cout << systemPtr->number_part << endl;


  //============================================================================
  // with particles vector now fully initialized, specify or initialize any
  // remaining quantities which depend on this

  // assign particles IDs
  for ( int i = 0; i < systemPtr->particles.size(); i++ )
    systemPtr->particles[i].ID = i;

  // set particles to print
  for ( int & p : settingsPtr->particles_to_print )
    systemPtr->particles[p].print_this_particle = true;
}







































/////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::advance_timestep_rk2( double dt )
{
  systemPtr->rk2 = 1;
  double t0      = systemPtr->t;
  double E0      = systemPtr->Ez;

  // initialize quantities at current time step
  systemPtr->set_current_timestep_quantities();

  ////////////////////////////////////////////
  //    first step
  ////////////////////////////////////////////

  // compute derivatives
  get_time_derivatives();

  // update quantities
  {
    for (int i = 0; i < (int)systemPtr->particles.size(); i++)
    {
      auto & p    = systemPtr->particles[i];

      p.r = systemPtr->r0[i] + 0.5*dt*p.hydro.v;
      if ( p.Freeze < 5 )
      {
        p.hydro.u            = systemPtr->u0[i]        + 0.5*dt*p.hydro.du_dt;
        p.eta_sigma    = systemPtr->etasigma0[i] + 0.5*dt*p.hydro.detasigma_dt;
        p.hydro.Bulk         = systemPtr->Bulk0[i]     + 0.5*dt*p.hydro.dBulk_dt;
        tmini( p.hydro.shv,    systemPtr->shv0[i]      + 0.5*dt*p.hydro.dshv_dt );

        p.contribution_to_total_Ez = systemPtr->particles_E0[i]
                                      + 0.5*dt*p.contribution_to_total_dEz;

        // regulate updated results if necessary
        if ( REGULATE_LOW_T && p.eta_sigma < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.eta_sigma    = systemPtr->etasigma0[i];
      }
    }
  }

  systemPtr->Ez = E0 + 0.5*dt*systemPtr->dEz;
  systemPtr->t  = t0 + 0.5*dt;

  ////////////////////////////////////////////
  //    second step
  ////////////////////////////////////////////

  // compute derivatives
  get_time_derivatives();

  // update quantities
  {
    for (int i = 0; i < (int)systemPtr->particles.size(); i++)
    {
      auto & p    = systemPtr->particles[i];

      p.r = systemPtr->r0[i] + dt*p.hydro.v;
      if ( p.Freeze < 5 )
      {
        p.hydro.u            = systemPtr->u0[i]        + dt*p.hydro.du_dt;
        p.eta_sigma    = systemPtr->etasigma0[i] + dt*p.hydro.detasigma_dt;
        p.hydro.Bulk         = systemPtr->Bulk0[i]     + dt*p.hydro.dBulk_dt;
        tmini( p.hydro.shv,    systemPtr->shv0[i]      + dt*p.hydro.dshv_dt );

        p.contribution_to_total_Ez = systemPtr->particles_E0[i]
                                      + dt*p.contribution_to_total_dEz;

//cout << "CHECK energies: " << i << "   " << t0+dt << "   " << p.r << "   " << p.e() << "   "
//      << systemPtr->particles_E0[i] << "   "
//      << p.contribution_to_total_E << endl;

        // regulate updated results if necessary
        if ( REGULATE_LOW_T && p.eta_sigma < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.eta_sigma    = systemPtr->etasigma0[i];
      }
    }
  }

  systemPtr->Ez = E0 + dt*systemPtr->dEz;
  systemPtr->t  = t0 + dt;

  return;
}




void SPHWorkstation::advance_timestep_rk4( double dt )
{
//  double ets1 = 0.0, ets2 = 0.0, ets3 = 0.0, ets4 = 0.0;
//  double b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;
//  Vector<double,2> k1, k2, k3, k4;
//  Vector<double,2> r1, r2, r3, r4;


  systemPtr->rk2 = 1;
  double t0      = systemPtr->t;
  double E0      = systemPtr->Ez;

  double E1 = 0.0, E2 = 0.0, E3 = 0.0, E4 = 0.0;

  // initialize quantities at current time step
  systemPtr->set_current_timestep_quantities();

  ////////////////////////////////////////////
  //    first step
  ////////////////////////////////////////////

  // compute derivatives
  get_time_derivatives();

  for (int i = 0; i < (int)systemPtr->particles.size(); i++)
  {
    auto & p    = systemPtr->particles[i];

    // store increments
    p.k1        = dt*p.hydro.du_dt;
    p.r1        = dt*p.hydro.v;
    p.ets1      = dt*p.hydro.detasigma_dt;
    p.b1        = dt*p.hydro.dBulk_dt;
    p.hydro.shv1      = dt*p.hydro.dshv_dt;

    // implement increments with appropriate coefficients
    p.hydro.u         = systemPtr->u0[i]        + 0.5*p.k1;
    p.r         = systemPtr->r0[i]        + 0.5*p.r1;
    p.eta_sigma = systemPtr->etasigma0[i] + 0.5*p.ets1;
    p.hydro.Bulk      = systemPtr->Bulk0[i]     + 0.5*p.b1;
    tmini(p.hydro.shv,  systemPtr->shv0[i]      + 0.5*p.hydro.shv1);

    // regulate updated results if necessary
    if ( REGULATE_LOW_T && p.eta_sigma < 0.0
          && p.T() < 50.0/constants::hbarc_MeVfm )
      p.eta_sigma    = systemPtr->etasigma0[i];
  }

  E1           = dt*systemPtr->dEz;


  ////////////////////////////////////////////
  //    second step
  ////////////////////////////////////////////

  systemPtr->t = t0 + 0.5*dt;
  get_time_derivatives();

  for (int i = 0; i < (int)systemPtr->particles.size(); i++)
  {
    auto & p    = systemPtr->particles[i];

    p.k2        = dt*p.hydro.du_dt;
    p.r2        = dt*p.hydro.v;
    p.ets2      = dt*p.hydro.detasigma_dt;
    p.b2        = dt*p.hydro.dBulk_dt;
    p.hydro.shv2      = dt*p.hydro.dshv_dt;

    p.hydro.u         = systemPtr->u0[i]        + 0.5*p.k2;
    p.r         = systemPtr->r0[i]        + 0.5*p.r2;
    p.eta_sigma = systemPtr->etasigma0[i] + 0.5*p.ets2;
    p.hydro.Bulk      = systemPtr->Bulk0[i]     + 0.5*p.b2;
    tmini(p.hydro.shv,  systemPtr->shv0[i]      + 0.5*p.hydro.shv2);

    // regulate updated results if necessary
    if ( REGULATE_LOW_T && p.eta_sigma < 0.0
          && p.T() < 50.0/constants::hbarc_MeVfm )
      p.eta_sigma    = systemPtr->etasigma0[i];
  }

  E2           = dt*systemPtr->dEz;


  ////////////////////////////////////////////
  //    third step
  ////////////////////////////////////////////

  get_time_derivatives();

  for (int i = 0; i < (int)systemPtr->particles.size(); i++)
  {
    auto & p    = systemPtr->particles[i];

    p.k3        = dt*p.hydro.du_dt;
    p.r3        = dt*p.hydro.v;
    p.ets3      = dt*p.hydro.detasigma_dt;
    p.b3        = dt*p.hydro.dBulk_dt;
    p.hydro.shv3      = dt*p.hydro.dshv_dt;


    p.hydro.u         = systemPtr->u0[i]        + p.k3;
    p.r         = systemPtr->r0[i]        + p.r3;
    p.eta_sigma = systemPtr->etasigma0[i] + p.ets3;
    p.hydro.Bulk      = systemPtr->Bulk0[i]     + p.b3;
    tmini(p.hydro.shv,  systemPtr->shv0[i]      + p.hydro.shv3);

    // regulate updated results if necessary
    if ( REGULATE_LOW_T && p.eta_sigma < 0.0
          && p.T() < 50.0/constants::hbarc_MeVfm )
      p.eta_sigma    = systemPtr->etasigma0[i];
  }

  E3           = dt*systemPtr->dEz;

  ////////////////////////////////////////////
  //    fourth step
  ////////////////////////////////////////////

  systemPtr->t = t0 + dt;
  get_time_derivatives();

  constexpr double w1 = 1.0/6.0, w2 = 1.0/3.0;
  for (int i = 0; i < (int)systemPtr->particles.size(); i++)
  {
    auto & p    = systemPtr->particles[i];

    p.k4        = dt*p.hydro.du_dt;
    p.r4        = dt*p.hydro.v;
    p.ets4      = dt*p.hydro.detasigma_dt;
    p.b4        = dt*p.hydro.dBulk_dt;
    p.hydro.shv4      = dt*p.hydro.dshv_dt;

    // sum the weighted steps into yf and return the final y values
    p.hydro.u         = systemPtr->u0[i]        + w1*p.k1   + w2*p.k2   + w2*p.k3   + w1*p.k4;
    p.r         = systemPtr->r0[i]        + w1*p.r1   + w2*p.r2   + w2*p.r3   + w1*p.r4;
    p.eta_sigma = systemPtr->etasigma0[i] + w1*p.ets1 + w2*p.ets2 + w2*p.ets3 + w1*p.ets4;
    p.hydro.Bulk      = systemPtr->Bulk0[i]     + w1*p.b1   + w2*p.b2   + w2*p.b3   + w1*p.b4;
    tmini(p.hydro.shv,  systemPtr->shv0[i]      + w1*p.hydro.shv1 + w2*p.hydro.shv2 + w2*p.hydro.shv3 + w1*p.hydro.shv4);

    // regulate updated results if necessary
    if ( REGULATE_LOW_T && p.eta_sigma < 0.0
          && p.T() < 50.0/constants::hbarc_MeVfm )
      p.eta_sigma    = systemPtr->etasigma0[i];
  }

  E4            = dt*systemPtr->dEz;
  systemPtr->Ez = E0 + w1*E1 + w2*E2 + w2*E3 + w1*E4;


}




////////////////////////////////////////////////////////////////////////////////
int SPHWorkstation::do_freezeout_checks()
{
  int curfrz = 0;
  if ( systemPtr->cfon == 1 )
  {
    // freeze-out checks for all particles
    for ( auto & p : systemPtr->particles )
      fo.frzcheck( p, systemPtr->t, curfrz, systemPtr->n() );

    // update global quantities accordingly
    systemPtr->number_part += curfrz;
    systemPtr->list.resize(curfrz);
  }
  return curfrz;
}



////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::update_freeze_out_lists()
{
  int m = 0;
  for ( auto & p : systemPtr->particles )
    if ( (p.Freeze==3) && (systemPtr->cfon==1) )
    {
      systemPtr->list[m++] = p.ID;
      p.Freeze         = 4;
    }
}


void SPHWorkstation::finalize_freeze_out(int curfrz)
{
  if (systemPtr->cfon==1)
    fo.bsqsvfreezeout( curfrz );


  // keep track of which particles have left EoS grid completely
  // (reset list at end of each timestep)
  systemPtr->particles_out_of_grid.clear();
  for ( auto & p : systemPtr->particles )
    if ( p.Freeze == 5 )
      systemPtr->particles_out_of_grid.push_back( p.ID );

  std::cout << "Summary at t = " << systemPtr->t << ": "
        << systemPtr->particles_out_of_grid.size()
        << " particles have gone out of the EoS grid." << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::get_time_derivatives()
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

  // freeze-out checks here
  int curfrz = do_freezeout_checks();

  //Computes gradients to obtain dsigma/dt
  smooth_all_particle_gradients();

  //calculate time derivatives needed for equations of motion
  evaluate_all_particle_time_derivatives();

  // update freeze out status/lists
  update_freeze_out_lists();

  // finalize frozen out particles
  finalize_freeze_out( curfrz );

  // check/update conserved quantities
  systemPtr->conservation_energy();

  return;
}










////////////////////////////////////////////////////////////////////////////////
double SPHWorkstation::locate_phase_diagram_point_eBSQ( Particle & p,
                 double e_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
  double current_sVal = p.s();
  if ( p.Freeze == 5 )
    return current_sVal;
  else
  {
    // default: use particle's current location as initial guess (pass in corresponding EoS as well!)
    eos.tbqs( p.T(), p.muB(), p.muQ(), p.muS(), p.get_current_eos_name() );
    bool solution_found = false;
    double sVal = eos.s_out( e_In, rhoB_In, rhoS_In, rhoQ_In, solution_found );

    // save results if either default or conformal EoS returned a result
    // (assumes latter always succeeds)
    if ( solution_found )
      eos.set_thermo( p.thermo );
    else
      p.Freeze = 5; // new label for (totally decoupled) particles which go outside grid

    return sVal;
  }
}

////////////////////////////////////////////////////////////////////////////////
double SPHWorkstation::locate_phase_diagram_point_eBSQ(Particle & p, double e_In)
                 { return locate_phase_diagram_point_eBSQ( p, e_In, 0.0, 0.0, 0.0 ); }






////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::locate_phase_diagram_point_sBSQ( Particle & p,
                 double s_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
  if ( p.Freeze != 5 )
  {
    // default: use particle's current location as initial guess
    eos.tbqs( p.T(), p.muB(), p.muQ(), p.muS(), p.get_current_eos_name() );
    bool update_s_success = eos.update_s( s_In, rhoB_In, rhoS_In, rhoQ_In );

    if ( update_s_success )
      eos.set_thermo( p.thermo );
    else
      p.Freeze = 5; // new label for (totally decoupled) particles which go outside grid
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::locate_phase_diagram_point_sBSQ(Particle & p, double s_In) // previously update_s
               { locate_phase_diagram_point_sBSQ( p, s_In, 0.0, 0.0, 0.0 ); }



////////////////////////////////////////////////////////////////////////////////
//  Computes gamma and velocity
void SPHWorkstation::calcbsq(Particle & p)
{
  p.hydro.gamma     = p.gamcalc();
  p.hydro.v         = (1.0/p.hydro.gamma)*p.hydro.u;
  double s_lab      = p.smoothed.s/p.hydro.gamma/systemPtr->t;
  double rhoB_lab   = p.rhoB_sub/p.hydro.gamma/systemPtr->t;
  double rhoS_lab   = p.rhoS_sub/p.hydro.gamma/systemPtr->t;
  double rhoQ_lab   = p.rhoQ_sub/p.hydro.gamma/systemPtr->t;
	locate_phase_diagram_point_sBSQ( p, s_lab, rhoB_lab, rhoS_lab, rhoQ_lab );
}




///////////////////////////////////////////////////////////////////////////////////
double SPHWorkstation::gradPressure_weight(const int a, const int b)
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

  return pb.sigmaweight * pa.hydro.sigma
        * ( pb.p() / (pb.hydro.sigma*pb.hydro.sigma)
          + pa.p() / (pa.hydro.sigma*pb.hydro.sigma) - innerp );
}






////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::setvisc( Particle & p )
{
  tc.setTherm( p.thermo );
  p.hydro.setas     = tc.eta();
  p.hydro.stauRelax = tc.tau_pi();
  p.hydro.zeta      = tc.zeta();
  p.hydro.tauRelax  = tc.tau_Pi();
}