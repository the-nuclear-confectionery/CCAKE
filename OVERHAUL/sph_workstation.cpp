#include <algorithm>
#include <memory>

#include "constants.h"
#include "eos_conformal_diagonal.h"
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


/* not sure it makes sense for this to be its own
separate function? */
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

    systemPtr->_n = systemPtr->particles.size();
    cout << "systemPtr->_n = " << systemPtr->_n << "\n";

    for (int i=0; i<systemPtr->_n; i++)
    {
	cout << "----------------------------------------"
			"----------------------------------------" << "\n";


      auto & p = systemPtr->particles[i];


		if (settingsPtr->gtyp!=5)
		{
			sw.Start();
			cout << setprecision(12) << "Doing this particle: "
					<< p.r(0) << "   " << p.r(1) << "\n";

      // solve for the entropy density
			p.s_an = locate_phase_diagram_point_eBSQ( p,
                    p.e_sub, p.rhoB_an, p.rhoS_an, p.rhoQ_an );

			sw.Stop();
			string successString = (p.s_an < 0.0) ?
									"unsuccessfully" : "successfully";
      cout << "Print particle info:\n";
			cout << "    SPH particle " << i << ", locate_phase_diagram_point_eBSQ: completed "
					<< successString << " in " << sw.printTime() << "s." << "\n";

		}

		///////////////////////////////////////////////////////////////////////////
		// for now, if we failed to find a real entropy density for this
		// point, just freeze it out, set its entropy to the freeze-out value,
		// and continue without setting anything else
		if (p.s_an < 0.0)
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
                << p.e_sub*hbarc_MeVfm << "   "
                << p.rhoB_an << "   "
                << p.rhoS_an << "   "
                << p.rhoQ_an << "\n";
			cout << "\t\t - output densities (eBSQ):    "
                << p.e()*hbarc_MeVfm << "   "
                << p.rhoB() << "   "
                << p.rhoS() << "   "
                << p.rhoQ() << "\n";
			cout << "\t\t - freeze-out status:   " << p.Freeze << "\n";
		}

    if (settingsPtr->gtyp==5) p.e_sub = p.e();

    p.gamma=p.gamcalc();

    p.sigmaweight *= p.s_an*p.gamma*settingsPtr->t0;	  // sigmaweight is constant after this
    p.rhoB_weight *= p.gamma*settingsPtr->t0; // rhoB_weight is constant after this
    p.rhoS_weight *= p.gamma*settingsPtr->t0; // rhoS_weight is constant after this
    p.rhoQ_weight *= p.gamma*settingsPtr->t0; // rhoQ_weight is constant after this

		p.B *= p.gamma*settingsPtr->t0;	// B does not evolve in ideal case
		p.S *= p.gamma*settingsPtr->t0;	// S does not evolve in ideal case
		p.Q *= p.gamma*settingsPtr->t0;	// Q does not evolve in ideal case

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
////////////////////////////////////////////////////////////////////////////////
/* Does this need to be its own separate function? Why not call
"p.reset_pi_tensor(...)" from BSQHydro and then call smoothing? linklist can
also be set from BSQHydro. It might make more sense to defne a
general function in workstation the loops over SPH paticles and
smooths them using its own various smoothing methods.*/
void SPHWorkstation::initial_smoothing()  // formerly BSQguess()
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
	/* this might be the only part that needs to stay... but the
	function name should change if it's only setting s_sub and
	some freeze=out checks.. maybe this can all be combined with 
	initialize_entropy_and_charge_densities(...) into one 
	intialize_particle_quantities function*/
	for ( auto & p : systemPtr->particles )
	{
		p.s_sub = p.sigma/p.gamma/settingsPtr->t0;

    // must reset smoothed charge densities also
		double smoothed_rhoB_lab = p.rhoB_sub/p.gamma/settingsPtr->t0;
		double smoothed_rhoS_lab = p.rhoS_sub/p.gamma/settingsPtr->t0;
		double smoothed_rhoQ_lab = p.rhoQ_sub/p.gamma/settingsPtr->t0;

		p.sigsub = 0;
		p.frzcheck(settingsPtr->t0, count1, systemPtr->_n);
	}

	return;

}
///////////////////////////////////////////////////////////////////////////////
// smoothing routines: first smoothing covers all hydrodyanmical fields
void SPHWorkstation::smooth_fields(Particle & pa)
{
  int a = pa.ID;

  pa.sigma           = 0.0;
  pa.eta             = 0.0;
  pa.rhoB_sub        = 0.0;
  pa.rhoS_sub        = 0.0;
  pa.rhoQ_sub        = 0.0;
  //int neighbor_count = 0;

  Vector<int,2> i;
  for ( i(0) = -2; i(0) <= 2; i(0)++ )
  for ( i(1) = -2; i(1) <= 2; i(1)++ )
  {
    int b = systemPtr->linklist.lead[
              systemPtr->linklist.triToSum(
                systemPtr->linklist.dael[a] + i,
                  systemPtr->linklist.size ) ];
    while ( b != -1 )
    {
      const auto & pb = systemPtr->particles[b];
      double kern     = kernel::kernel( pa.r - pb.r, settingsPtr->_h );
      pa.sigma       += pb.sigmaweight*kern;
      pa.eta         += pb.sigmaweight*pb.eta_sigma*kern;
      pa.rhoB_sub    += pb.B*kern;
      pa.rhoS_sub    += pb.S*kern;
      pa.rhoQ_sub    += pb.Q*kern;

      //if (kern>0.0) neighbor_count++;

      //===============
      // print status
      if ( ( VERBOSE > 2
              && settingsPtr->particles_to_print.size() > 0
              && settingsPtr->print_particle(a) )
            || pa.eta < 0 || isnan( pa.eta ) )
        std::cout << __FUNCTION__ << "(SPH particle == " << a << "): "
                  << systemPtr->t << "   "
                  << b << "   " << pa.r
                  << "   " << pa.sigma
                  << "   " << pa.eta
                  << "   " << pb.r
                  << "   " << pb.sigmaweight
                  << "   " << pb.eta_sigma
                  << "   " << pb.rhoB_an
                  << "   " << pa.rhoB_sub
                  << "   " << pb.rhoS_an
                  << "   " << pa.rhoS_sub
                  << "   " << pb.rhoQ_an
                  << "   " << pa.rhoQ_sub
                  << "   " << kern << "\n";

      b = systemPtr->linklist.link[b];
    }
  }

  // check if particle has gone nan or negative entropy
  if ( (pa.eta<0) || isnan(pa.eta) )
  {
    std::cout << pa.ID <<  " has invalid entropy "
              << pa.T()*hbarc << " " << pa.eta << endl;
    pa.eta = TOLERANCE;
  }


  return;
}




///////////////////////////////////////////////////////////////////////////////
//Second smoothing smoothes the gradients after constructing all the fields 
//and derivatives using the equation of state
void SPHWorkstation::smooth_gradients( Particle & pa, double tin, int & count )
{
  int a = pa.ID;

  pa.gradP     = 0.0;
  pa.gradBulk  = 0.0;
  pa.gradV     = 0.0;
  pa.gradshear = 0.0;
  pa.divshear  = 0.0;

  Vector<int,2> i;

  if ( pa.btrack != -1 ) pa.btrack = 0;

  double rdis = 0;

  for ( i(0) = -2; i(0) <= 2; i(0)++ )
  for ( i(1) = -2; i(1) <= 2; i(1)++ )
  {

    int b = systemPtr->linklist.lead[
            systemPtr->linklist.triToSum(
              systemPtr->linklist.dael[a] + i,
              systemPtr->linklist.size ) ];

    while( b != -1 )
    {

      auto & pb          = systemPtr->particles[b];

      Vector<double,2> rel_sep = pa.r - pb.r;
      double rel_sep_norm      = Norm( rel_sep );
      Vector<double,2> gradK   = kernel::gradKernel( rel_sep, rel_sep_norm, settingsPtr->_h );
      //Vector<double,2> gradK   = kernel::gradKernel( rel_sep, settingsPtr->_h );
      Vector<double,2> va      = rowp1(0, pa.shv);
      Vector<double,2> vb      = rowp1(0, pb.shv);
      Matrix<double,2,2> vminia, vminib;
      mini(vminia, pa.shv);
      mini(vminib, pb.shv);

      double sigsqra           = 1.0/(pa.sigma*pa.sigma);
      double sigsqrb           = 1.0/(pb.sigma*pb.sigma);
      Vector<double,2> sigsigK = pb.sigmaweight * pa.sigma * gradK;

      pa.gradP                += ( sigsqrb*pb.p() + sigsqra*pa.p() ) * sigsigK;

      //===============
      // print status
      if ( VERBOSE > 2
            && settingsPtr->particles_to_print.size() > 0
            && settingsPtr->print_particle(a) )
        std::cout << "CHECK grads: " << tin << "   "
                  << pa.gradP << "   " << a << "   " << b << "   "
                  << sigsqra << "   " << sigsqrb
                  << "   " << pa.p() << "   " << pb.p()
                  << "   " << pa.get_current_eos_name()
                  << "   " << pb.get_current_eos_name()
                  << "   " << gradK << "   " << sigsigK
                  << "   " << pa.sigma << "\n";

      double relative_distance_by_h = rel_sep_norm / settingsPtr->_h;
      if ( ( relative_distance_by_h <= 2.0 ) && ( a != b ) )
      {
        if ( pa.btrack != -1 ) pa.btrack++;
        if ( pa.btrack ==  1 ) rdis = relative_distance_by_h;
      }

      pa.gradBulk             += ( pb.Bulk/pb.sigma/pb.gamma
                                    + pa.Bulk/pa.sigma/pa.gamma)/tin*sigsigK;
      pa.gradV                += (pb.sigmaweight/pa.sigma)*( pb.v -  pa.v )*gradK;

      //===============
      // print status
      if ( VERBOSE > 2
            && settingsPtr->particles_to_print.size() > 0
            && settingsPtr->print_particle(a) )
          std::cout << "CHECK gradV: " << tin << "   " << a << "   " << b << "   "
                    << pb.sigmaweight/pa.sigma << "   " << pb.v -  pa.v
                    << "   " << gradK << "   " << pa.gradV << "\n";

      //===============
      // add shear terms
      if ( settingsPtr->using_shear )
      {
        pa.gradshear            += inner(sigsigK, pa.v)*( sigsqrb*vb + sigsqra*va );
        pa.divshear             += sigsqrb*sigsigK*transpose(vminib)
                                    + sigsqra*sigsigK*transpose(vminia);
      }

      //===============
      // check for nan pressure gradients
      if ( isnan( pa.gradP(0) ) )
      {
        cout << "gradP stopped working" << endl;
        cout << systemPtr->t <<" "  << pa.gradP << " " << a << " " << b << endl;
        cout << pb.sigmaweight << " " << pa.sigma << " " << pb.p() << endl;
        cout << systemPtr->linklist.Size << " " << pb.s() << " " << pa.s() << endl;

        cout << pa.r << endl;
        cout << pb.r << endl;
        cout << kernel::kernel( pa.r - pb.r, settingsPtr->_h ) << endl;
      }
      else if ( isnan( pa.gradP(1) ) )
        cout << "1 " << systemPtr->linklist.gradPressure_weight(systemPtr->particles, a, b)
             << " " << a << " " << b << endl;

      b = systemPtr->linklist.link[b];
    }
  }

//if (tin > 1.425) exit(8);

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
      [](Particle const & p) { return p.e_sub <= 0.00301 / hbarc_GeVfm; } ),
      systemPtr->particles.end() );



    std::cout << "Length of particles at line " << __LINE__
              << " is " << systemPtr->particles.size() << std::endl;



    //==========================================================================
    // cut out particles whose energy density is too small for charge densities
    // remove particles with no possible solutions
    systemPtr->particles.erase( std::remove_if(
      systemPtr->particles.begin(),
      systemPtr->particles.end(),
      [this](Particle const & p)
        { return !((this->eos.conformal_diagonal_EoS).eBSQ_has_solution(
                    p.e_sub, p.rhoB_an, p.rhoS_an, p.rhoQ_an ) );
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
		p.Bulk            = 0.0;
		p.B               = p.rhoB_an*dA;
		p.S               = p.rhoS_an*dA;
		p.Q               = p.rhoQ_an*dA;
		p.transverse_area = dA;

		// make educated initial guess here for this particle's (T, mu_i) coordinates
		// (improve this in the future)
		p.thermo.T        = 1000.0/hbarc_MeVfm;	// rootfinder seems to work better going downhill than "uphill"
		p.thermo.muB      = 0.0/hbarc_MeVfm;
		p.thermo.muS      = 0.0/hbarc_MeVfm;
		p.thermo.muQ      = 0.0/hbarc_MeVfm;
		p.thermo.eos_name = "default";  // uses whatever the default EoS is

		if ( p.e_sub > systemPtr->efcheck )	// impose freeze-out check for e, not s
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
  settingsPtr->is_printable.resize( systemPtr->particles.size(), false );
  for ( int & p : settingsPtr->particles_to_print )
  {
    settingsPtr->is_printable[ p ] = true;
    systemPtr->particles[p].print_this_particle = true;
  }
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

      p.r = systemPtr->r0[i] + 0.5*dt*p.v;
      if ( p.Freeze < 5 )
      {
        p.u            = systemPtr->u0[i]        + 0.5*dt*p.du_dt;
        p.eta_sigma    = systemPtr->etasigma0[i] + 0.5*dt*p.detasigma_dt;
        p.Bulk         = systemPtr->Bulk0[i]     + 0.5*dt*p.dBulk_dt;
        tmini( p.shv,    systemPtr->shv0[i]      + 0.5*dt*p.dshv_dt );

        p.contribution_to_total_Ez = systemPtr->particles_E0[i]
                                      + 0.5*dt*p.contribution_to_total_dEz;

        // regulate updated results if necessary
        if ( REGULATE_LOW_T && p.eta_sigma < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.eta_sigma    = systemPtr->etasigma0[i];
        //if ( 0.5*dt*p.detasigma_dt > 10.0*systemPtr->etasigma0[i]
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

      p.r = systemPtr->r0[i] + dt*p.v;
      if ( p.Freeze < 5 )
      {
        p.u            = systemPtr->u0[i]        + dt*p.du_dt;
        p.eta_sigma    = systemPtr->etasigma0[i] + dt*p.detasigma_dt;
        p.Bulk         = systemPtr->Bulk0[i]     + dt*p.dBulk_dt;
        tmini( p.shv,    systemPtr->shv0[i]      + dt*p.dshv_dt );

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
    p.k1        = dt*p.du_dt;
    p.r1        = dt*p.v;
    p.ets1      = dt*p.detasigma_dt;
    p.b1        = dt*p.dBulk_dt;
    p.shv1      = dt*p.dshv_dt;

    // implement increments with appropriate coefficients
    p.u         = systemPtr->u0[i]        + 0.5*p.k1;
    p.r         = systemPtr->r0[i]        + 0.5*p.r1;
    p.eta_sigma = systemPtr->etasigma0[i] + 0.5*p.ets1;
    p.Bulk      = systemPtr->Bulk0[i]     + 0.5*p.b1;
    tmini(p.shv,  systemPtr->shv0[i]      + 0.5*p.shv1);

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

    p.k2        = dt*p.du_dt;
    p.r2        = dt*p.v;
    p.ets2      = dt*p.detasigma_dt;
    p.b2        = dt*p.dBulk_dt;
    p.shv2      = dt*p.dshv_dt;

    p.u         = systemPtr->u0[i]        + 0.5*p.k2;
    p.r         = systemPtr->r0[i]        + 0.5*p.r2;
    p.eta_sigma = systemPtr->etasigma0[i] + 0.5*p.ets2;
    p.Bulk      = systemPtr->Bulk0[i]     + 0.5*p.b2;
    tmini(p.shv,  systemPtr->shv0[i]      + 0.5*p.shv2);

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

    p.k3        = dt*p.du_dt;
    p.r3        = dt*p.v;
    p.ets3      = dt*p.detasigma_dt;
    p.b3        = dt*p.dBulk_dt;
    p.shv3      = dt*p.dshv_dt;


    p.u         = systemPtr->u0[i]        + p.k3;
    p.r         = systemPtr->r0[i]        + p.r3;
    p.eta_sigma = systemPtr->etasigma0[i] + p.ets3;
    p.Bulk      = systemPtr->Bulk0[i]     + p.b3;
    tmini(p.shv,  systemPtr->shv0[i]      + p.shv3);

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

    p.k4        = dt*p.du_dt;
    p.r4        = dt*p.v;
    p.ets4      = dt*p.detasigma_dt;
    p.b4        = dt*p.dBulk_dt;
    p.shv4      = dt*p.dshv_dt;

    // sum the weighted steps into yf and return the final y values
    p.u         = systemPtr->u0[i]        + w1*p.k1   + w2*p.k2   + w2*p.k3   + w1*p.k4;
    p.r         = systemPtr->r0[i]        + w1*p.r1   + w2*p.r2   + w2*p.r3   + w1*p.r4;
    p.eta_sigma = systemPtr->etasigma0[i] + w1*p.ets1 + w2*p.ets2 + w2*p.ets3 + w1*p.ets4;
    p.Bulk      = systemPtr->Bulk0[i]     + w1*p.b1   + w2*p.b2   + w2*p.b3   + w1*p.b4;
    tmini(p.shv,  systemPtr->shv0[i]      + w1*p.shv1 + w2*p.shv2 + w2*p.shv3 + w1*p.shv4);

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
      p.frzcheck( systemPtr->t, curfrz, systemPtr->n() );

    // update global quantities accordingly
    systemPtr->number_part += curfrz;
    systemPtr->list.resize(curfrz);
  }
  return curfrz;
}



//==============================================================================
void SPHWorkstation::update_all_particles_dsigma_dt()
{
  for ( auto & p : systemPtr->particles )
  {
    p.dsigma_dt = -p.sigma * ( p.gradV(0,0) + p.gradV(1,1) );

    //===============
    // print status
    if ( VERBOSE > 2
          && settingsPtr->particles_to_print.size() > 0
          && settingsPtr->print_particle(p.ID) )
      std::cout << "CHECK dsigma_dt: " << p.ID << "   " << systemPtr->t << "   "
                << p.dsigma_dt << "   " << p.sigma << "   " << p.gradV << "\n";
  }
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
    bsqsvfreezeout( curfrz );


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
  smooth_all_particle_gradients( curfrz );

  // set dsigma_dt
  update_all_particles_dsigma_dt();

  // update fluid quantities
  update_all_particle_fluid_quantities();

  // update freeze out status/lists
  update_freeze_out_lists();

  // check/update conserved quantities
  systemPtr->conservation_energy();

  //calculate time derivatives needed for equations of motion
  evaluate_all_particle_time_derivatives();

  // finalize frozen out particles
  finalize_freeze_out( curfrz );

  return;
}










////////////////////////////////////////////////////////////////////////////////
double SPHWorkstation::locate_phase_diagram_point_eBSQ( Particle & p,
                 double e_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
  double current_sVal = s();
  if ( Freeze == 5 )
    return current_sVal;
  else
  {
    // default: use particle's current location as initial guess (pass in corresponding EoS as well!)
    eos.tbqs( p.T, p.muB, p.muQ, p.muS, p.eos_name );
    bool solution_found = false;
    double sVal = eos.s_out( e_In, rhoB_In, rhoS_In, rhoQ_In, solution_found );

    // save results if either default or conformal EoS returned a result
    // (assumes latter always succeeds)
    if ( solution_found )
      set_particle_thermo(p);
    else
      Freeze = 5; // new label for (totally decoupled) particles which go outside grid

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
  if ( Freeze != 5 )
  {
    // default: use particle's current location as initial guess
    eos.tbqs( p.T, p.muB, p.muQ, p.muS, p.eos_name );
    bool update_s_success = eos.update_s( s_In, rhoB_In, rhoS_In, rhoQ_In );

    if ( update_s_success )
      set_particle_thermo(p);
    else
      Freeze = 5; // new label for (totally decoupled) particles which go outside grid
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::locate_phase_diagram_point_sBSQ(Particle & p, double s_In) // previously update_s
               { locate_phase_diagram_point_sBSQ( p, s_In, 0.0, 0.0, 0.0 ); }



void SPHWorkstation::set_particle_thermo(Particle & p)
{
  p.eos_name = eos.get_current_eos_name();

  p.T    = eos.T();
  p.muB  = eos.muB();
  p.muS  = eos.muS();
  p.muQ  = eos.muQ();

  p.p    = eos.p();
  p.s    = eos.s();
  p.rhoB = eos.B();
  p.rhoS = eos.S();
  p.rhoQ = eos.Q();
  p.e    = eos.e();
  p.w    = eos.w();
  p.A    = eos.A();
  p.cs2  = eos.cs2();
  p.dwds = eos.dwds();
  p.dwdB = eos.dwdB();
  p.dwdS = eos.dwdS();
  p.dwdQ = eos.dwdQ();
}


////////////////////////////////////////////////////////////////////////////////
//  Computes gamma and velocity
void SPHWorkstation::calcbsq(Particle & p)
{
  p.gamma           = p.gamcalc();
  p.v               = (1.0/p.gamma)*p.u;
  double s_lab    = p.eta/p.gamma/systemPtr->t;
  double rhoB_lab = p.rhoB_sub/p.gamma/systemPtr->t;
  double rhoS_lab = p.rhoS_sub/p.gamma/systemPtr->t;
  double rhoQ_lab = p.rhoQ_sub/p.gamma/systemPtr->t;
  p.qmom            = ( (p.e()+p.p())*p.gamma/p.sigma )*p.u;
	locate_phase_diagram_point_sBSQ( p, s_lab, rhoB_lab, rhoS_lab, rhoQ_lab );
}





//////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::bsqsvfreezeout(int curfrz)
{
  cout << "CHECK BSQSVFREEZEOUT: " << systemPtr->frzc
        << "   " << systemPtr->tau << "   " << systemPtr->taup
        << "   " << systemPtr->taupp << "   " << systemPtr->cfon << endl;

  if (systemPtr->frzc==0)
  {
    systemPtr->taupp = systemPtr->t;
    systemPtr->frzc  = 1;
    for (auto & p : systemPtr->particles)
    {
      p.frz2.r       = p.r;
      p.frz2.u       = p.u;
      p.frz2.sigma   = p.sigma;
      p.frz2.T       = p.T();
      p.frz2.muB     = p.muB();
      p.frz2.muS     = p.muS();
      p.frz2.muQ     = p.muQ();
      p.frz2.e       = p.e();
      p.frz2.rhoB    = p.rhoB();
      p.frz2.rhoS    = p.rhoS();
      p.frz2.rhoQ    = p.rhoQ();
      p.frz2.bulk    = p.bigPI;
      p.frz2.theta   = p.div_u + p.gamma/t;
      p.frz2.gradP   = p.gradP;
      p.frz2.shear   = p.shv;
      p.frz2.shear33 = p.shv33;
      p.frz2.inside  = p.inside;
    }

  }
  else if (systemPtr->frzc==1)
  {
    systemPtr->taup = systemPtr->t;
    systemPtr->frzc = 2;
    for (auto & p : systemPtr->particles)
    {
      p.frz1.r       = p.r;
      p.frz1.u       = p.u;
      p.frz1.sigma   = p.sigma;
      p.frz1.T       = p.T();
      p.frz1.muB     = p.muB();
      p.frz1.muS     = p.muS();
      p.frz1.muQ     = p.muQ();
      p.frz1.e       = p.e();
      p.frz1.rhoB    = p.rhoB();
      p.frz1.rhoS    = p.rhoS();
      p.frz1.rhoQ    = p.rhoQ();
      p.frz1.bulk    = p.bigPI;
      p.frz1.theta   = p.div_u + p.gamma/t;
      p.frz1.gradP   = p.gradP;
      p.frz1.shear   = p.shv;
      p.frz1.shear33 = p.shv33;
      p.frz1.inside  = p.inside;
    }

    systemPtr->divTtemp.resize( curfrz );
    systemPtr->divT.resize( curfrz );
    systemPtr->gsub.resize( curfrz );
    systemPtr->uout.resize( curfrz );
    systemPtr->swsub.resize( curfrz );
    systemPtr->bulksub.resize( curfrz );
    systemPtr->shearsub.resize( curfrz );
    systemPtr->shear33sub.resize( curfrz );
    systemPtr->tlist.resize( curfrz );
    systemPtr->rsub.resize( curfrz );

    if ( curfrz > 0 )
      bsqsvinterpolate( curfrz );
    else
      systemPtr->cf = 0;
  }
  else
  {
    int i_local = 0;
    for (auto & p : systemPtr->particles)
    {
      if ( p.Freeze < 4 )
      {
        if ( ( p.btrack <= 3 ) && ( p.btrack > 0 ) )
        {
          p.fback4 = p.fback2;
          p.fback3 = p.fback;
          p.fback2 = p.frz2;
          p.fback  = p.frz1;
        }
        else if ( p.btrack == 0 )
        {
          if ( p.fback.gradP(0) != 0 )
          {
            p.frz2 = p.fback2;
            p.frz1 = p.fback;
          }
          else
          {
            p.frz2 = p.fback4;
            p.frz1 = p.fback3;
            cout << "back second"  << endl;
          }


          curfrz++;
          systemPtr->list.push_back( i_local );
          p.Freeze = 4;
          p.btrack = -1;
        }
      }

      i_local++;
    }

    systemPtr->tau = systemPtr->t;

    // resize vectors
    systemPtr->divTtemp.resize( curfrz );
    systemPtr->divT.resize( curfrz );
    systemPtr->gsub.resize( curfrz );
    systemPtr->uout.resize( curfrz );
    systemPtr->swsub.resize( curfrz );
    systemPtr->bulksub.resize( curfrz );
    systemPtr->shearsub.resize( curfrz );
    systemPtr->shear33sub.resize( curfrz );
    systemPtr->tlist.resize( curfrz );
    systemPtr->rsub.resize( curfrz );

    if ( curfrz > 0 )
      bsqsvinterpolate( curfrz );
    else
      systemPtr->cf = 0;


    //sets up the variables for the next time step
    for (auto & p : systemPtr->particles)
    {
      p.frz2         = p.frz1;

      p.frz1.r       = p.r;
      p.frz1.u       = p.u;
      p.frz1.sigma   = p.sigma;
      p.frz1.T       = p.T();
      p.frz1.muB     = p.muB();
      p.frz1.muS     = p.muS();
      p.frz1.muQ     = p.muQ();
      p.frz1.e       = p.e();
      p.frz1.rhoB    = p.rhoB();
      p.frz1.rhoS    = p.rhoS();
      p.frz1.rhoQ    = p.rhoQ();
      p.frz1.bulk    = p.bigPI ;
      p.frz1.theta   = p.div_u+p.gamma/systemPtr->t;
      p.frz1.gradP   = p.gradP;
      p.frz1.shear   = p.shv;
      p.frz1.shear33 = p.shv33;
      p.frz1.inside  = p.inside;
    }

    systemPtr->taupp = systemPtr->taup;
    systemPtr->taup  = systemPtr->tau;
  }

  systemPtr->cfon = 0;
}


//////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::bsqsvinterpolate(int curfrz)
{
  systemPtr->sFO.resize( curfrz, 0 );
  systemPtr->Tfluc.resize( curfrz, 0 );

  for (int j=0; j<curfrz; j++)
  {
    int i    = systemPtr->list[j];
    auto & p = systemPtr->particles[i];

    int swit = 0;
    if ( abs( p.frz1.e - systemPtr->efcheck ) < abs( p.frz2.e - systemPtr->efcheck ) )
      swit   = 1;
    else
      swit   = 2;

    double sigsub = 0.0, thetasub = 0.0, inside = 0.0;
    Vector<double,2> gradPsub;
    if ( swit == 1 )
    {
      if ( p.btrack != -1 )
        systemPtr->tlist[j]    = systemPtr->taup;
      else
        systemPtr->tlist[j]    = systemPtr->taup - systemPtr->dt;

      systemPtr->rsub[j]       = p.frz1.r;
      systemPtr->uout[j]       = p.frz1.u;
      systemPtr->bulksub[j]    = p.frz1.bulk;
      systemPtr->shearsub[j]   = p.frz1.shear;
      systemPtr->shear33sub[j] = p.frz1.shear33;

      gradPsub      = p.frz1.gradP;
      inside        = p.frz1.inside;
      sigsub        = p.frz1.sigma;
      thetasub      = p.frz1.theta;
      systemPtr->Tfluc[j]      = p.frz1.T;             // replace with e
    }
    else if ( swit == 2 )
    {
      if ( p.btrack != -1 )
        systemPtr->tlist[j]    = systemPtr->taupp;
      else
        systemPtr->tlist[j]    = systemPtr->taupp - systemPtr->dt;

      systemPtr->rsub[j]       = p.frz2.r;
      systemPtr->uout[j]       = p.frz2.u;
      systemPtr->bulksub[j]    = p.frz2.bulk;
      systemPtr->shearsub[j]   = p.frz2.shear;
      systemPtr->shear33sub[j] = p.frz2.shear33;

      gradPsub      = p.frz2.gradP;
      inside        = p.frz2.inside;
      sigsub        = p.frz2.sigma;
      thetasub      = p.frz2.theta;
      systemPtr->Tfluc[j]      = p.frz2.T;           // replace with e
    }
    else
    {
      cout << __PRETTY_FUNCTION__ << ": Not at freeze-out temperature" << endl;
    }

    systemPtr->sFO[j]       = eos.s_terms_T( Tfluc[j] );  // replace with e, BSQ

    systemPtr->gsub[j]      = sqrt( Norm2(systemPtr->uout[j]) + 1 );


    sigsub      /= systemPtr->gsub[j]*systemPtr->tlist[j];
    swsub[j]     = p.sigmaweight/sigsub;

    systemPtr->divT[j]      = (1.0/systemPtr->sFO[j])*gradPsub;
    systemPtr->divTtemp[j]  = -(1.0/(systemPtr->gsub[j]*systemPtr->sFO[j]))
                      *( cs2 * (systemPtr->wfz+bulksub[j]) * thetasub
                        - cs2*inside+inner(uout[j], gradPsub) );
//THIS NEEDS TO BE RESET


    double insub = systemPtr->divTtemp[j]*systemPtr->divTtemp[j] - Norm2(systemPtr->divT[j]);
    double norm  = -sqrt(abs(insub));
    systemPtr->divTtemp[j] /= norm;
    systemPtr->divT[j]      = (1.0/norm)*systemPtr->divT[j];


    if ( systemPtr->divTtemp[j] == 1 )
    {
      cout << "track sph=" << p.btrack << " " << i << endl;
      cout << systemPtr->divTtemp[j] << " " << systemPtr->divT[j] << " " << norm << endl;
      cout << gradPsub << " " << thetasub << endl;
      cout << systemPtr->tlist[j] << " " << p.r << endl;
      cout << p.frz1.gradP << " " << p.frz2.gradP << endl;
      cout << p.frz1.T*197.3<< " " << p.frz2.T*197.3 << endl;
      getchar();
    }

    avgetasig += systemPtr->sFO[j]/sigsub;

    if(isnan(systemPtr->divTtemp[j]))
    {
      cout << "divtemp" << endl;
      cout << systemPtr->divTtemp[j] << " " << systemPtr->divT[j] << " " << norm << endl;
      cout << gradPsub << " " << thetasub << endl;
      cout << bulksub[j] << endl;
      cout << gsub[j] << endl;
      cout << systemPtr->tlist[j] << " " << p.r << endl;
      cout << p.frz1.T*0.1973<< " " << p.frz2.T*0.1973<< endl;
    }

    systemPtr->sFO[j]   *= pow(systemPtr->Tfluc[j]*0.1973, 3);
    systemPtr->Tfluc[j] *= 0.1973;

  }

  systemPtr->cf = curfrz;
}
///////////////////////////////////////////////////////////////////////////////////
