#include <algorithm>
#include <memory>

#include "constants.h"
#include "eos_conformal_diagonal.h"
#include "sph_workstation.h"
#include "Stopwatch.h"

using namespace constants;

////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::set_EquationOfStatePtr( EquationOfState * eosPtr_in )
{
  eosPtr = eosPtr_in;
}

//void SPHWorkstation::set_EquationsOfMotionPtr( EquationsOfMotion * eomPtr_in )
//{
//  eomPtr = eomPtr_in;
//}

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
void SPHWorkstation::setshear(bool is_first_timestep)
{
    for ( auto & p : systemPtr->particles )
      p.sets(systemPtr->t*systemPtr->t, is_first_timestep);

//if (true) exit(1);
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
    cout << "systemPtr->_n = " << systemPtr->_n << endl;

    for (int i=0; i<systemPtr->_n; i++)
    {
	cout << "----------------------------------------"
			"----------------------------------------" << endl;


      auto & p = systemPtr->particles[i];


		if (settingsPtr->gtyp!=5)
		{
			sw.Start();
			cout << setprecision(12) << "Doing this particle: "
					<< p.r.x[0] << "   " << p.r.x[1] << "\n";

      // solve for the entropy density
			p.s_an = p.locate_phase_diagram_point_eBSQ(
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
                << p.r.x[0] << "   " << p.r.x[1] << "   "
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
			"----------------------------------------" << endl;

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
"p.sets(...)" from BSQHydro and then call smoothing? linklist can
also be set from BSQHydro. It might make more sense to defne a
general function in workstation the loops over SPH paticles and
smooths them using its own various smoothing methods.*/
void SPHWorkstation::initial_smoothing()  // formerly BSQguess()
{
	cout << "setshear..." << endl;
  setshear(true);
	cout << "reset..." << endl;
  systemPtr->reset_linklist();

	cout << "bsqsvoptimization..." << endl;
	bool initialization_mode = true;
	for (int i=0; i<systemPtr->_n; i++)
	{
    auto & p = systemPtr->particles[i];

		smooth_fields(i, initialization_mode);
	}
	cout << "One more loop!" << endl;

	int count1=0;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	/* this might be the only part that needs to stay... but the
	function name should change if it's only setting s_sub and
	some freeze=out checks.. maybe this can all be combined with 
	initialize_entropy_and_charge_densities(...) into one 
	intialize_particle_quantities function*/
	for (int i=0; i<systemPtr->_n; i++)
	{
    auto & p = systemPtr->particles[i];
		p.s_sub = p.sigma/p.gamma/settingsPtr->t0;

    // must reset smoothed charge densities also
		double smoothed_rhoB_lab = p.rhoB_sub/p.gamma/settingsPtr->t0;
		double smoothed_rhoS_lab = p.rhoS_sub/p.gamma/settingsPtr->t0;
		double smoothed_rhoQ_lab = p.rhoQ_sub/p.gamma/settingsPtr->t0;

//if (i==0)
//	cout << "SPH checkpoint c(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
//			<< p.sigmaweight << "   " << p.s_sub << "   "
//			<< p.T() << "   " << p.e() << "   "
//			<< p.p() << "   " << p.s_an << endl;
//		p.locate_phase_diagram_point_sBSQ(
//      p.s_sub, smoothed_rhoB_lab, smoothed_rhoS_lab, smoothed_rhoQ_lab );
//if (i==0)
//	cout << "SPH checkpoint c(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
//			<< p.sigmaweight << "   " << p.s_sub << "   "
//			<< p.T() << "   " << p.e() << "   "
//			<< p.p() << "   " << p.s_an << endl;

		p.sigsub = 0;
		p.frzcheck(settingsPtr->t0, count1, systemPtr->_n);
if (i==0)
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	}

//if (1) exit(1);

	return;

}
///////////////////////////////////////////////////////////////////////////////
// smoothing routines: first smoothing covers all hydrodyanmical fields
void SPHWorkstation::smooth_fields(int a, bool init_mode /*== false*/)
{
  auto & pa    = systemPtr->particles[a];

  pa.sigma           = 0.0;
  pa.eta             = 0.0;
  pa.rhoB_sub        = 0.0;
  pa.rhoS_sub        = 0.0;
  pa.rhoQ_sub        = 0.0;
  //int neighbor_count = 0;

  Vector<int,2> i;
  for ( i.x[0] = -2; i.x[0] <= 2; i.x[0]++ )
  for ( i.x[1] = -2; i.x[1] <= 2; i.x[1]++ )
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
      //if (abs(pa.r.x[0])<0.000001 && abs(pa.r.x[1])<0.000001)
      if ( isnan( pa.eta ) || pa.eta < 0 )
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
        << "   " << kern << std::endl;

      b = systemPtr->linklist.link[b];
    }
  }



/*
  //cout << "Check neighbor count: " << a << "   " << neighbor_count << endl;

  // reset total B, S, and Q charge of each SPH particle to reflect
  // smoothing from kernel function (ONLY ON FIRST TIME STEP)
  //cout << "-----------------------------------------------------------------" << endl;
  if ( init_mode )
  {
    //cout << "BEFORE: " << a << "   " << pa.B << "   "
    //    << pa.S << "   " << pa.Q << endl;
    //cout << pa.rho_weight << "   " << pa.rhoB_an << "   "
    //    << pa.rhoS_an << "   " << pa.rhoQ_an << endl;
    pa.B = pa.rhoB_sub * pa.rhoB_weight;
    pa.S = pa.rhoS_sub * pa.rhoS_weight;
    pa.Q = pa.rhoQ_sub * pa.rhoQ_weight;
    //cout << "AFTER: " << a << "   " << pa.B << "   "
    //    << pa.S << "   " << pa.Q << endl;
    //cout << pa.rho_weight << "   " << pa.rhoB_sub << "   "
    //    << pa.rhoS_sub << "   " << pa.rhoQ_sub << endl;
    //cout << "-----------------------------------------------------------------" << endl;
  }

*/

  return;
}




///////////////////////////////////////////////////////////////////////////////
//Second smoothing smoothes the gradients after constructing all the fields 
//and derivatives using the equation of state
void SPHWorkstation::smooth_gradients( int a, double tin, int & count )
{
  auto & pa    = systemPtr->particles[a];

  pa.gradP     = 0.0;
  pa.gradBulk  = 0.0;
//  pa.gradrhoB  = 0.0;
//  pa.gradrhoS  = 0.0;
//  pa.gradrhoQ  = 0.0;
  pa.gradV     = 0.0;
  pa.gradshear = 0.0;
  pa.divshear  = 0.0;

  Vector<int,2> i;

  if ( pa.btrack != -1 ) pa.btrack = 0;

  double rdis = 0;

  for ( i.x[0] = -2; i.x[0] <= 2; i.x[0]++ )
  for ( i.x[1] = -2; i.x[1] <= 2; i.x[1]++ )
  {

    int b=systemPtr->linklist.lead[
            systemPtr->linklist.triToSum(
              systemPtr->linklist.dael[a] + i, systemPtr->linklist.size ) ];

    while( b != -1 )
    {
      auto & pb          = systemPtr->particles[b];

      Vector<double,2> gradK   = kernel::gradKernel( pa.r - pb.r, settingsPtr->_h );
      Vector<double,2> va      = rowp1(0, pa.shv);
      Vector<double,2> vb      = rowp1(0, pb.shv);
      Matrix<double,2,2> vminia, vminib;

      mini(vminia, pa.shv);
      mini(vminib, pb.shv);

      double sigsqra           = 1.0/(pa.sigma*pa.sigma);
      double sigsqrb           = 1.0/(pb.sigma*pb.sigma);
      Vector<double,2> sigsigK = pb.sigmaweight * pa.sigma * gradK;

      pa.gradP                += ( sigsqrb*pb.p() + sigsqra*pa.p() ) * sigsigK;

//      if (abs(pa.r.x[0])<0.000001 && abs(pa.r.x[1])<0.000001)
//if (a==310||a==3000)
//  cout << "CHECK grads: " << tin << "   " << a << "   " << b << "   " << sigsqra << "   " << sigsqrb
//        << "   " << pa.p() << "   " << pb.p()
//        << "   " << pa.get_current_eos_name()
//        << "   " << pb.get_current_eos_name()
//        << "   " << gradK << "   " << sigsigK
//        << "   " << pa.sigma << endl;

      if ( ( ( Norm( pa.r - pb.r ) / settingsPtr->_h ) <= 2 ) && ( a != b ) )
      {
        if ( pa.btrack != -1 ) pa.btrack++;
        if ( pa.btrack ==  1 ) rdis = Norm(pa.r-pb.r)/settingsPtr->_h;
      }

      pa.gradBulk             += ( pb.Bulk/pb.sigma/pb.gamma
                                    + pa.Bulk/pa.sigma/pa.gamma)/tin*sigsigK;
      //pa.gradrhoB             += ( pb.rhoB/pb.sigma/pb.gamma
      //                            + pa.rhoB/pa.sigma/pa.gamma)/tin*sigsigK;
      //pa.gradrhoS             += ( pb.rhoS/pb.sigma/pb.gamma
      //                            + pa.rhoS/pa.sigma/pa.gamma)/tin*sigsigK;
      //pa.gradrhoQ             += ( pb.rhoQ/pb.sigma/pb.gamma
      //                            + pa.rhoQ/pa.sigma/pa.gamma)/tin*sigsigK;
      pa.gradV                += (pb.sigmaweight/pa.sigma)*( pb.v -  pa.v )*gradK;
//if (a==310||a==3000)
//{
//cout << "CHECK gradV: " << a << "   " << tin << "   " << pa.sigma << "   " << pa.v
//		<< "   " << gradK << "   " << pa.get_current_eos_name()
//    << "   " << pb.get_current_eos_name()
//    << "   " << b << "   " << pb.sigmaweight << "   " << pb.v
//		<< "   " << pb.v -  pa.v << "   " << pa.Bulk << "   " << pb.Bulk << "   " << sigsigK << endl;
//}

      if ( settingsPtr->using_shear )
      {
        pa.gradshear            += inner(sigsigK, pa.v)*( sigsqrb*vb + sigsqra*va );
        pa.divshear             += sigsqrb*sigsigK*transpose(vminib)
                                    + sigsqra*sigsigK*transpose(vminia);
      }

      if ( isnan( pa.gradP.x[0] ) )
      {
        cout << "gradP stopped working" << endl;
        cout << systemPtr->t <<" "  << pa.gradP << " " << a << " " << b << endl;
        cout << pb.sigmaweight << " " << pa.sigma << " " << pb.p() << endl;
        cout << systemPtr->linklist.Size << " " << pb.s() << " " << pa.s() << endl;

        cout << pa.r << endl;
        cout << pb.r << endl;
        cout << kernel::kernel( pa.r - pb.r, settingsPtr->_h ) << endl;
      }
      else if ( isnan( pa.gradP.x[1] ) )
        cout << "1 " << systemPtr->linklist.gradPressure_weight(systemPtr->particles, a, b)
             << " " << a << " " << b << endl;
      else if ( isnan( pa.gradP.x[2] ) )
        cout << "2 " << systemPtr->linklist.gradPressure_weight(systemPtr->particles, a, b)
             << " " << a << " " << b << endl;

      b=systemPtr->linklist.link[b];
    }
  }

  if ( ( pa.btrack == 1 )
        && ( ( pa.T()*197.3 ) >= 150 ) )
    pa.frz2.t=tin;
  else if ( ( pa.btrack == 0 )
            && ( ( pa.T()*197.3 ) >= 150 )
            && ( pa.Freeze < 4 ) )
    cout << "Missed " << a << " " << tin << "  "
         << pa.T()*197.3 << " "
         << rdis << " " << systemPtr->cfon <<  endl;

  return;
}
///////////////////////////////////////////////////////////////////////////////
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
        { return !((this->eosPtr->conformal_diagonal_EoS).eBSQ_has_solution(
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
    // set area element for each SPH particle (for Polar, depends on particle!)
    double dA = 0.0;
    if (settingsPtr->initial_coordinate_distribution == "Cartesian")
    {
      //cout << "stepX = " << settingsPtr->stepx << endl;
      //cout << "stepY = " << settingsPtr->stepy << endl;
      dA = (settingsPtr->stepx)*(settingsPtr->stepy);
    }
    else if (settingsPtr->initial_coordinate_distribution == "Polar")
    {
      //cout << "stepr = " << settingsPtr->stepr << endl;
      //cout << "stepphi = " << settingsPtr->stepphi << endl;
      
      dA = (settingsPtr->stepr)*(settingsPtr->stepphi)*Norm(p.r);
    }
    else
    {
      std::cerr << "ERROR: initial_coordinate_distribution = "
                << settingsPtr->initial_coordinate_distribution
                << " not supported!" << endl;
      exit(1);
    }

//cout << "CHECK PARTICLES: " << p.r.x[0] << "   " << p.r.x[0] << "   "
//      << eLocal << "   " << rhoBLocal << "   " << rhoSLocal << "   "
//      << rhoQLocal << "   " << ux << "   " << uy << endl;

    // Set the rest of particle elements using area element
		//p.u.x[0]          = 0.0;  // flow must be set in Particle constructor!!!
		//p.u.x[1]          = 0.0;  // flow must be set in Particle constructor!!!
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

		/*if (TMP_particle_count++==0)
      cout << "readICs_iccing(" << __LINE__ << "): "
        << "SPH particles: "
        << p.r.x[0] << "   " << p.r.x[1] << "   "
        << p.e_sub << "   " << p.rhoB_an << "   "
        << p.rhoS_an << "   " << p.rhoQ_an << "   "
        << p.sigmaweight << endl;*/

		// make educated initial guess here for this particle's (T, mu_i) coordinates
		// (improve this in the future)
		p.thermo.T        = 50.0/hbarc_MeVfm;	// rootfinder seems to work better going downhill than "uphill"
		p.thermo.muB      = 0.0/hbarc_MeVfm;
		p.thermo.muS      = 0.0/hbarc_MeVfm;
		p.thermo.muQ      = 0.0/hbarc_MeVfm;
		p.thermo.eos_name = "default";  // uses whatever the default EoS is

		if ( p.e_sub > systemPtr->efcheck )	// impose freeze-out check for e, not s
    {
      //cout << "Found " << p.e_sub << " greater than " << systemPtr->efcheck << endl;
			p.Freeze=0;
		}
		else
		{
      //cout << "Found " << p.e_sub << " less than " << systemPtr->efcheck << endl;
			p.Freeze=4;
			systemPtr->number_part++;
    //   cout << "number_part = " << systemPtr->number_part << endl;
		}
  }

  cout << "After freezeout (redundant): size = "
      << systemPtr->particles.size() - systemPtr->number_part << endl;
  cout << systemPtr->number_part << endl;

//if (1) exit(1);
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
  eomPtr->BSQshear(*systemPtr, *this); // will this compile?

  // update quantities
  {
    for (auto & p : systemPtr->particles)
    {
      p.r = systemPtr->r0[i] + 0.5*dx*p.v;
      if ( p.Freeze < 5 )
      {
        p.u            = systemPtr->u0[i]        + 0.5*dx*p.du_dt;
        p.eta_sigma    = systemPtr->etasigma0[i] + 0.5*dx*p.detasigma_dt;
        p.Bulk         = systemPtr->Bulk0[i]     + 0.5*dx*p.dBulk_dt;
        tmini( p.shv,    systemPtr->shv0[i]      + 0.5*dx*p.dshv_dt );
      }
    }
  }

  systemPtr->Ez = E0 + 0.5*dx*systemPtr->dEz;
  systemPtr->t  = t0 + 0.5*dx;

  ////////////////////////////////////////////
  //    second step
  ////////////////////////////////////////////

  // compute derivatives
  eomPtr->BSQshear(system, *this); // will this compile?

  // update quantities
  {
    for (auto & p : systemPtr->particles)
    {
      p.r = systemPtr->r0[i] + dx*p.v;
      if ( p.Freeze < 5 )
      {
        p.u            = systemPtr->u0[i]        + dx*p.du_dt;
        p.eta_sigma    = systemPtr->etasigma0[i] + dx*p.detasigma_dt;
        p.Bulk         = systemPtr->Bulk0[i]     + dx*p.dBulk_dt;
        tmini( p.shv,    systemPtr->shv0[i]      + dx*p.dshv_dt );
      }
    }
  }

  systemPtr->Ez = E0 + dx*systemPtr->dEz;
  systemPtr->t  = t0 + dx;

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
  eomPtr->BSQshear(*systemPtr, *this); // will this compile?

  for (int i = 0; i < (int)particles.size(); i++)
  {
    auto & p    = systemPtr->particles[i];

    // store increments
    p.k1        = dx*p.du_dt;
    p.r1        = dx*p.v;
    p.ets1      = dx*p.detasigma_dt;
    p.b1        = dx*p.dBulk_dt;
    p.shv1      = dx*p.dshv_dt;

    // implement increments with appropriate coefficients
    p.u         = systemPtr->u0[i]        + 0.5*p.k1;
    p.r         = systemPtr->r0[i]        + 0.5*p.r1;
    p.eta_sigma = systemPtr->etasigma0[i] + 0.5*p.ets1;
    p.Bulk      = systemPtr->Bulk0[i]     + 0.5*p.b1;
    tmini(p.shv,  systemPtr->shv0[i]      + 0.5*p.shv1);
  }

  E1           = dx*systemPtr->dEz;


  ////////////////////////////////////////////
  //    second step
  ////////////////////////////////////////////

  systemPtr->t = t0 + 0.5*dx;
  eomPtr->BSQshear(systemPtr, *this);

  for (int i = 0; i < (int)particles.size(); i++)
  {
    auto & p    = systemPtr->particles[i];

    p.k2        = dx*p.du_dt;
    p.r2        = dx*p.v;
    p.ets2      = dx*p.detasigma_dt;
    p.b2        = dx*p.dBulk_dt;
    p.shv2      = dx*p.dshv_dt;

    p.u         = systemPtr->u0[i]        + 0.5*p.k2;
    p.r         = systemPtr->r0[i]        + 0.5*p.r2;
    p.eta_sigma = systemPtr->etasigma0[i] + 0.5*p.ets2;
    p.Bulk      = systemPtr->Bulk0[i]     + 0.5*p.b2;
    tmini(p.shv,  systemPtr->shv0[i]      + 0.5*p.shv2);
  }

  E2           = dx*systemPtr->dEz;


  ////////////////////////////////////////////
  //    third step
  ////////////////////////////////////////////

  eomPtr->BSQshear(systemPtr, *this);

  for (int i = 0; i < (int)particles.size(); i++)
  {
    auto & p    = systemPtr->particles[i];

    p.k3        = dx*p.du_dt;
    p.r3        = dx*p.v;
    p.ets3      = dx*p.detasigma_dt;
    p.b3        = dx*p.dBulk_dt;
    p.shv3      = dx*p.dshv_dt;


    p.u         = systemPtr->u0[i]        + p.k3;
    p.r         = systemPtr->r0[i]        + p.r3;
    p.eta_sigma = systemPtr->etasigma0[i] + p.ets3;
    p.Bulk      = systemPtr->Bulk0[i]     + p.b3;
    tmini(p.shv,  systemPtr->shv0[i]      + p.shv3);
  }

  E3           = dx*systemPtr->dEz;

  ////////////////////////////////////////////
  //    fourth step
  ////////////////////////////////////////////

  systemPtr->t = t0 + dx;
  eomPtr->BSQshear(systemPtr, *this);

  constexpr double w1 = 1.0/6.0, w2 = 1.0/3.0;
  for (int i = 0; i < (int)particles.size(); i++)
  {
    auto & p    = systemPtr->particles[i];

    p.k4        = dx*p.du_dt;
    p.r4        = dx*p.v;
    p.ets4      = dx*p.detasigma_dt;
    p.b4        = dx*p.dBulk_dt;
    p.shv4      = dx*p.dshv_dt;

    // sum the weighted steps into yf and return the final y values
    p.u         = systemPtr->u0[i]        + w1*p.k1   + w2*p.k2   + w2*p.k3   + w1*p.k4;
    p.r         = systemPtr->r0[i]        + w1*p.r1   + w2*p.r2   + w2*p.r3   + w1*p.r4;
    p.eta_sigma = systemPtr->etasigma0[i] + w1*p.ets1 + w2*p.ets2 + w2*p.ets3 + w1*p.ets4;
    p.Bulk      = systemPtr->Bulk0[i]     + w1*p.b1   + w2*p.b2   + w2*p.b3   + w1*p.b4;
    tmini(p.shv,  systemPtr->shv0[i]      + w1*p.shv1 + w2*p.shv2 + w2*p.shv3 + w1*p.shv4);
  }

  E4            = dx*systemPtr->dEz;
  systemPtr->Ez = E0 + w1*E1 + w2*E2 + w2*E3 + w1*E4;


}






////////////////////////////////////////////////////////////////////////////////
// The structure here is temporary until we set the mode for different terms 
// which will be shear, bulk, diffusion, and coupling terms, 
// current equations are only set up for 2+1d.
void SPHWorkstation::BSQshear()
{
  // which particles print extra info for
  constexpr int ic = -1;
  bool printAll = false;


  // not initial call to setshear(bool is_first_timestep)
  setshear(false);
  systemPtr->reset_linklist();
  /* the above two functions should be called in BSQHydro
  (although it's not clear setshear needs to exist) */

  for (int i = 0; i < systemPtr->n(); i++)
  {
    auto & p = systemPtr->particles[i];

//if ( abs(p.r.x[0]) < 0.000001 && abs(p.r.x[1]) < 0.000001 )
//  cout << "CHECK CENTER: " << systemPtr->t << "   " << i << "   " << p.T()*hbarc << "   "
//        << p.eta/p.gamma/systemPtr->t << "   " << p.s() << endl;

    smooth_fields(i);

//if ( abs(p.r.x[0]) < 0.000001 && abs(p.r.x[1]) < 0.000001 )
//  cout << "CHECK CENTER: " << systemPtr->t << "   " << i << "   " << p.T()*hbarc << "   "
//        << p.eta/p.gamma/systemPtr->t << "   " << p.s() << endl;

    if ( (p.eta<0) || isnan(p.eta) )
    {
      cout << i <<  " has invalid entropy " <<  p.T()*hbarc << " " << p.eta << endl;
      p.eta = 0;
      printAll = true;  // turn on verbosity
    }
  }
  /* the above should be called in BSQHydro via something like
  smooth_all_particle_fields() */

  cout << "Finished first loop over SPH particles" << endl;

  int curfrz = 0;
  for ( int i = 0; i < systemPtr->n(); i++ )
  {
    auto & p = systemPtr->particles[i];

    //  Computes gamma and velocity
    p.calcbsq( systemPtr->t ); //resets EOS!!
    /* would be nice to remove the above from eom,
    need to think about where to put it */

    /*N.B. - eventually extend to read in viscosities from table, etc.*/
    p.setvisc( systemPtr->etaconst, systemPtr->bvf, systemPtr->svf,
               systemPtr->zTc,      systemPtr->sTc, systemPtr->zwidth,
               systemPtr->visc );
    /* the above is obsolete when including
     transport_coefficients class */

    if ( systemPtr->cfon == 1 )
      p.frzcheck( systemPtr->t, curfrz, systemPtr->n() );
  /* not sure what cfon is but I'm sure it doesn't need to be here */


  }


  cout << "Finished second loop over SPH particles" << endl;

  if ( systemPtr->cfon == 1 )
  {
    systemPtr->number_part += curfrz;
    systemPtr->list.resize(curfrz);
  }
  /* not sure what cfon is but I'm sure it doesn't need to be here */

  int m=0;
  for ( int i=0; i<systemPtr->n(); i++ )
  {
    auto & p = systemPtr->particles[i];

    //Computes gradients to obtain dsigma/dt
    smooth_gradients( i, systemPtr->t, curfrz );

    p.dsigma_dt = -p.sigma * ( p.gradV.x[0][0] + p.gradV.x[1][1] );
if (i==ic || printAll)
cout << "CHECK dsigma_dt: " << i << "   " << systemPtr->t << "   " << p.dsigma_dt << "   " << p.sigma
		<< "   " << p.gradV << endl;

    p.bsqsvsigset( systemPtr->t, i );

    if ( (p.Freeze==3) && (systemPtr->cfon==1) )
    {
      systemPtr->list[m++] = i;
      p.Freeze         = 4;
    }

  }
  /* the above should probably be put into something like
  smooth_all_particle_gradients() */

  if (systemPtr->rk2==1)
  systemPtr->bsqsvconservation();

  systemPtr->bsqsvconservation_Ez();

  /* conservation can definitely be called in BSQHydro */

//TRAVIS: ALL OF THE ABOVE SHOULD BE SPLIT OFF INTO DIFFERENT FUNCTIONS
// AND CALLED IN BSQHYDRO E.G. smooth_gradients, systemPtr->freeze_out_check, etc

  //calculate matrix elements
  for ( int i=0; i<systemPtr->n(); i++ )
  {
    auto & p = systemPtr->particles[i];

    double gamt = 0.0, pre = 0.0, p1 = 0.0;
    if ( settingsPtr->using_shear )
    {
      gamt = 1.0/p.gamma/p.stauRelax;
      pre  = p.eta_o_tau/p.gamma;
      p1   = gamt - 4.0/3.0/p.sigma*p.dsigma_dt + 1.0/systemPtr->t/3.0;
    }

    Vector<double,2> minshv   = rowp1(0, p.shv);
    Matrix <double,2,2> partU = p.gradU + transpose( p.gradU );

if (i==ic || printAll)
cout << "CHECK misc1: " << i << "   " << systemPtr->t << "   " << gamt << "   " << p.sigma
		<< "   " << p.dsigma_dt << endl;

if (i==ic || printAll)
cout << "CHECK minshv: " << i << "   " << systemPtr->t << "   " << minshv << endl;

if (i==ic || printAll)
cout << "CHECK partU: " << i << "   " << systemPtr->t << "   " << partU << endl;


    // set the Mass and the Force
    Matrix <double,2,2> M = p.Msub(i);
    Vector<double,2> F    = p.Btot*p.u + p.gradshear
                            - ( p.gradP + p.gradBulk + p.divshear );
/* Might make more sense for M and F to be members of particle? Then the
further above loop could be done in workstation and M and F could be set
at the same time... */

if (i==ic || printAll)
cout << "CHECK M: " << i << "   " << systemPtr->t << "   " << M << endl;



if (i==ic || printAll)
cout << "CHECK F: " << i << "   " << systemPtr->t << "   " << F << "   "
		<< p.Btot << "   " << p.u << "   "
		<< p.gradshear << "   " << p.gradP << "   "
		<< p.gradBulk << "   " << p.divshear << endl;

    // shear contribution
    if ( settingsPtr->using_shear )
      F += pre*p.v*partU + p1*minshv;

if (i==ic || printAll)
cout << "CHECK F(again): " << i << "   " << systemPtr->t << "   " << F << "   "
		<< pre << "   " << p.v << "   " << partU << "   "
		<< p1 << "   " << minshv << endl;

    double det=deter(M);


if (i==ic || printAll)
cout << "CHECK det: " << i << "   " << systemPtr->t << "   " << M << "   " << det << endl;


    Matrix <double,2,2> MI;
    MI.x[0][0]=M.x[1][1]/det;
    MI.x[0][1]=-M.x[0][1]/det;
    MI.x[1][0]=-M.x[1][0]/det;
    MI.x[1][1]=M.x[0][0]/det;
  /* This notation is still a bit weird.. but also
  MI should be a member of particle as well */

if (i==ic || printAll)
cout << "CHECK MI: " << i << "   " << systemPtr->t << "   " << MI << endl;


    p.du_dt.x[0]=F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1];
    p.du_dt.x[1]=F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1];

    Matrix <double,2,2> ulpi  = p.u*colp1(0, p.shv);

    double vduk               = inner( p.v, p.du_dt );

    Matrix <double,2,2> Ipi   = -2.0*p.eta_o_tau/3.0 * ( p.Imat + p.uu ) + 4./3.*p.pimin;

    p.div_u                   = (1./ p.gamma)*inner( p.u, p.du_dt)
                                  - ( p.gamma/ p.sigma ) * p.dsigma_dt;

    p.bigtheta                = p.div_u*systemPtr->t+p.gamma;
        /* the above lines could automaticlaly be set in particle after 
        calculating the matrix elements above */

if (i==ic || printAll)
cout << "CHECK div_u: " << i
		<< "   " << systemPtr->t
		<< "   " << p.div_u
		<< "   " << p.gamma
		<< "   " << p.u
		<< "   " << p.du_dt
		<< "   " << inner( p.u, p.du_dt)
		<< "   " << p.sigma 
		<< "   " << p.dsigma_dt << endl;
if (i==ic || printAll)
cout << "CHECK bigtheta: " << i
		<< "   " << systemPtr->t
		<< "   " << p.bigtheta
		<< "   " << p.gamma << endl;

    // this term occurs in Eqs. (250) and (251) of Jaki's long notes
    // translation: pi^{ij} + pi^{00} v^i v^j - pi^{i0} v^j - pi^{0j} v^i
    Matrix <double,2,2> sub   = p.pimin + (p.shv.x[0][0]/p.g2)*p.uu -1./p.gamma*p.piutot;

    // minshv = pi^{0i}                   (i   = 1,2)
    // pimin  = pi^{ij}                   (i,j = 1,2)
    // uu     = u^i u^j                   (i,j = 1,2)
    // piu    = pi^{0i} u^j               (i,j = 1,2)
    // piutot = pi^{0i} u^j + pi^{0j} u^i (i,j = 1,2)
    // gradU  = du_i/dx^j                 (i,j = 1,2)

    if ( settingsPtr->using_shear )
      p.inside                  = systemPtr->t*(
                                inner( -minshv+p.shv.x[0][0]*p.v, p.du_dt )
                                - con2(sub, p.gradU)
                                - p.gamma*systemPtr->t*p.shv33 );


if (i==ic || printAll)
std::cout << "CHECK inside: " << i << "   "
			<< systemPtr->t << "   "
			<< p.inside << "   "
			<< minshv << ";   "
			<< p.shv.x[0][0]*p.v << ";   "
			<< p.du_dt << ";   "
			<< sub << "   "
			<< p.gradU << ";   "
			<< p.gamma*systemPtr->t*p.shv33 << std::endl;



    p.detasigma_dt            = 1./p.sigma/p.T()*( -p.bigPI*p.bigtheta + p.inside );


if (i==ic || printAll)
std::cout << "CHECK detasigma_dt: " << i << "   "
			<< systemPtr->t << "   "
			<< p.detasigma_dt << "   "
			<< p.sigma << "   "
			<< p.T()*hbarc_MeVfm << "   "
			<< p.bigPI << "   "
			<< p.bigtheta << "   "
			<< p.inside << std::endl;


    // N.B. - ADD EXTRA TERMS FOR BULK EQUATION
    p.dBulk_dt = ( -p.zeta/p.sigma*p.bigtheta - p.Bulk/p.gamma )/p.tauRelax;

    Matrix <double,2,2> ududt = p.u*p.du_dt;

    // N.B. - ADD READABLE TERM NAMES
    if ( settingsPtr->using_shear )
      p.dshv_dt                 = - gamt*( p.pimin + p.setas*partU )
                               - p.eta_o_tau*( ududt + transpose(ududt) )
                               + p.dpidtsub() + p.sigl*Ipi
                               - vduk*( ulpi + transpose(ulpi) + (1/p.gamma)*Ipi );

  }


  if (systemPtr->cfon==1)
    systemPtr->bsqsvfreezeout( curfrz );


  // keep track of which particles have left EoS grid completely
  // (reset list at end of each timestep)
  systemPtr->particles_out_of_grid.clear();
  for ( int i = 0; i < systemPtr->n(); i++ )
    if ( systemPtr->particles[i].Freeze == 5 )
      systemPtr->particles_out_of_grid.push_back( i );

  std::cout << "Summary at t = " << systemPtr->t << ": "
        << systemPtr->particles_out_of_grid.size()
        << " particles have gone out of the EoS grid." << std::endl;


  /* Not sure what any of the above does but I'm certain it can be
  done somehwere else */


  return;
}
