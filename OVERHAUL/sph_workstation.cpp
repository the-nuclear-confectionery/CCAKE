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

	if (true)
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
      /*if (abs(pa.r.x[0])<0.000001 && abs(pa.r.x[1])<0.000001)
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
        << "   " << kern << std::endl;*/

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
//  cout << "CHECK grads: " << tin << "   " << a << "   " << b << "   " << sigsqra << "   " << sigsqrb
//        << "   " << pa.p() << "   " << pb.p() << "   " << gradK << "   " << sigsigK
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
if (a==15962)
{
cout << "CHECK gradV: " << a << "   " << tin << "   " << pa.sigma << "   " << pa.v
		<< "   " << gradK << "   " << b << "   " << pb.sigmaweight << "   " << pb.v
		<< "   " << pb.v -  pa.v << "   " << pa.Bulk << "   " << pb.Bulk << "   " << sigsigK << endl;
}

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
      [](Particle const & p)
        { return !((eosPtr->conformal_diagonal_EoS)->eBSQ_has_solution(
                    p.e_sub, p.rhoB_an, p.rhoS_an, p.rhoQ_an ) );
        } ),
      systemPtr->particles.end() );



    std::cout << "Length of particles at line " << __LINE__
              << " is " << systemPtr->particles.size() << std::endl;

  }



  cout << "After e-cutoff and freeze-out: size = " << systemPtr->particles.size() << endl;

if (1) exit(8);


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
