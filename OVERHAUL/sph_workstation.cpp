#include "constants.h"
#include "sph_workstation.h"
#include "Stopwatch.h"

// functions calls to static EoS C library
#include <lib.h>
#include "eos_delaunay/eos_delaunay.h"

using namespace constants;

void SPHWorkstation::set_SystemStatePtr( SystemState * systemPtr_in )
{
  systemPtr = systemPtr_in;
}


void SPHWorkstation::set_SettingsPtr( Settings * settingsPtr_in )
{
  settingsPtr = settingsPtr_in;
}


void SPHWorkstation::setshear()
{
    for ( auto & p : systemPtr->particles )
      p.sets(systemPtr->t*systemPtr->t);
}


////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::initialize_entropy_and_charge_densities() // formerly updateIC
{
	// set up EoS C library
	initialize_thermodynamics("/projects/jnorhos/BSQ/EoS_BQS_Derivatives/Coefficients_Parameters.dat");
	Stopwatch sw, swTotal;
	swTotal.Start();
	long long failCounter = 0;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;

    cout << "systemPtr->_n = " << systemPtr->_n << endl;

    for (int i=0; i<systemPtr->_n; i++)
    {
      auto & p = systemPtr->particles[i];



		cout << "----------------------------------------"
				"----------------------------------------" << endl;

		if (settingsPtr->gtyp!=5)
		{
			sw.Start();
			cout << "Doing this particle: "
					<< p.r.x[0] << "   " << p.r.x[1] << "\n";

cout << "SPH checkpoint(" << __LINE__ << "): " << endl;
cout << i << endl;
cout << systemPtr->t << endl;
cout << p.sigmaweight << endl;
cout << p.e_sub << endl;
cout << p.eosPtr->T() << endl;
cout << p.eosPtr->e() << endl;
cout << p.eosPtr->p() << endl;
cout << p.s_an << endl;

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
			p.s_an = p.locate_phase_diagram_point_eBSQ(
                    p.e_sub, p.rhoB_an, p.rhoS_an, p.rhoQ_an );
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;

			if (true || settingsPtr->VERBOSE>5)
			{
				if (p.s_an>0.0)
				{
					double phase_diagram_point[4]
						= { p.eosPtr->T()*197.3,
							  p.eosPtr->muB()*197.3,
							  p.eosPtr->muS()*197.3,
							  p.eosPtr->muQ()*197.3 };
					double densities_at_point[4];
					get_eBSQ_densities(phase_diagram_point, densities_at_point);
					cout << i << ":   " << p.e_sub*197.3
						<< "   " << p.rhoB_an
						<< "   " << p.rhoS_an
						<< "   " << p.rhoQ_an
						<< "   " << p.eosPtr->T()*197.3
						<< "   " << p.eosPtr->muB()*197.3
						<< "   " << p.eosPtr->muS()*197.3
						<< "   " << p.eosPtr->muQ()*197.3;
						for (int iii = 0; iii < 4; iii++)
							cout << "   " << densities_at_point[iii];		
					cout << "\n";
				}
				else
					cout << i << ":   " << p.e_sub*197.3
						<< "   " << p.rhoB_an
						<< "   " << p.rhoS_an
						<< "   " << p.rhoQ_an
						<< "   " << 0.0 << "   " << 0.0
						<< "   " << 0.0 << "   " << 0.0
						<< "   nan   nan   nan   nan\n";
			}

			sw.Stop();
			string successString = (p.s_an < 0.0) ?
									"unsuccessfully" : "successfully";
			cout << "SPH particle " << i << ", locate_phase_diagram_point_eBSQ: completed "
					<< successString << " in " << sw.printTime() << "s." << "\n";

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
		}

		///////////////////////////////////////////////////////////////////////////
		// for now, if we failed to find a real entropy density for this
		// point, just freeze it out, set its entropy to the freeze-out value,
		// and continue without setting anything else
		if (p.s_an < 0.0)
		{

			/////////////////////////////////////////////////////////////////////////
			// if failed with charge densities, set them to zero and re-solve;
			// if that fails too, guesstimate an answer
			cout << "\t --> Densities not found in EoS table (setting BSQ --> 0): "
					<< p.r.x[0] << "   " << p.r.x[1] << "\n"
					<< "\t\t - densities: "
					<< p.e_sub*197.3 << "   " << p.rhoB_an << "   "
					<< p.rhoS_an << "   " << p.rhoQ_an << "\n";

			// set charge densities to zero and re-solve
			p.rhoB_an = 0.0;
			p.rhoS_an = 0.0;
			p.rhoQ_an = 0.0;

			p.s_an = p.locate_phase_diagram_point_eBSQ(
                    p.e_sub, p.rhoB_an, p.rhoS_an, p.rhoQ_an );

			// if this fails too...
			if (p.s_an < 0.0)
			{
				double scale_factor = std::min( 1.0, p.e_sub / systemPtr->efcheck );
	
				cout << "\t\t - scaling e to get s: "
						<< systemPtr->efcheck*0.1973 << "   "
						<< systemPtr->sfcheck << "   "
						<< scale_factor << "   "
						<< scale_factor * systemPtr->sfcheck << "\n";
	
				p.s_an = scale_factor * systemPtr->sfcheck;
			}
			else	// if a solution was found
			{
				cout << "\t\t - phase diagram point: "
						<< p.eosPtr->T()*197.3 << "   "
						<< p.eosPtr->muB()*197.3 << "   "
						<< p.eosPtr->muS()*197.3 << "   "
						<< p.eosPtr->muQ()*197.3 << "\n";
			}

			// freeze this particle out!
			p.Freeze = 4;
			systemPtr->number_part++;
			////////////////////////////////////////////////////////
		}
		else
		{
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;

			cout << "\t --> Densities found in EoS table: "
				<< p.r.x[0] << "   " << p.r.x[1] << "\n";
			cout << "\t\t - phase diagram point: "
					<< p.eosPtr->T()*197.3 << "   "
					<< p.eosPtr->muB()*197.3 << "   "
					<< p.eosPtr->muS()*197.3 << "   "
					<< p.eosPtr->muQ()*197.3 << "\n";
			cout << "\t\t - densities: "
					<< p.e_sub*197.3 << "   "
					<< p.rhoB_an << "   "
					<< p.rhoS_an << "   "
					<< p.rhoQ_an << "\n";
			
			cout << "\t --> Exact:\n";
			double phase_diagram_point[4] = { p.eosPtr->T()*197.3,
											  p.eosPtr->muB()*197.3,
											  p.eosPtr->muS()*197.3,
											  p.eosPtr->muQ()*197.3 };
			double densities_at_point[4];
			get_eBSQ_densities(phase_diagram_point, densities_at_point);
			cout << "\t\t - phase diagram point:";
			for (int iii = 0; iii < 4; iii++) cout << "   " << phase_diagram_point[iii];
			cout << "\n\t\t - densities:";
			for (int iii = 0; iii < 4; iii++) cout << "   " << densities_at_point[iii];
			cout << "\n";
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;

		}

    if (settingsPtr->gtyp==5) p.e_sub = p.eosPtr->e();

    p.gamma=p.gamcalc();

    p.sigmaweight *= p.s_an*p.gamma*settingsPtr->t0;	// sigmaweight is constant after this
    //p.rho_weight *= p.gamma*t0;				// rho_weight is constant after this

		p.B *= p.gamma*settingsPtr->t0;	// B does not evolve in ideal case (confirm with Jaki)
		p.S *= p.gamma*settingsPtr->t0;	// S does not evolve in ideal case (confirm with Jaki)
		p.Q *= p.gamma*settingsPtr->t0;	// Q does not evolve in ideal case (confirm with Jaki)

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;

    }
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;

	swTotal.Stop();
	cout << "Finished function call to updateIC(...) in "
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
void SPHWorkstation::initial_smoothing()  // formerly BSQguess()
{
	cout << "setshear..." << endl;
  setshear();
	cout << "initiate..." << endl;
  //initiate();

	cout << "bsqsvoptimization..." << endl;
	bool initialization_mode = true;
	for (int i=0; i<systemPtr->_n; i++)
	{
    auto & p = systemPtr->particles[i];

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.s_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
		smooth_fields(i, initialization_mode);
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.s_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
	}
	cout << "One more loop!" << endl;

	int count1=0;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	for (int i=0; i<systemPtr->_n; i++)
	{
    auto & p = systemPtr->particles[i];
		p.s_sub = p.sigma/p.gamma/settingsPtr->t0;

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.s_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
		p.locate_phase_diagram_point_sBSQ(
      p.s_sub, p.rhoB_sub, p.rhoS_sub, p.rhoQ_sub );
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << systemPtr->t << "   "
			<< p.sigmaweight << "   " << p.s_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;

		p.sigsub = 0;
		p.frzcheck(settingsPtr->t0, count1, systemPtr->_n);
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	}
	return;




  cout << "BSQ-SV simulation completed!" << endl;

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
  int neighbor_count = 0;

  Vector<int,2> i;
  for ( i.x[0] =- 2; i.x[0] <= 2; i.x[0]++ )
  for ( i.x[1] =- 2; i.x[1] <= 2; i.x[1]++ )
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
      pa.rhoB_sub    += pb.rho_weight*pb.rhoB_an*kern;    //confirm with Jaki
      pa.rhoS_sub    += pb.rho_weight*pb.rhoS_an*kern;    //confirm with Jaki
      pa.rhoQ_sub    += pb.rho_weight*pb.rhoQ_an*kern;    //confirm with Jaki

      if (kern>0.0) neighbor_count++;
      if (false)
        std::cout << __FUNCTION__ << "(SPH particle == " << a << " ): "
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
    pa.B = pa.rhoB_sub * pa.rho_weight;
    pa.S = pa.rhoS_sub * pa.rho_weight;
    pa.Q = pa.rhoQ_sub * pa.rho_weight;
    //cout << "AFTER: " << a << "   " << pa.B << "   "
    //    << pa.S << "   " << pa.Q << endl;
    //cout << pa.rho_weight << "   " << pa.rhoB_sub << "   "
    //    << pa.rhoS_sub << "   " << pa.rhoQ_sub << endl;
    //cout << "-----------------------------------------------------------------" << endl;
  }

  return;
}




///////////////////////////////////////////////////////////////////////////////
//Second smoothing smoothes the gradients after constructing all the  fields 
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

  for ( i.x[0] =- 2; i.x[0] <= 2; i.x[0]++ )
  for ( i.x[1] =- 2; i.x[1] <= 2; i.x[1]++ )
  {

    int b=systemPtr->linklist.lead[
            systemPtr->linklist.triToSum(
              systemPtr->linklist.dael[a] + i, systemPtr->linklist.size ) ];

    while( b != -1 )
    {
      const auto & pb          = systemPtr->particles[b];

      Vector<double,2> gradK   = kernel::gradKernel( pa.r - pb.r, settingsPtr->_h );
      Vector<double,2> va      = rowp1(0, pa.shv);
      Vector<double,2> vb      = rowp1(0, pb.shv);
      Matrix<double,2,2> vminia, vminib;

      mini(vminia, pa.shv);
      mini(vminib, pb.shv);

      double sigsqra           = 1.0/(pa.sigma*pa.sigma);
      double sigsqrb           = 1.0/(pb.sigma*pb.sigma);
      Vector<double,2> sigsigK = pb.sigmaweight * pa.sigma * gradK;

      pa.gradP                += ( sigsqrb*pb.eosPtr->p()
                                  + sigsqra*pa.eosPtr->p() ) * sigsigK;

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

      pa.gradshear            += inner(sigsigK, pa.v)*( sigsqrb*vb + sigsqra*va );
      pa.divshear             += sigsqrb*sigsigK*transpose(vminib)
                                  + sigsqra*sigsigK*transpose(vminia);

      if ( isnan( pa.gradP.x[0] ) )
      {
        cout << "gradP stopped working" << endl;
        cout << systemPtr->t <<" "  << pa.gradP << " " << a << " " << b << endl;
        cout << pb.sigmaweight << " " << pa.sigma << " " << pb.eosPtr->p() << endl;
        cout << systemPtr->linklist.Size << " " << pb.eosPtr->s() << " " << pa.eosPtr->s() << endl;

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
        && ( ( pa.eosPtr->T()*197.3 ) >= 150 ) )
    pa.frz2.t=tin;
  else if ( ( pa.btrack == 0 )
            && ( ( pa.eosPtr->T()*197.3 ) >= 150 )
            && ( pa.Freeze < 4 ) )
    cout << "Missed " << a << " " << tin << "  "
         << pa.eosPtr->T()*197.3 << " "
         << rdis << " " << settingsPtr->cfon <<  endl;

  return;
}
///////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::process_initial_conditions()
{
  // impose the energy cut-off before the initial time step of hydro
  // the original cut should be 0.15
  vector<Particle> threshold_particles;
  for (auto & p : systemPtr->particles)
  {
	  if (p.esub>0.00301) // this will be changed, NOT HARDCODED!!
	  {
		 threshold_particles.push_back(p);
	  }
  }
  systemPtr->particles = threshold_particles;


  // fill out initial particle information
	double stepX = settings->stepx;
	double stepY = settings->stepy;
	for (auto & p : systemPtr->particles)
  {
    //double stepx = 0.05, stepy = 0.05;
		p.u.x[0]          = 0.0;
		p.u.x[1]          = 0.0;
		p.eta_sigma       = 1.0;
		p.sigmaweight     = stepX*stepY;
		p.rho_weight      = stepX*stepY;
		p.Bulk            = 0.0;
		p.B               = p.rhoB_an*stepX*stepY;
		p.S               = p.rhoS_an*stepX*stepY;
		p.Q               = p.rhoQ_an*stepX*stepY;
		p.transverse_area = stepX*sstepY;

		//if (j==0)
		cout << "readICs_iccing(" << __LINE__ << "): "
		  << "SPH particles: "
		  << p.r.x[0] << "   " << p.r.x[1] << "   "
		  << p.e_sub << "   " << p.rhoB_an << "   "
		  << p.rhoS_an << "   " << p.rhoQ_an << "   "
		  << p.sigmaweight << endl;

		// make educated initial guess here for this particle's (T, mu_i) coordinates
		// (improve this in the future)
		p.thermo.T   = 500.0/hbarc_MeVfm;	// rootfinder seems to work better going downhill than "uphill"
		p.thermo.muB = 0.0/hbarc_MeVfm;
		p.thermo.muS = 0.0/hbarc_MeVfm;
		p.thermo.muQ = 0.0/hbarc_MeVfm;

		if (p.e_sub>systemPtr->efcheck)	// impose freeze-out check for e, not s
	    {
			p.Freeze=0;
		}
		else
		{
			p.Freeze=4;
			systemPtr->number_part += 1;
		}
  }
}
/////////////////////////////////////////////////////////////////////////////////
