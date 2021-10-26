#include "system_state.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

//using namespace std;
using std::cout;
using std::endl;
using std::string;

#include "constants.h"
#include "vector.h"
#include "particle.h"
#include "runge_kutta.h"
#include "eos.h"

using namespace constants;


////////////////////////////////////////////////////////////////////////////////
void SystemState::initialize()  // formerly called "manualenter"
{
  double h, factor;
  double it0;
  int start, end;

  int df;

  linklist.setv( fvisc );
  linklist.eost   = eostype;
  linklist.cevent = 0;
  std::cout << fvisc << " hydro, h=" << h <<  " dimensions=" << D
            << " dt=" << ics.dt << " QM fluc:  " << linklist.qmf << "\n";

  //////////////////////////////////////////////////////////////////////////////
  // SET EQUATION OF STATE
  // rewrite by C. Plumberg: allow for different EOS format if using BSQ
  double efcheck = 0.0, sfcheck = 0.0;
  eos EOS0;	// now declared globally
  if ( linklist.visc == 4 )	//if we're running BSQ (table is only option)
  {
    bool using_HDF = false;
    if (using_HDF)
    {
      string quantityFile   = ifolder + std::string("quantityFile.h5");
      string derivativeFile = ifolder + std::string("derivFile.h5");
      std::cout << "Using BSQ Equation of State table from: "
                << quantityFile << " and " << derivativeFile << "\n";

      EOS0.init( quantityFile, derivativeFile );
    }
    else
    {
      string quantityFilename   = "EoS_Taylor_AllMu_T0_1200.dat";
      string derivativeFilename = "EoS_Taylor_AllMu_Derivatives_T0_1200.dat";
      string quantityFile       = ifolder + quantityFilename;
      string derivativeFile     = ifolder + derivativeFilename;
      std::cout << "Using BSQ Equation of State table from: "
                << quantityFile << " and " << derivativeFile << "\n";

      EOS0.init( quantityFile, derivativeFile );
    }

    EOS0.eosin( eostype );			// does nothing!
    const double freeze_out_T_at_mu_eq_0
                  = 0.15/hbarc_GeVfm;	//1/fm
    efcheck       = EOS0.efreeze( freeze_out_T_at_mu_eq_0 );
    sfcheck       = EOS0.sfreeze( freeze_out_T_at_mu_eq_0 );
    //efcheck = 0.266112/0.1973;
    //sfcheck = 2.05743;

    std::cout << "efcheck = " << efcheck*hbarc_GeVfm << " GeV/fm^3\n";
    std::cout << "sfcheck = " << sfcheck << " 1/fm^3\n";
  }
  else
  {
    std::cerr << "This EoS model not currently supported!" << std::endl;
  }

  linklist.efcheck = efcheck;
  linklist.sfcheck = sfcheck;
  linklist.fcount  = 0;
  linklist.average = 0;
  //       Start reading ICs          //

  int numpart, _Ntable3;

  //  cout << "setting up SPH" << endl;

  cout << "Initial conditions type: " << ictype << endl;

  linklist.gtyp=0;
  if ( ictype == iccing )
  {

    int count           = 1;
    linklist.ebe_folder = outf;
    vector<string>        filelist( count );

    int j               = 0;
    filelist[j]         = ic + "/ic0.dat"; // only doing single event
    linklist.filenames  = filelist;
    linklist.fcount     = count;
    linklist.fnum       = linklist.start;

    // already done
    //readICs_iccing(linklist.filenames[0], _Ntable3, _p, factor, efcheck, numpart, EOS0);

    ////////////////////////////////////////////////////////////////////////////
    //assign thermodynamic quatities to sph particles
    //Particle::set_equation_of_state( EOS0 );
    Particle::set_equation_of_state( &EOS0 );

    // initialize 0th particle
    //_p[0].start(eostype, EOS0);

    // assume initial conditions have been read in from file
    
    linklist.initialize( it0, _Ntable3, h, particles, ics.dt, numpart );

    cout << "number of sph particles=" << _Ntable3 << endl;
    linklist.gtyp=6;

  }

  /*if ( ictype == iccing )
  {
    linklist.updateIC();
    cout << "bsq optimization done" << endl;
    linklist.bsqsvfreezeset();
  }*/

  // formerly bsqsv_set in this loop
  for (auto & p : particles)
  {
    double gg = p.gamcalc();
    p.g2      = gg*gg;
    p.shv33   = 0.0;
  }

  return;
}
///////////////////////////////////////////////////////////////////////////////
//Dekra: BSQsimulation was moved from this location to BSQhydro
///////////////////////////////////////////////////////////////////////////////
//Start routines for checking conservation of energy density, entropy density and BSQ
// charges at each time step of the simulation
///////////////////////////////////////
void SystemState::check_BSQ_energy_conservation()
{
  E=0.0;
  for (int i=0; i<_n; i++)
    E += ( _p[i].C* _p[i].g2 - _p[i].EOSp() - _p[i].bigPI + _p[i].shv.x[0][0] )
          *_p[i].sigmaweight*t/_p[i].sigma;

  if (first == 1)
  {
    first = 0;
    E0    = E;
  }

  return;
}
////////////////////////////////////////
void SystemState::check_BSQ_charge_conservation()
{
  Btotal = 0.0;
  Stotal = 0.0;
  Qtotal = 0.0;

  for (int i=0; i<_n; i++)
  {
    //Btotal += _p[i].B;
    //Stotal += _p[i].S;
    //Qtotal += _p[i].Q;
    Btotal += _p[i].rhoB_sub*_p[i].rho_weight;
    Stotal += _p[i].rhoS_sub*_p[i].rho_weight;
    Qtotal += _p[i].rhoQ_sub*_p[i].rho_weight;
  }

  if (first==1)
  {
    Btotal0 = Btotal;
    Stotal0 = Stotal;
    Qtotal0 = Qtotal;
  }

	return;
}
///////////////////////////////////////
void SystemState::bsqsvconservation()
{
    bsqsvconservation_E();
    Etot  = E + Ez;
    Eloss = (E0-Etot)/E0*100;
    rk2   = 0;
}
///////////////////////////////////////
void SystemState::conservation_entropy()
{
  S=0.0;

  for (int i=0; i<_n; i++)
  {
    S += particles[i].eta_sigma*particles[i].sigmaweight;
    if (i==0)
    std::cout << "\t\t --> " << i << "   " << particles[i].eta_sigma << "   "
              << particles[i].sigmaweight << "   " << S << endl;
  }

  if (first==1)
    S0=S;
}
///////////////////////////////////////
// function to check conservation of B, S, and Q
void SystemState::conservation_BSQ()
{
    Btotal = 0.0;
    Stotal = 0.0;
    Qtotal = 0.0;

    for (int i=0; i<_n; i++)
	{
        //Btotal += particles[i].B;
        //Stotal += particles[i].S;
        //Qtotal += particles[i].Q;
        Btotal += particles[i].rhoB_sub*particles[i].rho_weight;
        Stotal += particles[i].rhoS_sub*particles[i].rho_weight;
        Qtotal += particles[i].rhoQ_sub*particles[i].rho_weight;
    }

    if (first==1)
    {
        Btotal0 = Btotal;
        Stotal0 = Stotal;
        Qtotal0 = Qtotal;
    }
	return;
}
///////////////////////////////////////
void SystemState::bsqsvconservation_E()
{

    E=0.;
    for (int i=0; i<_n; i++)
    {
      const auto & p = particles[i];

        E += ( p.C*p.g2 - p.eosPtr->p() - p.bigPI + p.shv.x[0][0] )
              / p.sigma*p.sigmaweight*t;
        if (i==0)
          std::cout << "E: " << i << "   " << t
              << "   " << p.eosPtr->T()
              << "   " << p.EOSe()
              << "   " << p.C
              << "   " << p.g2
              << "   " << p.eosPtr->p()
              << "   " << p.bigPI
              << "   " << p.shv.x[0][0]
              << "   " << p.sigma
              << "   " << p.sigmaweight << endl;    }

    if (first==1)
    {
      first=0;
      E0=E;
    }
}
///////////////////////////////////////
void SystemState::bsqsvconservation_Ez()
{
  dEz=0.;

  double t2=t*t;
  for (int i=0; i<_n; i++)
  {
    const auto & p = particles[i];

    dEz += ( p.eosPtr->p() + p.bigPI + p.shv33*t2 ) / p.sigma*p.sigmaweight;

    if (false)
      std::cout << "dEz: " << i << "   " << t
        << "   " << p.eosPtr->p()
        << "   " << p.bigPI
        << "   " << p.shv33*t2
        << "   " << p.sigma
        << "   " << p.sigmaweight << endl;
  }
}
// first smoothing routine covers all hydrodyanmical fields
void SystemState::smooth_fields(int a, bool init_mode /*== false*/)
{
  const auto & pa    = _p[a];
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
    int b = lead[ triToSum( dael[a] + i, size ) ];
    while ( b != -1 )
    {
      const auto & pb = _p[b];
      double kern     = kernel( pa.r - pb.r );
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

      b = link[b];
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
//Second smoothing smoothes the gradients after reading in the updated fields 
//from the equation of state
void SystemState::smooth_gradients( int a, double tin, int & count )
{
  const auto & pa    = _p[a];

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

    int b=lead[ triToSum( dael[a] + i, size ) ];

    while( b != -1 )
    {
      const auto & pb          = _p[b];

      Vector<double,2> gradK   = gradKernel( pa.r - pb.r,
                                  static_cast<bool>( a == 30 && b == 43 ) );
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

      if ( ( ( Norm( pa.r - pb.r ) / _h ) <= 2 ) && ( a != b ) )
      {
        if ( pa.btrack != -1 ) pa.btrack++;
        if ( pa.btrack ==  1 ) rdis = Norm(pa.r-pb.r)/_h;
      }

      pa.gradBulk             += ( pb.Bulk/pb.sigma/pb.gamma
                                    + pa.Bulk/pa.sigma/pa.gamma)/tin*sigsigK;
      //pa.gradrhoB             += ( pb.rhoB/pb.sigma/pb.gamma
      //                            + pa.rhoB/pa.sigma/pa.gamma)/tin*sigsigK;
      //pa.gradrhoS             += ( pb.rhoS/pb.sigma/pb.gamma
      //                            + pa.rhoS/pa.sigma/pa.gamma)/tin*sigsigK;
      //pa.gradrhoQ             += ( pb.rhoQ/pb.sigma/pb.gamma
      //                            + pa.rhoQ/pa.sigma/pa.gamma)/tin*sigsigK;
      pa.gradV                += pb.sigmaweight*( pb.v -  pa.v )*gradK/pa.sigma;

      pa.gradshear            += inner(sigsigK, pa.v)*( sigsqrb*vb + sigsqra*va );
      pa.divshear             += sigsqrb*sigsigK*transpose(vminib)
                                  + sigsqra*sigsigK*transpose(vminia);

      if ( isnan( pa.gradP.x[0] ) )
      {
        cout << "gradP stopped working" << endl;
        cout << t <<" "  << pa.gradP << " " << a << " " << b << endl;
        cout << pb.sigmaweight << " " << pa.sigma << " " << pb.eosPtr->p() << endl;
        cout << Size << " " << pb.eosPtr->s() << " " << pa.eosPtr->s() << endl;

        cout << pa.r << endl;
        cout << pb.r << endl;
        cout << kernel( pa.r - pb.r ) << endl;
      }
      else if ( isnan( pa.gradP.x[1] ) )
        cout << "1 " << gradPressure_weight(pa, pb)
             << " " << a << " " << b << endl;
      else if ( isnan( pa.gradP.x[2] ) )
        cout << "2 " << gradPressure_weight(pa, pb)
             << " " << a << " " << b << endl;

      b=link[b];
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
         << rdis << " " << cfon <<  endl;

  return;
}
///////////////////////////////////////
void SystemState::setshear()
{
    for ( auto & p : particles ) p.sets(t*t);
}
///////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////
void SystemState::initialize_entropy_and_charge_densities() // formerly updateIC
{
	// set up EoS C library
	initialize("/projects/jnorhos/BSQ/EoS_BQS_Derivatives/Coefficients_Parameters.dat");
	Stopwatch sw, swTotal;
	swTotal.Start();
	long long failCounter = 0;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
    for (int i=0; i<_n; i++)
    {
      auto & p = particles[i];



		cout << "----------------------------------------"
				"----------------------------------------" << endl;

		if (settings.gtyp!=5)
		{
			sw.Start();
			cout << "Doing this particle: "
					<< p.r.x[0] << "   " << p.r.x[1] << "\n";
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
			p.s_an = p.locate_phase_diagram_point_eBSQ(
                    p.e_sub, p.rhoB_an, p.rhoS_an, p.rhoQ_an );
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;

			if (true || settings.VERBOSE>5)
			{
				if (p.s_an>0.0)
				{
					double phase_diagram_point[4]
						= { p.SPH_cell.T*197.3,
							  p.SPH_cell.muB*197.3,
							  p.SPH_cell.muS*197.3,
							  p.SPH_cell.muQ*197.3 };
					double densities_at_point[4];
					get_eBSQ_densities(phase_diagram_point, densities_at_point);
					cout << i << ":   " << p.e_sub*197.3
						<< "   " << p.rhoB_an
						<< "   " << p.rhoS_an
						<< "   " << p.rhoQ_an
						<< "   " << p.SPH_cell.T*197.3
						<< "   " << p.SPH_cell.muB*197.3
						<< "   " << p.SPH_cell.muS*197.3
						<< "   " << p.SPH_cell.muQ*197.3;
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
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
		}

		////////////////////////////////////////////////////////////////////////
		// for now, if we failed to find a real entropy density for this
		// point, just freeze it out, set its entropy to the freeze-out value,
		// and continue without setting anything else
		if (p.s_an < 0.0)
		{

			////////////////////////////////////////////////////////
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
				double scale_factor = std::min( 1.0, p.e_sub / efcheck );
	
				cout << "\t\t - scaling e to get s: "
						<< efcheck*0.1973 << "   "
						<< sfcheck << "   "
						<< scale_factor << "   "
						<< scale_factor * sfcheck << "\n";
	
				p.s_an = scale_factor * sfcheck;
			}
			else	// if a solution was found
			{
				cout << "\t\t - phase diagram point: "
						<< p.SPH_cell.T*197.3 << "   "
						<< p.SPH_cell.muB*197.3 << "   "
						<< p.SPH_cell.muS*197.3 << "   "
						<< p.SPH_cell.muQ*197.3 << "\n";
			}

			// freeze this particle out!
			p.Freeze = 4;
			number_part++;
			////////////////////////////////////////////////////////
		}
		else
		{
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;

			cout << "\t --> Densities found in EoS table: "
				<< p.r.x[0] << "   " << p.r.x[1] << "\n";
			cout << "\t\t - phase diagram point: "
					<< p.SPH_cell.T*197.3 << "   "
					<< p.SPH_cell.muB*197.3 << "   "
					<< p.SPH_cell.muS*197.3 << "   "
					<< p.SPH_cell.muQ*197.3 << "\n";
			cout << "\t\t - densities: "
					<< p.e_sub*197.3 << "   "
					<< p.rhoB_an << "   "
					<< p.rhoS_an << "   "
					<< p.rhoQ_an << "\n";
			
			cout << "\t --> Exact:\n";
			double phase_diagram_point[4] = { p.SPH_cell.T*197.3,
											  p.SPH_cell.muB*197.3,
											  p.SPH_cell.muS*197.3,
											  p.SPH_cell.muQ*197.3 };
			double densities_at_point[4];
			get_eBSQ_densities(phase_diagram_point, densities_at_point);
			cout << "\t\t - phase diagram point:";
			for (int iii = 0; iii < 4; iii++) cout << "   " << phase_diagram_point[iii];
			cout << "\n\t\t - densities:";
			for (int iii = 0; iii < 4; iii++) cout << "   " << densities_at_point[iii];
			cout << "\n";
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< p.sigmaweight << "   " << p.e_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;

		}

    if (settings.gtyp==5) p.e_sub = p.eosPtr->e();

    p.gamma=p.gamcalc();

    p.sigmaweight *= p.s_an*p.gamma*t0;	// sigmaweight is constant after this
    //p.rho_weight *= p.gamma*t0;				// rho_weight is constant after this

		p.B *= p.gamma*t0;	// B does not evolve in ideal case (confirm with Jaki)
		p.S *= p.gamma*t0;	// S does not evolve in ideal case (confirm with Jaki)
		p.Q *= p.gamma*t0;	// Q does not evolve in ideal case (confirm with Jaki)

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
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
void SystemState::initial_smoothing()  // formerly BSQguess()
{
	cout << "setshear..." << endl;
  setshear();
	cout << "initiate..." << endl;
  //initiate();

	cout << "bsqsvoptimization..." << endl;
	bool initialization_mode = true;
	for (int i=0; i<_n; i++)
	{
    auto & p = particles[i];

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< p.sigmaweight << "   " << p.s_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
		smooth_fields(i, initialization_mode);
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< p.sigmaweight << "   " << p.s_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
	}
	cout << "One more loop!" << endl;

	int count1=0;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	for (int i=0; i<_n; i++)
	{
		p.s_sub = p.sigma/p.gamma/settings.t0;

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< p.sigmaweight << "   " << p.s_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;
		p.locate_phase_diagram_point_sBSQ(
      p.s_sub, p.rhoB_sub, p.rhoS_sub, p.rhoQ_sub );
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< p.sigmaweight << "   " << p.s_sub << "   "
			<< p.eosPtr->T() << "   " << p.eosPtr->e() << "   "
			<< p.eosPtr->p() << "   " << p.s_an << endl;

		p.sigsub = 0;
		p.frzcheck(t0, count1, _n);
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	}
	return;




  cout << "BSQ-SV simulation completed!" << endl;

  return;
}




void SystemState::set_current_timestep_quantities();
{
  for (int i=0; i<N; ++i)
  {
    const auto & p = particles[i];
    u0[i]        = p.u;
    r0[i]        = p.r;
    etasigma0[i] = p.eta_sigma;
    Bulk0[i]     = p.Bulk;
    mini( shv0[i], p.shv );
  }
}

void SystemState::set_current_timestep_quantities()
{
  etasigma0.resize(N, 0.0);
  Bulk0.resize(N, 0.0);

  u0.resize(N, 0.0);
  r0.resize(N, 0.0, 0.0);

  shv0.resize(N, 0.0);

  for (int i=0; i<N; ++i)
  {
    const auto & p = particles[i];
    u0[i]        = p.u;
    r0[i]        = p.r;
    etasigma0[i] = p.eta_sigma;
    Bulk0[i]     = p.Bulk;
    mini( shv0[i], p.shv );
  }
}

void SystemState::get_derivative_halfstep()
{
  for (int i=0; i<N; ++i)
  {
    const auto & p = particles[i];
    p.u            = u0[i]        + 0.5*dx*p.du_dt;
    p.r            = r0[i]        + 0.5*dx*p.v;
    p.eta_sigma    = etasigma0[i] + 0.5*dx*p.detasigma_dt;
    p.Bulk         = Bulk0[i]     + 0.5*dx*p.dBulk_dt;
    tmini( p.shv,    shv0[i]      + 0.5*dx*p.dshv_dt );
  }
}


void SystemState::get_derivative_fullstep()
{
  for (int i=0; i<N; ++i)
  {
    const auto & p = particles[i];
    p.u            = u0[i]        + dx*p.du_dt;
    p.r            = r0[i]        + dx*p.v;
    p.eta_sigma    = etasigma0[i] + dx*p.detasigma_dt;
    p.Bulk         = Bulk0[i]     + dx*p.dBulk_dt;
    tmini( p.shv,    shv0[i]      + dx*p.dshv_dt );
  }
}
