#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

#include "../include/constants.h"
#include "../include/formatted_output.h"
#include "../include/sph_workstation.h"
#include "../include/stopwatch.h"

using std::to_string;

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

  if ( VERBOSE > 5 )
  {
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
  }

    systemPtr->n_particles = systemPtr->particles.size();
    //cout << "systemPtr->n_particles = " << systemPtr->n_particles << "\n";

    for (int i=0; i<systemPtr->n_particles; i++)
    {
      auto & p = systemPtr->particles[i];


  if ( VERBOSE > 5 || p.print_this_particle )
  {
      cout << "----------------------------------------"
          "----------------------------------------" << "\n";
  }

		
		{
			sw.Start();
  if ( VERBOSE > 5 || p.print_this_particle )
  {
			cout << setprecision(12) << "Doing this particle: "
					<< p.r(0) << "   " << p.r(1) << "\n";
  }
      // solve for the entropy density
			p.input.s = locate_phase_diagram_point_eBSQ( p,
                    p.input.e, p.input.rhoB, p.input.rhoS, p.input.rhoQ );

			sw.Stop();
  if ( VERBOSE > 5 || p.print_this_particle )
  {
			string successString = (p.input.s < 0.0) ?
									"unsuccessfully" : "successfully";
      cout << "Print particle info:\n";
			cout << "    SPH particle " << i << ", locate_phase_diagram_point_eBSQ: completed "
					<< successString << " in " << sw.printTime() << "s." << "\n";
  }
		}

		///////////////////////////////////////////////////////////////////////////
		// for now, if we failed to find a real entropy density for this
		// point, just freeze it out, set its entropy to the freeze-out value,
		// and continue without setting anything else
		if (p.input.s < 0.0)
		{
			// freeze this particle out!
      cout << "This shouldn't have happened!" << endl;
      failCounter++;
      exit(8);
			p.Freeze = 4;
			////////////////////////////////////////////////////////
		}
		else
		{
  if ( VERBOSE > 5 || p.print_this_particle )
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
		}

    p.hydro.gamma = p.gamcalc();

    p.norm_spec.s    *= p.input.s*p.hydro.gamma*settingsPtr->t0;	  // constant after this
    p.norm_spec.rhoB *= p.input.rhoB*p.hydro.gamma*settingsPtr->t0; // constant after this
    p.norm_spec.rhoS *= p.input.rhoS*p.hydro.gamma*settingsPtr->t0; // constant after this
    p.norm_spec.rhoQ *= p.input.rhoQ*p.hydro.gamma*settingsPtr->t0; // constant after this

  if ( VERBOSE > 5 || p.print_this_particle )
  {
	cout << "----------------------------------------"
			"----------------------------------------" << "\n";
  }

	if (true)
	{
		cout << "Exiting prematurely from " << __PRETTY_FUNCTION__
			<< "::" << __LINE__ << "!" << endl;
		cerr << "Exiting prematurely from " << __PRETTY_FUNCTION__
			<< "::" << __LINE__ << "!" << endl;
		exit(8);
	}

    }

  if ( VERBOSE > 5 )
  {
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
  }

	swTotal.Stop();
	//cout << "Finished function call to " << __FUNCTION__ << "(...) in "
	//		<< swTotal.printTime() << " s." << endl;
  formatted_output::update("finished initializing particle densities in "
                              + to_string(swTotal.printTime()) + " s");

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
  if (VERBOSE > 5)
    cout << "----------------------------------------"
        "----------------------------------------" << endl;

	for ( auto & p : systemPtr->particles )
	{
    // must reset smoothed charge densities also
    double smoothed_s_lab    = p.hydro.sigma/p.hydro.gamma/settingsPtr->t0;
		double smoothed_rhoB_lab = p.smoothed.rhoB/p.hydro.gamma/settingsPtr->t0;
		double smoothed_rhoS_lab = p.smoothed.rhoS/p.hydro.gamma/settingsPtr->t0;
		double smoothed_rhoQ_lab = p.smoothed.rhoQ/p.hydro.gamma/settingsPtr->t0;

    // UNCOMMENT THIS AND DOCUMENT OUTPUT AS REFERENCE
    locate_phase_diagram_point_sBSQ( p, smoothed_s_lab, smoothed_rhoB_lab,
                                      smoothed_rhoS_lab, smoothed_rhoQ_lab );

		fo.check_freeze_out_status(p, settingsPtr->t0, count1, systemPtr->n_particles);
	}

	return;

}
///////////////////////////////////////////////////////////////////////////////
// smoothing routines: first smoothing covers all hydrodyanmical fields
void SPHWorkstation::smooth_fields(Particle & pa)
{
  int a = pa.ID;

  pa.hydro.sigma          = 0.0;
  pa.smoothed.s           = 0.0;
  pa.smoothed.rhoB        = 0.0;
  pa.smoothed.rhoS        = 0.0;
  pa.smoothed.rhoQ        = 0.0;

  auto & a_neighbors = systemPtr->linklist.all_neighbors[a];

  for ( int b : a_neighbors )
  {
    auto & pb       = systemPtr->particles[b];

    double kern     = kernel::kernel( pa.r - pb.r, settingsPtr->h );
    pa.hydro.sigma   += pb.norm_spec.s    * kern;
    pa.smoothed.s    += pb.norm_spec.s    * pb.specific.s    * kern;
    pa.smoothed.rhoB += pb.norm_spec.rhoB * pb.specific.rhoB * kern;
    pa.smoothed.rhoS += pb.norm_spec.rhoS * pb.specific.rhoS * kern;
    pa.smoothed.rhoQ += pb.norm_spec.rhoQ * pb.specific.rhoQ * kern;

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
                << "   " << pb.norm_spec.s
                << "   " << pb.specific.s
                << "   " << pb.input.rhoB
                << "   " << pa.smoothed.rhoB
                << "   " << pb.input.rhoS
                << "   " << pa.smoothed.rhoS
                << "   " << pb.input.rhoQ
                << "   " << pa.smoothed.rhoQ
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
    Vector<double,2> gradK   = kernel::gradKernel( rel_sep, rel_sep_norm, settingsPtr->h );
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

    double relative_distance_by_h = rel_sep_norm / settingsPtr->h;
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
      cout << kernel::kernel( pa.r - pb.r, settingsPtr->h ) << endl;
    }
    else if ( isnan( pah.gradP(1) ) )
      cout << "1 " << gradPressure_weight(a, b)
           << " " << a << " " << b << endl;
  }

  const double hc = constants::hbarc_MeVfm;

  if ( ( pa.btrack == 1 )                               // if particle a has only
            && ( ( pa.T()*hc ) >= 150 ) )               // one nearest neighbor (and T>150),
    fo.frz2[pa.ID].t=tin;                               // set penultimate timestep;
  else if ( ( pa.btrack == 0 )                          // otherwise, if a has no nearest neighbors
            && ( ( pa.T()*hc ) >= 150 )                 // but has T>150 and isn't frozen out,
            && ( pa.Freeze < 4 ) )                      // just print this warning message
    cout << boolalpha << "Missed " << a << " "
         << tin << " " << pa.T()*hc << " "
         << rdis << " " << systemPtr->do_freeze_out << endl;

  return;
}




////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::process_initial_conditions()
{
  formatted_output::report("Processing initial conditions");

  //============================================================================
  // IMPOSE ENERGY/CHARGE CUTOFFS TO REGULATE EVENT (NO CUTOFFS FOR GUBSER)
  if ( settingsPtr->IC_type != "Gubser"
        && settingsPtr->IC_type != "Gubser_with_shear")
  {

    formatted_output::update("Input number of particles: "
                              + to_string(systemPtr->particles.size()));

    //==========================================================================
    // impose the energy cut-off before the initial time step of hydro
    systemPtr->particles.erase( std::remove_if(
      systemPtr->particles.begin(),
      systemPtr->particles.end(),
      [this](Particle const & p){ return p.input.e <= settingsPtr->e_cutoff; }),
      systemPtr->particles.end() );



    formatted_output::update("Number of particles after e-cutoff: "
                              + to_string(systemPtr->particles.size()));



//    //==========================================================================
//    // add a buffer of particles around the edge to stabilize evolution in low-
//    // density regions
//    if ( settingsPtr->buffer_event ) add_buffer();
//
//
//    formatted_output::update("Number of particles after buffering: "
//                              + to_string(systemPtr->particles.size()));
//
//
//  //if (1) exit(8);



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


    formatted_output::update("Number of particles after solution checks: "
                              + to_string(systemPtr->particles.size()));
  }

  // fill out initial particle information
  //int TMP_particle_count = 0;
	for (auto & p : systemPtr->particles)
  {
    // set area element for each SPH particle
    double dA = (settingsPtr->stepx)*(settingsPtr->stepy);

    // Set the rest of particle elements using area element
		//p.u(0)          = 0.0;  // flow must be set in Particle constructor!!!
		//p.u(1)          = 0.0;  // flow must be set in Particle constructor!!!

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
		p.thermo.T        = 1000.0/hbarc_MeVfm;	// rootfinder seems to work better going downhill than "uphill"
		p.thermo.muB      = 0.0/hbarc_MeVfm;
		p.thermo.muS      = 0.0/hbarc_MeVfm;
		p.thermo.muQ      = 0.0/hbarc_MeVfm;
		p.thermo.eos_name = "default";  // uses whatever the default EoS is

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
                           + to_string(systemPtr->particles.size()
                                        - systemPtr->number_part) );


  //============================================================================
  // with particles vector now fully initialized, specify or initialize any
  // remaining quantities which depend on this

  // assign particles IDs
  for ( int i = 0; i < systemPtr->particles.size(); i++ )
  {
    systemPtr->particles[i].ID       = i;
    systemPtr->particles[i].hydro.ID = i;
  }

  // set particles to print
  for ( int & p : settingsPtr->particles_to_print )
  {
    systemPtr->particles[p].print_this_particle  = true;
    systemPtr->particles[p].hydro.print_particle = true;
  }
}


//==============================================================================
// initialize bulk Pi 
void SPHWorkstation::set_bulk_Pi()
{
  double t0 = settingsPtr->t0;
  if ( settingsPtr->initializing_with_full_Tmunu )
    for ( auto & p : systemPtr->particles )
      p.hydro.Bulk = ( p.hydro.varsigma - p.p() )
                      * p.hydro.gamma * t0 / p.hydro.sigma;

  return;
}


//==============================================================================
void SPHWorkstation::freeze_out_particles()
{
  //---------------------------------------
  // perform freeze out checks
  int n_freezing_out = 0;
  for ( auto & p : systemPtr->particles )
    fo.check_freeze_out_status( p, systemPtr->t, n_freezing_out,
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
  fo.bsqsvfreezeout( n_freezing_out );
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










////////////////////////////////////////////////////////////////////////////////
double SPHWorkstation::locate_phase_diagram_point_eBSQ( Particle & p,
                 double e_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
  // default: use particle's current location as initial guess
  // (pass in corresponding EoS as well!)
  eos.tbqs( p.T(), p.muB(), p.muQ(), p.muS(), p.get_current_eos_name() );

  bool solution_found = false;
  double sVal = eos.s_out( e_In, rhoB_In, rhoS_In, rhoQ_In, solution_found );

  if ( solution_found )
{
    eos.set_thermo( p.thermo );

if (p.print_this_particle)
  cout << "check thermo: " << systemPtr->t << "   "
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
  if (p.thermo.cs2<0)
  {
    cout << __LINE__ << "cs2 was negative!" << endl;
    exit(8);
  }
}

  return sVal;
}

////////////////////////////////////////////////////////////////////////////////
double SPHWorkstation::locate_phase_diagram_point_eBSQ(Particle & p, double e_In)
                 { return locate_phase_diagram_point_eBSQ( p, e_In, 0.0, 0.0, 0.0 ); }






////////////////////////////////////////////////////////////////////////////////
void SPHWorkstation::locate_phase_diagram_point_sBSQ( Particle & p,
                 double s_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
//  cout << "Rootfinder for p.ID = " << p.ID << endl;

  // default: use particle's current location as initial guess
  eos.tbqs( p.T(), p.muB(), p.muQ(), p.muS(), p.get_current_eos_name() );

  bool update_s_success = eos.update_s( s_In, rhoB_In, rhoS_In, rhoQ_In );

  if ( update_s_success )
{
    eos.set_thermo( p.thermo );

if (p.print_this_particle)
  cout << "check thermo: " << systemPtr->t << "   "
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
  if (p.thermo.cs2<0)
  {
    cout << __LINE__ << "cs2 was negative!" << endl;
    exit(8);
  }
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
  double rhoB_lab   = p.smoothed.rhoB/p.hydro.gamma/systemPtr->t;
  double rhoS_lab   = p.smoothed.rhoS/p.hydro.gamma/systemPtr->t;
  double rhoQ_lab   = p.smoothed.rhoQ/p.hydro.gamma/systemPtr->t;
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

  return pb.norm_spec.s * pa.hydro.sigma
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


