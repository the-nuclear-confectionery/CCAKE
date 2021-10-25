#include "BSQHydro.h"


// Constructors and destructors.
  BSQHydro::BSQHydro(){}
 ~BSQHydro::BSQHydro(){}


void BSQHydro::load_settings_file( string path_to_settings_file )
{
  Input_Parameters input_parameters;
  io.load_settings_file(path_to_settings_file); // sets the settings path in
  // input_output, then loads parameters into Input_parameters struct
  io.set_EoS_type() // input_output talks to EoS and
  // tells it where to find its tables

  return
}



void BSQHydro::set_results_directory( string path_to_results_directory )
{
  io.set_results_directory(path_to_results_directory) //set the results directory in input_output
  return
}



void BSQHydro::read_in_initial_conditions()
{
  io.read_in_initial_conditions(system) // tells input_output to talk to system state and set initial system state
  return
}



////////////////////////////////////////////////////////////////////////////////
void BSQHydro::initialize_hydrodynamics()
{
  Imat.identity();  // need to set the identity matrix for use in simulation

  system.initialize();

  initialize_entropy_and_charge_densities(); // this should be in a switch/if

  initial_smoothing();

  return;
}


////////////////////////////////////////////////////////////////////////////////
void BSQHydro::run()
{
  system.BSQSimulation(); // dt and LinkList now members of system

  return;
}

void BSQHydro::find_freeze_out_surface(){}


void BSQHydro::print_results(){}


////////////////////////////////////////////////////////////////////////////////
void BSQHydro::initialize_entropy_and_charge_densities() // formerly updateIC
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
		cout << "----------------------------------------"
				"----------------------------------------" << endl;

		if (gtyp!=5)
		{
			sw.Start();
			cout << "Doing this particle: "
					<< _p[i].r.x[0] << "   " << _p[i].r.x[1] << "\n";
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].e_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;
			_p[i].s_an = _p[i].locate_phase_diagram_point_eBSQ(
                    _p[i].e_sub, _p[i].rhoB_an, _p[i].rhoS_an, _p[i].rhoQ_an );
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].e_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;

			if (true || VERBOSE>5)
			{
				if (_p[i].s_an>0.0)
				{
					double phase_diagram_point[4]
						= { _p[i].SPH_cell.T*197.3,
							  _p[i].SPH_cell.muB*197.3,
							  _p[i].SPH_cell.muS*197.3,
							  _p[i].SPH_cell.muQ*197.3 };
					double densities_at_point[4];
					get_eBSQ_densities(phase_diagram_point, densities_at_point);
					cout << i << ":   " << _p[i].e_sub*197.3
						<< "   " << _p[i].rhoB_an
						<< "   " << _p[i].rhoS_an
						<< "   " << _p[i].rhoQ_an
						<< "   " << _p[i].SPH_cell.T*197.3
						<< "   " << _p[i].SPH_cell.muB*197.3
						<< "   " << _p[i].SPH_cell.muS*197.3
						<< "   " << _p[i].SPH_cell.muQ*197.3;
						for (int iii = 0; iii < 4; iii++)
							cout << "   " << densities_at_point[iii];		
					cout << "\n";
				}
				else
					cout << i << ":   " << _p[i].e_sub*197.3
						<< "   " << _p[i].rhoB_an
						<< "   " << _p[i].rhoS_an
						<< "   " << _p[i].rhoQ_an
						<< "   " << 0.0 << "   " << 0.0
						<< "   " << 0.0 << "   " << 0.0
						<< "   nan   nan   nan   nan\n";
			}

			sw.Stop();
			string successString = (_p[i].s_an < 0.0) ?
									"unsuccessfully" : "successfully";
			cout << "SPH particle " << i << ", locate_phase_diagram_point_eBSQ: completed "
					<< successString << " in " << sw.printTime() << "s." << "\n";

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].e_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;
		}

		////////////////////////////////////////////////////////////////////////
		// for now, if we failed to find a real entropy density for this
		// point, just freeze it out, set its entropy to the freeze-out value,
		// and continue without setting anything else
		if (_p[i].s_an < 0.0)
		{

			////////////////////////////////////////////////////////
			// if failed with charge densities, set them to zero and re-solve;
			// if that fails too, guesstimate an answer
			cout << "\t --> Densities not found in EoS table (setting BSQ --> 0): "
					<< _p[i].r.x[0] << "   " << _p[i].r.x[1] << "\n"
					<< "\t\t - densities: "
					<< _p[i].e_sub*197.3 << "   " << _p[i].rhoB_an << "   "
					<< _p[i].rhoS_an << "   " << _p[i].rhoQ_an << "\n";

			// set charge densities to zero and re-solve
			_p[i].rhoB_an = 0.0;
			_p[i].rhoS_an = 0.0;
			_p[i].rhoQ_an = 0.0;

			_p[i].s_an = _p[i].locate_phase_diagram_point_eBSQ(
                    _p[i].e_sub, _p[i].rhoB_an, _p[i].rhoS_an, _p[i].rhoQ_an );

			// if this fails too...
			if (_p[i].s_an < 0.0)
			{
				double scale_factor = std::min( 1.0, _p[i].e_sub / efcheck );
	
				cout << "\t\t - scaling e to get s: "
						<< efcheck*0.1973 << "   "
						<< sfcheck << "   "
						<< scale_factor << "   "
						<< scale_factor * sfcheck << "\n";
	
				_p[i].s_an = scale_factor * sfcheck;
			}
			else	// if a solution was found
			{
				cout << "\t\t - phase diagram point: "
						<< _p[i].SPH_cell.T*197.3 << "   "
						<< _p[i].SPH_cell.muB*197.3 << "   "
						<< _p[i].SPH_cell.muS*197.3 << "   "
						<< _p[i].SPH_cell.muQ*197.3 << "\n";
			}

			// freeze this particle out!
			_p[i].Freeze = 4;
			number_part++;
			////////////////////////////////////////////////////////
		}
		else
		{
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].e_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;

			cout << "\t --> Densities found in EoS table: "
				<< _p[i].r.x[0] << "   " << _p[i].r.x[1] << "\n";
			cout << "\t\t - phase diagram point: "
					<< _p[i].SPH_cell.T*197.3 << "   "
					<< _p[i].SPH_cell.muB*197.3 << "   "
					<< _p[i].SPH_cell.muS*197.3 << "   "
					<< _p[i].SPH_cell.muQ*197.3 << "\n";
			cout << "\t\t - densities: "
					<< _p[i].e_sub*197.3 << "   "
					<< _p[i].rhoB_an << "   "
					<< _p[i].rhoS_an << "   "
					<< _p[i].rhoQ_an << "\n";
			
			cout << "\t --> Exact:\n";
			double phase_diagram_point[4] = { _p[i].SPH_cell.T*197.3,
											  _p[i].SPH_cell.muB*197.3,
											  _p[i].SPH_cell.muS*197.3,
											  _p[i].SPH_cell.muQ*197.3 };
			double densities_at_point[4];
			get_eBSQ_densities(phase_diagram_point, densities_at_point);
			cout << "\t\t - phase diagram point:";
			for (int iii = 0; iii < 4; iii++) cout << "   " << phase_diagram_point[iii];
			cout << "\n\t\t - densities:";
			for (int iii = 0; iii < 4; iii++) cout << "   " << densities_at_point[iii];
			cout << "\n";
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].e_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;

		}

    if (gtyp==5) _p[i].e_sub=_p[i].eos.e();

    _p[i].gamma=_p[i].gamcalc();

    _p[i].sigmaweight *= _p[i].s_an*_p[i].gamma*t0;	// sigmaweight is constant after this
    //_p[i].rho_weight *= _p[i].gamma*t0;				// rho_weight is constant after this

		_p[i].B *= _p[i].gamma*t0;	// B does not evolve in ideal case (confirm with Jaki)
		_p[i].S *= _p[i].gamma*t0;	// S does not evolve in ideal case (confirm with Jaki)
		_p[i].Q *= _p[i].gamma*t0;	// Q does not evolve in ideal case (confirm with Jaki)

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].e_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;

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
void BSQHydro::initial_smoothing()  // formerly BSQguess()
{
	cout << "setshear..." << endl;
	setshear();
	cout << "initiate..." << endl;
    initiate();

	cout << "bsqsvoptimization..." << endl;
	bool initialization_mode = true;
	for (int i=0; i<_n; i++)
	{
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].s_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;
		bsqsvoptimization(i, initialization_mode);
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].s_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;
	}
	cout << "One more loop!" << endl;

	int count1=0;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	for (int i=0; i<_n; i++)
	{
		_p[i].s_sub = _p[i].sigma/_p[i].gamma/t0;

if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].s_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;
		_p[i].locate_phase_diagram_point_sBSQ(
      _p[i].s_sub, _p[i].rhoB_sub, _p[i].rhoS_sub, _p[i].rhoQ_sub );
if (i==0)
	cout << "SPH checkpoint(" << __LINE__ << "): " << i << "   " << t << "   "
			<< _p[i].sigmaweight << "   " << _p[i].s_sub << "   "
			<< _p[i].eos.T() << "   " << _p[i].eos.e() << "   "
			<< _p[i].eos.p() << "   " << _p[i].s_an << endl;

		_p[i].sigsub = 0;
		_p[i].frzcheck(t0, count1, _n);
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	}
	return;
}