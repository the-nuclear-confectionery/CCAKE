#include <functional>

#include "../include/rootfinder.h"

////////////////////////////////////////////////////////////////////////////////
//struct to pass the target (E, rhoB, rhoQ, rhoS) into the rootfinder function
struct rootfinder_parameters
{
  int e_or_entr_mode;
  double eorEntGiven;          //these are the desired e/s and BSQ
  double rhoBGiven;
  double rhoQGiven;
  double rhoSGiven;

  // this function should take (T,muX) and return (e/s,rhoX)
  std::function<void(double[], double[])> f;

  rootfinder_parameters();
  rootfinder_parameters( double setEorEntGiven, double setRhoBGiven,
                         double setRhoQGiven, double setRhoSGiven,
                         int set_e_or_entr_mode,
                         std::function<void(double[], double[])> f_in );
};

rootfinder_parameters::rootfinder_parameters() {}
rootfinder_parameters::rootfinder_parameters(
  double setEorEntGiven, double setRhoBGiven, double setRhoQGiven,
  double setRhoSGiven, int set_e_or_entr_mode,
  std::function<void(double[], double[])> f_in )
{
  e_or_entr_mode = set_e_or_entr_mode;
  eorEntGiven    = setEorEntGiven;
  rhoBGiven      = setRhoBGiven;
  rhoQGiven      = setRhoQGiven;
  rhoSGiven      = setRhoSGiven;

  f              = f_in;
}



////////////////////////////////////////////////////////////////////////////////
int rootfinder_f(const gsl_vector *x, void *params, gsl_vector *f)
{
cout << "Entered at line = " << __LINE__ << endl;
    //x contains the next (T, muB, muS) coordinate to test
    vector<double> tbqsToEval(4);
    tbqsToEval[0] = gsl_vector_get(x,0);
    tbqsToEval[1] = gsl_vector_get(x,1);	// convert x into densevector so it
    tbqsToEval[2] = gsl_vector_get(x,2);	// can be a BSpline evaluation point
    tbqsToEval[3] = gsl_vector_get(x,3);

    int e_or_entr_mode;
    double eorEntGiven, rhoBGiven, rhoQGiven, rhoSGiven, eorEnt, rhoB, rhoQ, rhoS;
    e_or_entr_mode  = ((rootfinder_parameters*)params)->e_or_entr_mode;
    eorEntGiven     = ((rootfinder_parameters*)params)->eorEntGiven;
    rhoBGiven       = ((rootfinder_parameters*)params)->rhoBGiven;
    rhoQGiven       = ((rootfinder_parameters*)params)->rhoQGiven;
    rhoSGiven       = ((rootfinder_parameters*)params)->rhoSGiven;
    std::function<void(double[], double[])>
      get_densities = ((rootfinder_parameters*)params)->f;

    // limit scope for readability
    {
      double phase_diagram_point[4]
              = { tbqsToEval[0], tbqsToEval[1], tbqsToEval[2], tbqsToEval[3] };
      double densities_at_point[4] = {-1.0, -1.0, -1.0, -1.0};

      // compute densities using passed-in function object
      get_densities( phase_diagram_point, densities_at_point );

cout << "PD point:";
for (int i = 0; i < 4; i++) cout << "   " << phase_diagram_point[i];
cout << endl << "Densities:";
for (int i = 0; i < 4; i++) cout << "   " << densities_at_point[i];
cout << endl;
//if (1) exit(1);

      // set densities (convert to powers of fm if necessary)
      eorEnt  = densities_at_point[0];
      rhoB    = densities_at_point[1];
      rhoS    = densities_at_point[2];
      rhoQ    = densities_at_point[3];
    }

    // set differences from zero
    gsl_vector_set(f, 0, (eorEnt - eorEntGiven));
    gsl_vector_set(f, 1, (rhoB   - rhoBGiven));
    gsl_vector_set(f, 2, (rhoQ   - rhoQGiven));
    gsl_vector_set(f, 3, (rhoS   - rhoSGiven));

cout << "e: " << eorEnt << "   " << eorEntGiven << "   " << eorEnt - eorEntGiven << endl;
cout << "B: " << rhoB << "   " << rhoBGiven << "   " << rhoB - rhoBGiven << endl;
cout << "Q: " << rhoQ << "   " << rhoQGiven << "   " << rhoQ - rhoQGiven << endl;
cout << "S: " << rhoS << "   " << rhoSGiven << "   " << rhoS - rhoSGiven << endl
      <<"------" << endl;
//if (1) exit(1);


cout << "Exited at line = " << __LINE__ << endl;
    return GSL_SUCCESS;
}




////////////////////////////////////////////////////////////////////////////////
// ROOTFINDER METHODS BELOW THIS LINE
void Rootfinder::tbqs( vector<double> & tbqsIn )
{
  tbqs( tbqsIn[0], tbqsIn[1], tbqsIn[2], tbqsIn[3] );
}


void Rootfinder::tbqs(double setT, double setmuB, double setmuQ, double setmuS)
{
  if(setT < minT || setT > maxT) {
    std::cout << "T = " << setT
      << " is out of range. Valid values are between ["
      << minT << "," << maxT << "]\n";
    return;
  }
  if(setmuB < minMuB || setmuB > maxMuB) {
    std::cout << "muB = " << setmuB
      << " is out of range. Valid values are between ["
      << minMuB << "," << maxMuB << "]\n";
    return;
  }
  if(setmuQ < minMuQ || setmuQ > maxMuQ) {
    std::cout << "muQ = " << setmuQ
      << " is out of range. Valid values are between ["
      << minMuQ << "," << maxMuQ << "]\n";
    return;
  }
  if(setmuS < minMuS || setmuS > maxMuS) {
    std::cout << "muS = " << setmuS
      << " is out of range. Valid values are between ["
      << minMuS << "," << maxMuS << "]\n";
    return;
  }

	tbqsPosition[0] = setT;
	tbqsPosition[1] = setmuB;
	tbqsPosition[2] = setmuQ;
	tbqsPosition[3] = setmuS;

}


////////////////////////////////////////////////////////////////////////////////
bool Rootfinder::rootfinder4D(double e_or_s_Given, int e_or_s_mode,
						double rhoBGiven, double rhoSGiven, double rhoQGiven,
						double error, size_t steps,
            std::function<void(double[], double[])> function_to_evaluate,
            vector<double> & updated_tbqs )
{
  /////////////////////////////////////////////
  // e_or_s_mode == 0: using entropy density //
  // e_or_s_mode == 1: using energy density  //
  /////////////////////////////////////////////

  std::cout << "Starting new rootfinder4D call\n";

  ////////////////////
  // set initial guess
  gsl_vector *x = gsl_vector_alloc(4);
  for (int iTBQS = 0; iTBQS < 4; iTBQS++)
    gsl_vector_set(x, iTBQS, tbqsPosition[iTBQS]);

  gsl_vector *chosen_densities = gsl_vector_alloc(4);
  gsl_vector_set(chosen_densities, 0, e_or_s_Given);
  gsl_vector_set(chosen_densities, 1, rhoBGiven);
  gsl_vector_set(chosen_densities, 2, rhoQGiven);
  gsl_vector_set(chosen_densities, 3, rhoSGiven);


  ////////////////////
  // pass relevant parameters to rootfinder
  rootfinder_parameters p( e_or_s_Given, rhoBGiven, rhoQGiven, rhoSGiven,
                           e_or_s_mode, function_to_evaluate );

  std::cout << __LINE__ << ": " << e_or_s_Given << std::endl;
  std::cout << __LINE__ << ": " << rhoBGiven << std::endl;
  std::cout << __LINE__ << ": " << rhoSGiven << std::endl;
  std::cout << __LINE__ << ": " << rhoQGiven << std::endl;
  std::cout << __LINE__ << ": " << e_or_s_mode << std::endl;

  ////////////////////
  // initialize multiroot solver
  gsl_multiroot_fsolver *solver;
  gsl_multiroot_function f = {&rootfinder_f, 4, &p};

  //f.n      = 4;
  //f.params = &p;
  //f.f      = &rootfinder_f;

  const gsl_multiroot_fsolver_type *TYPE = gsl_multiroot_fsolver_hybrids;
  solver = gsl_multiroot_fsolver_alloc(TYPE, 4);
  gsl_multiroot_fsolver_set(solver, &f, x);

  int status;
  size_t iter = 0;
  double previous_solver_step[4];

  ////////////////////
  // Loop.
  do
  {
    for (int iPrev = 0; iPrev < 4; iPrev++)
      previous_solver_step[iPrev] = gsl_vector_get(solver->x, iPrev);

    std::cout << "iter = " << iter << "\n";
std::cout << gsl_vector_get(solver->x, 0) << std::endl;
std::cout << gsl_vector_get(solver->x, 1) << std::endl;
std::cout << gsl_vector_get(solver->x, 2) << std::endl;
std::cout << gsl_vector_get(solver->x, 3) << std::endl;

    ++iter;
    status = gsl_multiroot_fsolver_iterate(solver);

    if ( status )
    {
      if ( VERBOSE > 5 )
      {
        if ( status == GSL_EBADFUNC )
          std::cout << "Error: something went to +/-Inf or NaN!\n";
        else if ( status == GSL_ENOPROG )
          std::cout << "Error: not making enough progress!\n";
        else if ( status == GSL_ENOPROGJ )
          std::cout << "Error: not making enough progress in Jacobian!\n";
        else
          std::cout << "Check: " << iter << "   " << steps << "   "
                    << status << std::endl;
      }

      //break if the rootfinder gets stuck
      break;
    }

    //break if the rootfinder goes out of bounds
    if(gsl_vector_get(solver->x, 0) < minT)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (T < minT)!\n";
      status = -10;
      break;
    }
    else if(gsl_vector_get(solver->x, 0) > maxT)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (T > maxT)!\n";
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 1) < minMuB)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuB < minMuB)!\n";
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 1) > maxMuB)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuB > maxMuB)!\n";
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 2) < minMuQ)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuQ < minMuQ)!\n";
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 2) > maxMuQ)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuQ > maxMuQ)!\n";
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 3) < minMuS)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuS < minMuS)!\n";
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 3) > maxMuS)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuS > maxMuS)!\n";
      status = -10;
      break;
    }

    if ( VERBOSE > 8 )
      std::cout << "\t --> Status: " << status << "   "
           << iter << "   " << error << "\n\t             ("
           << gsl_vector_get(solver->x, 0)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->x, 1)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->x, 2)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->x, 3)*hbarc_MeVfm << ")"
           << "\n\t             ("
           << gsl_vector_get(solver->dx, 0)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->dx, 1)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->dx, 2)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->dx, 3)*hbarc_MeVfm << ")"
           << "\n\t             ("
           << e_or_s_Given*hbarc_MeVfm << ","
           << rhoBGiven << ","
           << rhoQGiven << ","
           << rhoSGiven << ")"
           << "\n\t             ("
           << gsl_vector_get(solver->f, 0)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->f, 1) << ","
           << gsl_vector_get(solver->f, 2) << ","
           << gsl_vector_get(solver->f, 3) << ")" << std::endl;


    // test absolute error
    status = gsl_multiroot_test_residual(solver->f, error);

    // test relative error
    //status = gsl_multiroot_test_delta(solver->f, chosen_densities, 1e-15, 1e-15);

    // test convergence of root location
    //status = gsl_multiroot_test_delta(solver->dx, solver->x, 0.0, error);

//    // test all conditions at once
//    if (   gsl_multiroot_test_residual(solver->f, error) == GSL_CONTINUE
//        || gsl_multiroot_test_delta(solver->dx, solver->x, error, error) == GSL_CONTINUE
//        || gsl_multiroot_test_delta(solver->f, chosen_densities, error, error) == GSL_CONTINUE
//       )
//      status = GSL_CONTINUE;

  } while (status == GSL_CONTINUE && iter < steps);

  // check if a solution was found
  bool found = true;
  if ( iter >= steps || status != 0 )
  {
    if ( VERBOSE > 2 )
    {
      if ( status == GSL_EBADFUNC )
        std::cout << "Error: something went to +/-Inf or NaN!\n";
      else if ( status == GSL_ENOPROG )
        std::cout << "Error: not making enough progress!\n";
      else if ( status == GSL_ENOPROGJ )
        std::cout << "Error: not making enough progress in Jacobian!\n";
      else
        std::cout << "Check: " << iter << "   " << steps << "   "
                  << status << "\n";
    }
    found = false;
  }

//  std::cout << "Check: " << iter << "   " << steps << "   " << status << "\n";


  // if so, return the solution
  if ( found )
  {
    //DON'T NEED TO UPDATE TBQSPOSITION
    /*tbqs( gsl_vector_get(solver->x, 0),
          gsl_vector_get(solver->x, 1),
          gsl_vector_get(solver->x, 2),
          gsl_vector_get(solver->x, 3) );*/

    //UPDATE RESULT PASSED BACK TO EOS INSTEAD
    for (int iTBQS = 0; iTBQS < 4; iTBQS++)
      updated_tbqs[iTBQS] = gsl_vector_get(solver->x, iTBQS);
  }

  // memory deallocation
  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(x);

  return found;
}



////////////////////////////////////////////////////////////////////////////////
bool Rootfinder::find_root( const string & e_or_s, double ein_or_sin,
                            double Bin, double Sin, double Qin,
                            std::function<void(double[], double[])>
                                 function_to_evaluate,
                            vector<double> & tbqs_minima,
                            vector<double> & tbqs_maxima,
                            vector<double> & updated_tbqs )
{
    // set mode
    int e_or_s_mode = -1;
    if ( e_or_s == "e" || e_or_s == "energy" )
      e_or_s_mode = 1;
    else if ( e_or_s == "s" || e_or_s == "entropy" )
      e_or_s_mode = 0;
    else
    {
      std::cerr << "e_or_s = " << e_or_s << " is an invalid choice!\n";
      exit(8);
    }

    int number_of_attempts = 1;

    minT   = tbqs_minima[0]; maxT   = tbqs_maxima[0];
    minMuB = tbqs_minima[1]; maxMuB = tbqs_maxima[1];
    minMuQ = tbqs_minima[2]; maxMuQ = tbqs_maxima[2];
    minMuS = tbqs_minima[3]; maxMuS = tbqs_maxima[3];

    if ( VERBOSE > 6 )
      std::cout << "Using grid ranges: "
        << minT*hbarc_MeVfm   << "   " << maxT*hbarc_MeVfm << "   "
        << minMuB*hbarc_MeVfm << "   " << maxMuB*hbarc_MeVfm << "   "
        << minMuQ*hbarc_MeVfm << "   " << maxMuQ*hbarc_MeVfm << "   "
        << minMuS*hbarc_MeVfm << "   " << maxMuS*hbarc_MeVfm
        << "\n" << "TBQS seed values: "
        << updated_tbqs[0]*hbarc_MeVfm << "   "
        << updated_tbqs[1]*hbarc_MeVfm << "   "
        << updated_tbqs[2]*hbarc_MeVfm << "   "
        << updated_tbqs[3]*hbarc_MeVfm << std::endl;

    tbqsPosition = updated_tbqs;


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";


    //==========================================================================
    // try default seed point
    if (rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        return true;

    ///////////////////////////

    double t0 = tbqsPosition[0];
    double mub0 = tbqsPosition[1];
    double muq0 = tbqsPosition[2];
    double mus0 = tbqsPosition[3];
    double t10 = t0*.2;
    double muB10 = mub0*.2;
    double muQ10 = muq0*.2;
    double muS10 = mus0*.2;

    //==========================================================================
    // perturb T

    // perturb up
    number_of_attempts++;
    if (t0 + t10 > maxT)
      tbqs(maxT - 1, mub0, muq0, mus0);
    else
      tbqs(t0 + t10, mub0, muq0, mus0);


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check1: " << t0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << "\n";
    std::cout << "Check2: " << t10 << "   " << muB10 << "   " << muQ10 << "   " << muS10 << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";

    if ( rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                       function_to_evaluate, updated_tbqs ) ) 
        return true;

    // perturb down
    number_of_attempts++;
    if (t0 - t10 < minT)
      tbqs(minT + 1, mub0, muq0, mus0);
    else
      tbqs(t0 - t10, mub0, muq0, mus0);


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check1: " << t0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << "\n";
    std::cout << "Check2: " << t10 << "   " << muB10 << "   " << muQ10 << "   " << muS10 << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";

    if ( rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                       function_to_evaluate, updated_tbqs ) ) 
        return true;

    //==========================================================================
    // perturb muB

    // perturb up
    number_of_attempts++;
    if (mub0 + muB10 > maxMuB)
      tbqs(t0, maxMuB - 1, muq0, mus0);
    else
      tbqs(t0, mub0 + muB10, muq0, mus0);


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check1: " << t0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << "\n";
    std::cout << "Check2: " << t10 << "   " << muB10 << "   " << muQ10 << "   " << muS10 << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";

    if( rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs ) ) 
        return true;

    // perturb down
    number_of_attempts++;
    if (mub0 - muB10 < minMuB)
      tbqs(t0, minMuB + 1, muq0, mus0);
    else
      tbqs(t0, mub0 - muB10, muq0, mus0);


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check1: " << t0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << "\n";
    std::cout << "Check2: " << t10 << "   " << muB10 << "   " << muQ10 << "   " << muS10 << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";

    if ( rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs ) ) 
        return true;

    //==========================================================================
    // perturb muS

    // perturb up
    number_of_attempts++;
    if(mus0 + muS10 > maxMuS)
      tbqs(t0, mub0, muq0, maxMuS - 1);
    else
      tbqs(t0, mub0, muq0, mus0 + muS10);


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check1: " << t0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << "\n";
    std::cout << "Check2: " << t10 << "   " << muB10 << "   " << muQ10 << "   " << muS10 << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";

    if( rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs ) ) 
        return true;

    // perturb down
    number_of_attempts++;
    if(mus0 - muS10 < minMuS)
      tbqs(t0, mub0, muq0, minMuS + 1);
    else
      tbqs(t0, mub0, muq0, mus0 - muS10);


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check1: " << t0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << "\n";
    std::cout << "Check2: " << t10 << "   " << muB10 << "   " << muQ10 << "   " << muS10 << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";

    if( rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs ) ) 
        return true;

    //==========================================================================
    // perturb muQ

    // perturb up
    number_of_attempts++;
    if (muq0 + muQ10 > maxMuQ)
      tbqs(t0, mub0, maxMuQ - 1, mus0);
    else
      tbqs(t0, mub0, muq0 + muQ10, mus0);


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check1: " << t0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << "\n";
    std::cout << "Check2: " << t10 << "   " << muB10 << "   " << muQ10 << "   " << muS10 << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";

    if( rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs ) ) 
        return true;

    // perturb down
    number_of_attempts++;
    if(muq0 - muQ10 < minMuQ)
      tbqs(t0, mub0, minMuQ + 1, mus0);
    else
      tbqs(t0, mub0, muq0 - muQ10, mus0);


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check1: " << t0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << "\n";
    std::cout << "Check2: " << t10 << "   " << muB10 << "   " << muQ10 << "   " << muS10 << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";

    if( rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs ) ) 
        return true;

    //==========================================================================
    //check mu = 0
    tbqs(t0, 0, 0, 0);

    number_of_attempts++;


    std::cout << "Start of attempt #" << number_of_attempts << " at line " << __LINE__ << "\n";
    std::cout << "Check1: " << t0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << "\n";
    std::cout << "Check2: " << t10 << "   " << muB10 << "   " << muQ10 << "   " << muS10 << "\n";
    std::cout << "Check3: " << tbqsPosition[0] << "   " << tbqsPosition[1]
              << "   " << tbqsPosition[2] << "   " << tbqsPosition[3] << "\n";

    if( rootfinder4D( ein_or_sin, e_or_s_mode, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs ) ) 
        return true;

    if ( VERBOSE > 8 )
      std::cout << __PRETTY_FUNCTION__ << ": failed after "
                << number_of_attempts << "!\n";

    tbqs(t0, mub0, muq0, mus0);
    return false;
}

