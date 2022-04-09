#include <functional>

#include "rootfinder.h"

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
      double densities_at_point[4];

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
cout << "S: " << rhoS << "   " << rhoSGiven << "   " << rhoS - rhoSGiven << endl;
if (1) exit(1);


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
    std::cout << "T = " << setT << " is out of range. Valid values are between ["
      << minT << "," << maxT << "]" << std::endl;
    return;
  }
  if(setmuB < minMuB || setmuB > maxMuB) {
    std::cout << "muB = " << setmuB << " is out of range. Valid values are between ["
      << minMuB << "," << maxMuB << "]" << std::endl;
    return;
  }
  if(setmuQ < minMuQ || setmuQ > maxMuQ) {
    std::cout << "muQ = " << setmuQ << " is out of range. Valid values are between ["
      << minMuQ << "," << maxMuQ << "]" << std::endl;
    return;
  }
  if(setmuS < minMuS || setmuS > maxMuS) {
    std::cout << "muS = " << setmuS << " is out of range. Valid values are between ["
      << minMuS << "," << maxMuS << "]" << std::endl;
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


  ////////////////////
  // set initial guess
//cout << "Check initial guess:";
//for (int iTBQS = 0; iTBQS < 4; iTBQS++) cout << "   " << tbqsPosition[iTBQS];
//cout << endl;

  gsl_vector *x = gsl_vector_alloc(4);
  for (int iTBQS = 0; iTBQS < 4; iTBQS++)
    gsl_vector_set(x, iTBQS, tbqsPosition[iTBQS]);


  ////////////////////
  // pass relevant parameters to rootfinder
  rootfinder_parameters p( e_or_s_Given, rhoBGiven, rhoQGiven, rhoSGiven,
                           e_or_s_mode, function_to_evaluate );

  ////////////////////
  // initialize multiroot solver
  gsl_multiroot_fsolver *solver;
  gsl_multiroot_function f;

  f.n      = 4;
  f.params = &p;
  f.f      = &rootfinder_f;

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

    ++iter;
//cout << "Status before(1): " << status << "   " << GSL_CONTINUE << endl;
    status = gsl_multiroot_fsolver_iterate(solver);
//cout << "Status after(1): " << status << "   " << GSL_CONTINUE << endl;

    if(VERBOSE > 5 && status)
    {
      if ( status == GSL_EBADFUNC && VERBOSE > 5 )
        std::cout << "Error: something went to +/-Inf or NaN!" << std::endl;
      else if ( status == GSL_ENOPROG && VERBOSE > 5 )
        std::cout << "Error: not making enough progress!" << std::endl;
      else if ( status == GSL_ENOPROGJ && VERBOSE > 5 )
        std::cout << "Error: not making enough progress in Jacobian!" << std::endl;
      //break if the rootfinder gets stuck
      break;
    }

    //break if the rootfinder goes out of bounds
    if(gsl_vector_get(solver->x, 0) < minT)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (T < minT)!" << std::endl;
      status = -10;
      break;
    }
    else if(gsl_vector_get(solver->x, 0) > maxT)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (T > maxT)!" << std::endl;
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 1) < minMuB)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuB < minMuB)!" << std::endl;
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 1) > maxMuB)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuB > maxMuB)!" << std::endl;
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 2) < minMuQ)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuQ < minMuQ)!" << std::endl;
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 2) > maxMuQ)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuQ > maxMuQ)!" << std::endl;
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 3) < minMuS)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuS < minMuS)!" << std::endl;
      status = -10;
      break;
    }
    else if (gsl_vector_get(solver->x, 3) > maxMuS)
    {
      if ( VERBOSE > 5 )
        std::cout << "Error: out-of-bounds (MuS > maxMuS)!" << std::endl;
      status = -10;
      break;
    }

    if ( VERBOSE > 8 )
      std::cout << "\t --> Status: " << status << "   "
           << iter << "   " << error
           << std::endl << "\t             ("
           << gsl_vector_get(solver->x, 0)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->x, 1)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->x, 2)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->x, 3)*hbarc_MeVfm << ")"
           << std::endl << "\t             ("
           << e_or_s_Given*hbarc_MeVfm << ","
           << rhoBGiven << ","
           << rhoQGiven << ","
           << rhoSGiven << ")"
           << std::endl << "\t             ("
           << gsl_vector_get(solver->f, 0)*hbarc_MeVfm << ","
           << gsl_vector_get(solver->f, 1) << ","
           << gsl_vector_get(solver->f, 2) << ","
           << gsl_vector_get(solver->f, 3) << ")" << std::endl;

//cout << "Status before(2): " << status << "   " << GSL_CONTINUE << endl;
    status = gsl_multiroot_test_residual(solver->f, error);
//cout << "Status after(2): " << status << "   " << GSL_CONTINUE << endl;

  } while (status == GSL_CONTINUE && iter < steps);

  // check if a solution was found
  bool found = true;
  if ( iter >= steps || status != 0 )
  {
    if ( status == GSL_EBADFUNC )
      std::cout << "Error: something went to +/-Inf or NaN!" << std::endl;
    else if ( status == GSL_ENOPROG )
      std::cout << "Error: not making enough progress!" << std::endl;
    else if ( status == GSL_ENOPROGJ )
      std::cout << "Error: not making enough progress in Jacobian!" << std::endl;
    else
      std::cout << "Check: " << iter << "   " << steps << "   " << status << std::endl;
    found = false;
    //exit(8);
  }

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

/*cout << "\t --> " << __LINE__ << ": iter = " << iter << endl;
cout << "\t --> " << __LINE__ << ": "
      << gsl_vector_get(solver->x, 0) << "   "
      << gsl_vector_get(solver->x, 1) << "   "
      << gsl_vector_get(solver->x, 2) << "   "
      << gsl_vector_get(solver->x, 3) << endl;*/

  // memory deallocation
  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(x);

  return found;
}






////////////////////////////////////////////////////////////////////////////////
bool Rootfinder::find_eBSQ_root( double ein, double Bin, double Sin, double Qin,
                          std::function<void(double[], double[])> function_to_evaluate,
                           vector<double> & tbqs_minima,
                           vector<double> & tbqs_maxima,
                           vector<double> & updated_tbqs )
{
    int number_of_attempts = 1;
    minT   = tbqs_minima[0]; maxT   = tbqs_maxima[0];
    minMuB = tbqs_minima[1]; maxMuB = tbqs_maxima[1];
    minMuQ = tbqs_minima[2]; maxMuQ = tbqs_maxima[2];
    minMuS = tbqs_minima[3]; maxMuS = tbqs_maxima[3];
std::cout << "Using grid ranges: "
  << minT*hbarc_MeVfm << "   "
  << maxT*hbarc_MeVfm << "   "
  << minMuB*hbarc_MeVfm << "   "
  << maxMuB*hbarc_MeVfm << "   "
  << minMuQ*hbarc_MeVfm << "   "
  << maxMuQ*hbarc_MeVfm << "   "
  << minMuS*hbarc_MeVfm << "   "
  << maxMuS*hbarc_MeVfm << std::endl;
    tbqsPosition = updated_tbqs;

    if (rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { /*cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl;*/ return true; }

    ///////////////////////////

    double t0 = tbqsPosition[0];
    double mub0 = tbqsPosition[1];
    double muq0 = tbqsPosition[2];
    double mus0 = tbqsPosition[3];
    double t10 = t0*.2;
    double muB10 = mub0*.2;
    double muQ10 = muq0*.2;
    double muS10 = mus0*.2;

    //perturb T
    number_of_attempts++;

    if(t0 + t10 > maxT) {
        tbqs(maxT - 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 + t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs)) 
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    number_of_attempts++;

    if(t0 - t10 < minT) {
        tbqs(minT + 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 - t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs)) 
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    //perturb mub

    number_of_attempts++;

    if(mub0 + muB10 > maxMuB) {
        tbqs(t0, maxMuB - 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 + muB10, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs)) 
        { /*cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl;*/ return true; }

    number_of_attempts++;

    if(mub0 - muB10 < minMuB) {
        tbqs(t0, minMuB + 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 - muB10, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs)) 
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    //perturn muq

    number_of_attempts++;

    if(muq0 + muQ10 > maxMuQ) {
        tbqs(t0, mub0, maxMuQ - 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 + muQ10, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs)) 
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    number_of_attempts++;

    if(muq0 - muQ10 < minMuQ) {
        tbqs(t0, mub0, minMuQ + 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 - muQ10, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    //perturb mus

    number_of_attempts++;

    if(mus0 + muS10 > maxMuS) {
        tbqs(t0, mub0, muq0, maxMuS - 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 + muS10);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs)) 
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    number_of_attempts++;

    if(mus0 - muS10 < maxMuS) {
        tbqs(t0, mub0, muq0, minMuS + 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 - muS10);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs)) 
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    //check mu = 0
    tbqs(t0, 0, 0, 0);

    number_of_attempts++;

    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs)) 
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

cout << __PRETTY_FUNCTION__ << ": failed after " << number_of_attempts << "!" << endl;

    tbqs(t0, mub0, muq0, mus0);
    return false;
}






////////////////////////////////////////////////////////////////////////////////
bool Rootfinder::find_sBSQ_root( double sin, double Bin, double Sin, double Qin,
                           std::function<void(double[], double[])> function_to_evaluate,
                           vector<double> & tbqs_minima,
                           vector<double> & tbqs_maxima,
                           vector<double> & updated_tbqs )
{
    int number_of_attempts = 1;
    minT   = tbqs_minima[0]; maxT   = tbqs_maxima[0];
    minMuB = tbqs_minima[1]; maxMuB = tbqs_maxima[1];
    minMuQ = tbqs_minima[2]; maxMuQ = tbqs_maxima[2];
    minMuS = tbqs_minima[3]; maxMuS = tbqs_maxima[3];
    tbqsPosition = updated_tbqs;

    if (rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { /*cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl;*/ return true; }

	///////////////////////////
    double t0 = tbqsPosition[0];
    double mub0 = tbqsPosition[1];
    double muq0 = tbqsPosition[2];
    double mus0 = tbqsPosition[3];
    double t10 = t0*.2;
    double muB10 = mub0*.2;
    double muQ10 = muq0*.2;
    double muS10 = mus0*.2;

    //perturb T

    number_of_attempts++;

    if(t0 + t10 > maxT) {
        tbqs(maxT - 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 + t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    number_of_attempts++;

    if(t0 - t10 < minT) {
        tbqs(minT + 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 - t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    //perturb mub

    number_of_attempts++;

    if(mub0 + muB10 > maxMuB) {
        tbqs(t0, maxMuB - 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 + muB10, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    number_of_attempts++;

    if(mub0 - muB10 < minMuB) {
        tbqs(t0, minMuB + 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 - muB10, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    //perturn muq

    number_of_attempts++;

    if(muq0 + muQ10 > maxMuQ) {
        tbqs(t0, mub0, maxMuQ - 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 + muQ10, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    number_of_attempts++;

    if(muq0 - muQ10 < minMuQ) {
        tbqs(t0, mub0, minMuQ + 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 - muQ10, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    //perturb mus

    number_of_attempts++;

    if(mus0 + muS10 > maxMuS) {
        tbqs(t0, mub0, muq0, maxMuS - 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 + muS10);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    number_of_attempts++;

    if(mus0 - muS10 < maxMuS) {
        tbqs(t0, mub0, muq0, minMuS + 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 - muS10);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

    //check mu = 0
    tbqs(t0, 0, 0, 0);

    number_of_attempts++;

    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS,
                      function_to_evaluate, updated_tbqs))
        { cout << __PRETTY_FUNCTION__ << ": Completed in " << number_of_attempts
              << " attempts" << endl; return true; }

cout << __PRETTY_FUNCTION__ << ": failed after " << number_of_attempts << "!" << endl;

    tbqs(t0, mub0, muq0, mus0);
    return false;
}
