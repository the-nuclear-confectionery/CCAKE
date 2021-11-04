#include "rootfinder.h"




////////////////////////////////////////////////////////////////////////////////
//struct to pass the target (E, rhoB, rhoQ, rhoS) into the rootfinder function
struct rootfinder_parameters
  {
      double eorEntGiven;          //these are the desired s and BSQ
      double rhoBGiven;
      double rhoQGiven;
      double rhoSGiven;
      rootfinder_parameters();
      rootfinder_parameters(
      double seteorEntGiven, double setRhoBGiven,
      double setRhoQGiven, double setRhoSGiven);
  public:
      void set( double setEorEntGiven, double setRhoBGiven,
          double setRhoQGiven, double setRhoSGiven);
  };

  rootfinder_parameters::rootfinder_parameters() {}
  rootfinder_parameters::rootfinder_parameters(
    double setEorEntGiven, double setRhoBGiven,
    double setRhoQGiven, double setRhoSGiven
    )
  {
      eorEntGiven = setEorEntGiven;
      rhoBGiven = setRhoBGiven;
      rhoQGiven = setRhoQGiven;
      rhoSGiven = setRhoSGiven;
  }

  void rootfinder_parameters::set(
    double setEorEntGiven, double setRhoBGiven,
    double setRhoQGiven, double setRhoSGiven)
  {
      eorEntGiven = setEorEntGiven;
      rhoBGiven = setRhoBGiven;
      rhoQGiven = setRhoQGiven;
      rhoSGiven = setRhoSGiven;

  }



////////////////////////////////////////////////////////////////////////////////
int rootfinder_fsbqs(const gsl_vector *x, void *params, gsl_vector *f)
{
    //x contains the next (T, muB, muS) coordinate to test
    vector<double> tbqsToEval(4);
    tbqsToEval[0] = gsl_vector_get(x,0);
    tbqsToEval[1] = gsl_vector_get(x,1);	// convert x into densevector so it
    tbqsToEval[2] = gsl_vector_get(x,2);	// can be a BSpline evaluation point
    tbqsToEval[3] = gsl_vector_get(x,3);


    double entrGiven, rhoBGiven, rhoQGiven, rhoSGiven, entr, rhoB, rhoQ, rhoS;
    entrGiven = ((rootfinder_parameters*)params)->eorEntGiven;
    rhoBGiven = ((rootfinder_parameters*)params)->rhoBGiven;            //given variables contain the target point
    rhoQGiven = ((rootfinder_parameters*)params)->rhoQGiven;
    rhoSGiven = ((rootfinder_parameters*)params)->rhoSGiven;
{
	double phase_diagram_point[4] = {tbqsToEval[0]*hbarc_MeVfm, tbqsToEval[1]*hbarc_MeVfm,
					 tbqsToEval[3]*hbarc_MeVfm, tbqsToEval[2]*hbarc_MeVfm};	// NOTE: S <<-->> Q swapped!!!
	double densities_at_point[4];
	get_sBSQ_densities(phase_diagram_point, densities_at_point);
	entr = densities_at_point[0];
	rhoB = densities_at_point[1];
	rhoS = densities_at_point[2];
	rhoQ = densities_at_point[3];

}

    gsl_vector_set(f, 0, (entr - entrGiven)); //f[0] contains (s(T,muB,muQ,muS) - sGiven)
    gsl_vector_set(f, 1, (rhoB - rhoBGiven)); //f[1] contains (rhoB(T,muB,muQ,muS) - rhoBGiven)
    gsl_vector_set(f, 2, (rhoQ - rhoQGiven)); //f[2] contains (rhoQ(T,muB,muQ,muS) - rhoQGiven)
    gsl_vector_set(f, 3, (rhoS - rhoSGiven)); //f[2] contains (rhoS(T,muB,muQ,muS) - rhoSGiven)

    return GSL_SUCCESS;
}



////////////////////////////////////////////////////////////////////////////////
int rootfinder_febqs(const gsl_vector *x, void *params, gsl_vector *f)
{
    //x contains the next (T, muB, muQ, muS) coordinate to test
    vector<double> tbqsToEval(4);
    tbqsToEval[0] = gsl_vector_get(x,0);
    tbqsToEval[1] = gsl_vector_get(x,1);	// convert x into densevector so it
    tbqsToEval[2] = gsl_vector_get(x,2);	// can be a BSpline evaluation point
    tbqsToEval[3] = gsl_vector_get(x,3);

    double eGiven, rhoBGiven, rhoQGiven, rhoSGiven, e, rhoB, rhoQ, rhoS;
    eGiven = ((rootfinder_parameters*)params)->eorEntGiven;
    rhoBGiven = ((rootfinder_parameters*)params)->rhoBGiven;            //given variables contain the target point
    rhoQGiven = ((rootfinder_parameters*)params)->rhoQGiven;
    rhoSGiven = ((rootfinder_parameters*)params)->rhoSGiven;
{
	double phase_diagram_point[4] = {tbqsToEval[0]*hbarc_MeVfm, tbqsToEval[1]*hbarc_MeVfm,
					 tbqsToEval[3]*hbarc_MeVfm, tbqsToEval[2]*hbarc_MeVfm};	// NOTE: S <<-->> Q swapped!!!
	double densities_at_point[4];
	get_eBSQ_densities(phase_diagram_point, densities_at_point);
	e = densities_at_point[0]/hbarc_MeVfm;
	rhoB = densities_at_point[1];
	rhoS = densities_at_point[2];
	rhoQ = densities_at_point[3];
}

    gsl_vector_set(f, 0, (e - eGiven)); //f[0] contains (e(T,muB,muQ,muS) - eGiven)
    gsl_vector_set(f, 1, (rhoB - rhoBGiven)); //f[1] contains the (rhoB(T,muB,muQ,muS) - rhoBGiven)
    gsl_vector_set(f, 2, (rhoQ - rhoQGiven)); //f[2] contains the (rhoQ(T,muB,muQ,muS) - rhoQGiven)
    gsl_vector_set(f, 3, (rhoS - rhoSGiven)); //f[2] contains the (rho2(T,muB,muQ,muS) - rhoSGiven)

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
						double error, size_t steps)
{
  ////////////////////
  // set initial guess
  gsl_vector *x = gsl_vector_alloc(4);
  for (int iTBQS = 0; iTBQS < 4; iTBQS++)
    gsl_vector_set(x, iTBQS, tbqsPosition[iTBQS]);


  ////////////////////
  // decide if entropy or energy density passed in
  bool isEntropy = false;
  if ( e_or_s_mode == 0 ) isEntropy = true;

  ////////////////////
  // pass relevant parameters to rootfinder
  rootfinder_parameters p;
  p.set( e_or_s_Given, rhoBGiven, rhoQGiven, rhoSGiven);

  ////////////////////
  // initialize multiroot solver
  gsl_multiroot_fsolver *solver;
  gsl_multiroot_function f;

  if ( isEntropy )
    f.f = &rootfinder_fsbqs;
  else
    f.f = &rootfinder_febqs;

  f.n = 4;
  f.params = &p;

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
    status = gsl_multiroot_fsolver_iterate(solver);

    if(VERBOSE > 5 && status)
    {
      if ( status == GSL_EBADFUNC && VERBOSE > 5 )
        std::cout << "Error: something went to +/-Inf or NaN!" << std::endl;
      else if ( status == GSL_ENOPROG && VERBOSE > 5 )
        std::cout << "Error: not making enough progress!" << std::endl;
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

    status = gsl_multiroot_test_residual(solver->f, error);

  } while (status == GSL_CONTINUE && iter < steps);

  // check if a solution was found
  bool found = true;
  if ( iter >= steps || status != 0 )
  {
    if ( status == GSL_EBADFUNC )
      std::cout << "Error: something went to +/-Inf or NaN!" << std::endl;
    else if ( status == GSL_ENOPROG )
      std::cout << "Error: not making enough progress!" << std::endl;
    else
      std::cout << "Check: " << iter << "   " << steps << "   " << status << std::endl;
    found = false;
  }

  // if so, return the solution
  if ( found )
  {
    tbqs( gsl_vector_get(solver->x, 0),
          gsl_vector_get(solver->x, 1),
          gsl_vector_get(solver->x, 2),
          gsl_vector_get(solver->x, 3) );
  }

  // memory deallocation
  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(x);

  return found;
}






////////////////////////////////////////////////////////////////////////////////
bool Rootfinder::find_eBSQ_root( double ein, double Bin, double Sin, double Qin,
                          vector<double> & updated_tbqs )
{
    tbqsPosition = updated_tbqs;

    if (rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) { return true; }

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
    if(t0 + t10 > maxT) {
        tbqs(maxT - 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 + t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(t0 - t10 < minT) {
        tbqs(minT + 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 - t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //perturb mub
    if(mub0 + muB10 > maxMuB) {
        tbqs(t0, maxMuB - 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 + muB10, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(mub0 - muB10 < minMuB) {
        tbqs(t0, minMuB + 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 - muB10, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //perturn muq
    if(muq0 + muQ10 > maxMuQ) {
        tbqs(t0, mub0, maxMuQ - 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 + muQ10, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(muq0 - muQ10 < minMuQ) {
        tbqs(t0, mub0, minMuQ + 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 - muQ10, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //perturb mus
    if(mus0 + muS10 > maxMuS) {
        tbqs(t0, mub0, muq0, maxMuS - 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 + muS10);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(mus0 - muS10 < maxMuS) {


        tbqs(t0, mub0, muq0, minMuS + 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 - muS10);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //check mu = 0
    tbqs(t0, 0, 0, 0);
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    tbqs(t0, mub0, muq0, mus0);
    return false;
}






////////////////////////////////////////////////////////////////////////////////
bool Rootfinder::find_sBSQ_root( double sin, double Bin, double Sin, double Qin,
                           vector<double> & updated_tbqs )
{
    tbqsPosition = updated_tbqs;

    if (rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) { return true; }

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
    if(t0 + t10 > maxT) {
        tbqs(maxT - 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 + t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(t0 - t10 < minT) {
        tbqs(minT + 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 - t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //perturb mub
    if(mub0 + muB10 > maxMuB) {
        tbqs(t0, maxMuB - 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 + muB10, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(mub0 - muB10 < minMuB) {
        tbqs(t0, minMuB + 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 - muB10, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //perturn muq
    if(muq0 + muQ10 > maxMuQ) {
        tbqs(t0, mub0, maxMuQ - 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 + muQ10, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(muq0 - muQ10 < minMuQ) {
        tbqs(t0, mub0, minMuQ + 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 - muQ10, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //perturb mus
    if(mus0 + muS10 > maxMuS) {
        tbqs(t0, mub0, muq0, maxMuS - 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 + muS10);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(mus0 - muS10 < maxMuS) {
        tbqs(t0, mub0, muq0, minMuS + 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 - muS10);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //check mu = 0
    tbqs(t0, 0, 0, 0);
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    tbqs(t0, mub0, muq0, mus0);
    return false;
}
