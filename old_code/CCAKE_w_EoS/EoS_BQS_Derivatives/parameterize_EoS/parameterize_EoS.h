#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>	    // gsl random number generators
#include <gsl/gsl_randist.h>	    // gsl random number distributions
#include <gsl/gsl_vector.h>	    // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>	    // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_linalg.h>

using namespace std;

struct density_data
{
	bool e_or_charge;
	size_t data_length;
	vector<double> muB;
	vector<double> muS;
	vector<double> muQ;
	vector<double> y;
	vector<double> sigma;
};

inline double get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr)
{
	return gsl_vector_get (solver_ptr->x, i);
}

inline double get_fit_err (int i, gsl_matrix * covariance_ptr)
{
	return sqrt (gsl_matrix_get (covariance_ptr, i, i));
}


//*********************************************************************
//  Function returning the residuals for each point; that is, the 
//  difference of the fit function using the current parameters
//  and the data to be fit.
int Fittarget_density_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
	bool e_or_charge = ((struct density_data *) params_ptr)->e_or_charge;
	size_t n = ((struct density_data *) params_ptr)->data_length;
	vector<double> & muB = ((struct density_data *) params_ptr)->muB;
	vector<double> & muS = ((struct density_data *) params_ptr)->muS;
	vector<double> & muQ = ((struct density_data *) params_ptr)->muQ;
	vector<double> sigma = ((struct density_data *) params_ptr)->sigma;
	vector<double> y = ((struct density_data *) params_ptr)->y;

	//fit parameters
	double d0 = gsl_vector_get (xvec_ptr, 0);
	double a = gsl_vector_get (xvec_ptr, 1);
	double b = gsl_vector_get (xvec_ptr, 2);
	double c = gsl_vector_get (xvec_ptr, 3);

	size_t i;

	if ( e_or_charge )
		for (i = 0; i < n; i++)
		{
			double Yi = d0*cosh(a*muB[i]+b*muS[i]+c*muQ[i]);
			gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
		}
	else
		for (i = 0; i < n; i++)
		{
			double Yi = d0*sinh(a*muB[i]+b*muS[i]+c*muQ[i]);
			gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
		}

	return GSL_SUCCESS;
}

//*********************************************************************
//  Function returning the Jacobian of the residual function
int Fittarget_density_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
	bool e_or_charge = ((struct density_data *) params_ptr)->e_or_charge;
	size_t n = ((struct density_data *) params_ptr)->data_length;
	vector<double> & muB = ((struct density_data *) params_ptr)->muB;
	vector<double> & muS = ((struct density_data *) params_ptr)->muS;
	vector<double> & muQ = ((struct density_data *) params_ptr)->muQ;
	vector<double> sigma = ((struct density_data *) params_ptr)->sigma;
	vector<double> y = ((struct density_data *) params_ptr)->y;

	//fit parameters
	double d0 = gsl_vector_get (xvec_ptr, 0);
	double a = gsl_vector_get (xvec_ptr, 1);
	double b = gsl_vector_get (xvec_ptr, 2);
	double c = gsl_vector_get (xvec_ptr, 3);

	size_t i;

	if ( e_or_charge )
		for (i = 0; i < n; i++)
		{
			double sig = sigma[i];
	
			// derivatives
			double sinh_arg = sinh(a*muB[i]+b*muS[i]+c*muQ[i]);
			double cosh_arg = cosh(a*muB[i]+b*muS[i]+c*muQ[i]);
	      
			gsl_matrix_set (Jacobian_ptr, i, 0, cosh_arg/sig);
			gsl_matrix_set (Jacobian_ptr, i, 1, d0*muB[i]*sinh_arg/sig);
			gsl_matrix_set (Jacobian_ptr, i, 2, d0*muS[i]*sinh_arg/sig);
			gsl_matrix_set (Jacobian_ptr, i, 3, d0*muQ[i]*sinh_arg/sig);
		}
	else
		for (i = 0; i < n; i++)
		{
			double sig = sigma[i];
	
			// derivatives
			double sinh_arg = sinh(a*muB[i]+b*muS[i]+c*muQ[i]);
			double cosh_arg = cosh(a*muB[i]+b*muS[i]+c*muQ[i]);
	      
			gsl_matrix_set (Jacobian_ptr, i, 0, sinh_arg/sig);
			gsl_matrix_set (Jacobian_ptr, i, 1, d0*muB[i]*cosh_arg/sig);
			gsl_matrix_set (Jacobian_ptr, i, 2, d0*muS[i]*cosh_arg/sig);
			gsl_matrix_set (Jacobian_ptr, i, 3, d0*muQ[i]*cosh_arg/sig);
		}


	return GSL_SUCCESS;
}


//*********************************************************************
//  Function combining the residual function and its Jacobian
int Fittarget_density_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
	Fittarget_density_f(xvec_ptr, params_ptr, f_ptr);
	Fittarget_density_df(xvec_ptr, params_ptr, Jacobian_ptr);

	return GSL_SUCCESS;
}

void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
	return;
}

void fit( vector<double> & muBvec, vector<double> & muSvec, vector<double> & muQvec,
		  vector<double> & fvec, const size_t nmuB, const size_t nmuS, const size_t nmuQ,
		  bool e_or_charge )
{
	const int VERBOSE = 10;
	const int fit_max_iterations = 1000;
	const double fit_tolerance = 1e-6;
	const size_t data_length = nmuB*nmuS*nmuQ;  // # of points
	const size_t n_para = 4;  // # of parameters

	// allocate space for a covariance matrix of size p by p
	gsl_matrix *covariance_ptr = gsl_matrix_alloc (n_para, n_para);


	//set up test data
	struct density_data f_data;
	f_data.e_or_charge = e_or_charge;
	f_data.data_length = data_length;
	f_data.muB.resize(data_length);
	f_data.muS.resize(data_length);
	f_data.muQ.resize(data_length);
	f_data.y.resize(data_length);
	f_data.sigma.resize(data_length);

	cout << "Attempting to fit the following data: " << endl;

	int idx = 0;
	for (int i = 0; i < nmuB; i++)
	for (int j = 0; j < nmuS; j++)
	for (int k = 0; k < nmuQ; k++)
	{
		f_data.muB[idx] = muBvec[idx];
		f_data.muS[idx] = muSvec[idx];
		f_data.muQ[idx] = muQvec[idx];
		f_data.y[idx] = fvec[idx];
		f_data.sigma[idx] = 1e-1;
		cout << "   " << idx
			<< "   " << muBvec[idx]
			<< "   " << muSvec[idx]
			<< "   " << muQvec[idx]
			<< "   " << fvec[idx] << endl;
		idx++;
	}

	//double para_init[n_para] = { 1.0, 1.0, 1.0, 1.0 };  // initial guesses of parameters

	//gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);


const size_t n = data_length;
const size_t p = n_para;


gsl_vector *f;
  gsl_matrix *J;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double t[data_length], y[data_length], weights[data_length];
  //struct data d = { n, t, y };
  double x_init[4] = { 1.0, 0.01, 0.01, 0.01 }; /* starting values */
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  gsl_vector_view wts = gsl_vector_view_array(weights, n);
  double chisq, chisq0;
  int status, info;
  size_t i;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;






  
	// set up the function to be fit 
	/*gsl_multifit_function_fdf target_func;
	target_func.f = &Fittarget_density_f;        // the function of residuals
	target_func.df = &Fittarget_density_df;      // the gradient of this function
	target_func.fdf = &Fittarget_density_fdf;    // combined function and gradient
	target_func.n = data_length;              // number of points in the data set
	target_func.p = n_para;              // number of parameters in the fit function
	target_func.params = &f_data;  // structure with the data and error bars

	const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver_ptr = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
	gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);*/

	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_workspace *w;
	gsl_multifit_nlinear_fdf fdf;
	gsl_multifit_nlinear_parameters fdf_params =
	gsl_multifit_nlinear_default_parameters();


/* define the function to be minimized */
  fdf.f = Fittarget_density_f;
  fdf.df = Fittarget_density_df;   /* set to NULL for finite-difference Jacobian */
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = data_length;
  fdf.p = n_para;
  fdf.params = &f_data;


	/*
	size_t iteration = 0;         // initialize iteration counter
	//if (VERBOSE > 2) print_fit_state_3D (iteration, solver_ptr);
	int status;  		// return value from gsl function calls (e.g., error)
	do
	{
		iteration++;
      
		// perform a single iteration of the fitting routine
		status = gsl_multifit_fdfsolver_iterate (solver_ptr);

		// print out the status of the fit
		if (VERBOSE > 2) cout << "status = " << gsl_strerror (status) << endl;

		// customized routine to print out current parameters
		//if (VERBOSE > 2) print_fit_state_3D (iteration, solver_ptr);

		if (status)    // check for a nonzero status code
		{
			break;  // this should only happen if an error code is returned 
		}

		// test for convergence with an absolute and relative error (see manual)
		status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, fit_tolerance, fit_tolerance);
	}
	while (status == GSL_CONTINUE && iteration < fit_max_iterations);

	// calculate the covariance matrix of the best-fit parameters
	gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

	// print out the covariance matrix using the gsl function (not elegant!)
	if (VERBOSE > 2) cout << endl << "Covariance matrix: " << endl;
	if (VERBOSE > 2) gsl_matrix_fprintf (stdout, covariance_ptr, "%g");*/



/* allocate workspace with default parameters */
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

  /* compute initial cost function */
  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f, f, &chisq0);

  /* solve the system with a maximum of 100 iterations */
  status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol,
                                       callback, NULL, &info, w);

  /* compute covariance of best fit parameters */
  J = gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar (J, 0.0, covar);

  /* compute final cost */
  gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  fprintf(stderr, "summary from method '%s/%s'\n",
          gsl_multifit_nlinear_name(w),
          gsl_multifit_nlinear_trs_name(w));
  fprintf(stderr, "number of iterations: %zu\n",
          gsl_multifit_nlinear_niter(w));
  fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
  fprintf(stderr, "reason for stopping: %s\n",
          (info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
  fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));


{
    double dof = n - p;
    double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

    fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

    fprintf (stderr, "d0      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    fprintf (stderr, "a = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    fprintf (stderr, "b      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
    fprintf (stderr, "c      = %.5f +/- %.5f\n", FIT(3), c*ERR(3));
  }

  fprintf (stderr, "status = %s\n", gsl_strerror (status));

  gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);


/*
	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		                // # of digits in doubles

	double chi = gsl_blas_dnrm2(solver_ptr->f);
	double dof = data_length - n_para;
	double c = GSL_MAX_DBL(1, chi/sqrt(dof));

	cout << "Solution:" << endl;
	cout << "\t" << get_fit_results(0, solver_ptr) << " +/- " << get_fit_err(0, covariance_ptr)
			<< "   " << get_fit_results(1, solver_ptr) << " +/- " << get_fit_err(1, covariance_ptr)
			<< "   " << get_fit_results(2, solver_ptr) << " +/- " << get_fit_err(2, covariance_ptr)
			<< "   " << get_fit_results(3, solver_ptr) << " +/- " << get_fit_err(3, covariance_ptr)
			<< endl;

	//clean up
	gsl_matrix_free (covariance_ptr);
	gsl_rng_free (rng_ptr);

	gsl_multifit_fdfsolver_free (solver_ptr);  // free up the solver
*/
	
	return;
}
