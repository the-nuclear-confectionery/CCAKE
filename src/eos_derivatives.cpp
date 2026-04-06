#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "../include/eos.h"
using namespace ccake;

double EquationOfState::dentr_dt()   { return calc_term_1();        }
double EquationOfState::dentr_dmub() { return calc_term_2(ChargeIdx::B); }
double EquationOfState::dentr_dmuq() { return calc_term_2(ChargeIdx::Q); }
double EquationOfState::dentr_dmus() { return calc_term_2(ChargeIdx::S); }
double EquationOfState::db_dt()      { return calc_term_3(ChargeIdx::B); }
double EquationOfState::db_dmub()    { return calc_term_4(ChargeIdx::B, ChargeIdx::B); }
double EquationOfState::db_dmuq()    { return calc_term_4(ChargeIdx::B, ChargeIdx::Q); }
double EquationOfState::db_dmus()    { return calc_term_4(ChargeIdx::B, ChargeIdx::S); }
double EquationOfState::ds_dt()      { return calc_term_3(ChargeIdx::S); }
double EquationOfState::ds_dmub()    { return calc_term_4(ChargeIdx::S, ChargeIdx::B); }
double EquationOfState::ds_dmuq()    { return calc_term_4(ChargeIdx::S, ChargeIdx::Q); }
double EquationOfState::ds_dmus()    { return calc_term_4(ChargeIdx::S, ChargeIdx::S); }
double EquationOfState::dq_dt()      { return calc_term_3(ChargeIdx::Q); }
double EquationOfState::dq_dmub()    { return calc_term_4(ChargeIdx::Q, ChargeIdx::B); }
double EquationOfState::dq_dmuq()    { return calc_term_4(ChargeIdx::Q, ChargeIdx::Q); }
double EquationOfState::dq_dmus()    { return calc_term_4(ChargeIdx::Q, ChargeIdx::S); }

double EquationOfState::calc_term_1() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    // Reduced modes:
    // - T-only: no conserved charges -> term_1 reduces to dt2
    if (T_only_mode) return dt2;

    // - baryon-only: only B contributes -> 1D Schur complement
    if (baryon_only_mode)
    {
      const double den = db2 + 1e-30; // protect
      return dt2 - (dtdb*dtdb) / den;
    }
    gsl_vector_set(_dv,0,dtdb);
    gsl_vector_set(_dv,1,dtds);
    gsl_vector_set(_dv,2,dtdq);

/*if (verbose)
{
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "calc_term_1 check: v = ";
	for (int iGSL = 0; iGSL < 3; iGSL++)
		cout << gsl_vector_set(v,iGSL);
	cout << endl;
}*/

    gsl_matrix_set(_dm,0,0,db2);
    gsl_matrix_set(_dm,0,1,dbds);
    gsl_matrix_set(_dm,0,2,dbdq);
    gsl_matrix_set(_dm,1,0,dbds);
    gsl_matrix_set(_dm,1,1,ds2);
    gsl_matrix_set(_dm,1,2,dsdq);
    gsl_matrix_set(_dm,2,0,dbdq);
    gsl_matrix_set(_dm,2,1,dsdq);
    gsl_matrix_set(_dm,2,2,dq2);

    double toReturn = dt2 - deriv_mult_aTm_1b(_dv,_dm,_dv);

    return toReturn;
}

double EquationOfState::calc_term_2(ChargeIdx i) {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__
								 << ": i = " << static_cast<int>(i) << std::endl;
    if (T_only_mode) return 0.0;

    // baryon-only: only B active, so compute with 1x1 reduction
    if (baryon_only_mode)
    {
      if (i == ChargeIdx::B)
      {
        const double den = dtdb + 1e-30;
        return dtdb - (dt2 * db2) / den;
      }
      // muS/muQ derivatives don't exist in baryon-only
      return 0.0;
    }
    double toReturn = 0;

    if (i == ChargeIdx::B) {
        gsl_vector_set(_da,0,dt2);
        gsl_vector_set(_da,1,dtds);
        gsl_vector_set(_da,2,dtdq);

        gsl_vector_set(_db,0,db2);
        gsl_vector_set(_db,1,dbds);
        gsl_vector_set(_db,2,dbdq);

        gsl_matrix_set(_dm,0,0,dtdb);
        gsl_matrix_set(_dm,0,1,dbds);
        gsl_matrix_set(_dm,0,2,dbdq);
        gsl_matrix_set(_dm,1,0,dtds);
        gsl_matrix_set(_dm,1,1,ds2);
        gsl_matrix_set(_dm,1,2,dsdq);
        gsl_matrix_set(_dm,2,0,dtdq);
        gsl_matrix_set(_dm,2,1,dsdq);
        gsl_matrix_set(_dm,2,2,dq2);

        toReturn = dtdb - deriv_mult_aTm_1b(_da,_dm,_db);
    } else if (i == ChargeIdx::S) {
        gsl_vector_set(_da,0,dt2);
        gsl_vector_set(_da,1,dtdb);
        gsl_vector_set(_da,2,dtdq);

        gsl_vector_set(_db,0,dbds);
        gsl_vector_set(_db,1,ds2);
        gsl_vector_set(_db,2,dsdq);

        gsl_matrix_set(_dm,0,0,dtdb);
        gsl_matrix_set(_dm,0,1,db2);
        gsl_matrix_set(_dm,0,2,dbdq);
        gsl_matrix_set(_dm,1,0,dtds);
        gsl_matrix_set(_dm,1,1,dbds);
        gsl_matrix_set(_dm,1,2,dsdq);
        gsl_matrix_set(_dm,2,0,dtdq);
        gsl_matrix_set(_dm,2,1,dbdq);
        gsl_matrix_set(_dm,2,2,dq2);

        toReturn = dtds - deriv_mult_aTm_1b(_da,_dm,_db);
    } else if (i == ChargeIdx::Q) {
        gsl_vector_set(_da,0,dt2);
        gsl_vector_set(_da,1,dtdb);
        gsl_vector_set(_da,2,dtds);

        gsl_vector_set(_db,0,dbdq);
        gsl_vector_set(_db,1,dsdq);
        gsl_vector_set(_db,2,dq2);

        gsl_matrix_set(_dm,0,0,dtdb);
        gsl_matrix_set(_dm,0,1,db2);
        gsl_matrix_set(_dm,0,2,dbds);
        gsl_matrix_set(_dm,1,0,dtds);
        gsl_matrix_set(_dm,1,1,dbds);
        gsl_matrix_set(_dm,1,2,ds2);
        gsl_matrix_set(_dm,2,0,dtdq);
        gsl_matrix_set(_dm,2,1,dbdq);
        gsl_matrix_set(_dm,2,2,dsdq);

        toReturn = dtdq - deriv_mult_aTm_1b(_da,_dm,_db);
    } else {
        std::cout << "Error calculating derivative term 2" << std::endl;
    }

    return toReturn;
}

double EquationOfState::calc_term_3(ChargeIdx i) {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__
								 << ": i = " << static_cast<int>(i) << std::endl;
    double toReturn = 0;
    if (T_only_mode) return 0.0;

    if (baryon_only_mode)
    {
      if (i == ChargeIdx::B)
      {
        const double den = dtdb + 1e-30;
        return dtdb - (db2 * dt2) / den;
      }
      return 0.0;
    }
    if (i == ChargeIdx::B) {
        gsl_vector_set(_da,0,db2);
        gsl_vector_set(_da,1,dbds);
        gsl_vector_set(_da,2,dbdq);

        gsl_vector_set(_db,0,dt2);
        gsl_vector_set(_db,1,dtds);
        gsl_vector_set(_db,2,dtdq);

        gsl_matrix_set(_dm,0,0,dtdb);
        gsl_matrix_set(_dm,0,1,dtds);
        gsl_matrix_set(_dm,0,2,dtdq);
        gsl_matrix_set(_dm,1,0,dbds);
        gsl_matrix_set(_dm,1,1,ds2);
        gsl_matrix_set(_dm,1,2,dsdq);
        gsl_matrix_set(_dm,2,0,dbdq);
        gsl_matrix_set(_dm,2,1,dsdq);
        gsl_matrix_set(_dm,2,2,dq2);

        toReturn = dtdb - deriv_mult_aTm_1b(_da,_dm,_db);
    } else if (i == ChargeIdx::S) {
        gsl_vector_set(_da,0,dbds);
        gsl_vector_set(_da,1,ds2);
        gsl_vector_set(_da,2,dsdq);

        gsl_vector_set(_db,0,dt2);
        gsl_vector_set(_db,1,dtdb);
        gsl_vector_set(_db,2,dtdq);

        gsl_matrix_set(_dm,0,0,dtdb);
        gsl_matrix_set(_dm,0,1,dtds);
        gsl_matrix_set(_dm,0,2,dtdq);
        gsl_matrix_set(_dm,1,0,db2);
        gsl_matrix_set(_dm,1,1,dbds);
        gsl_matrix_set(_dm,1,2,dbdq);
        gsl_matrix_set(_dm,2,0,dbdq);
        gsl_matrix_set(_dm,2,1,dsdq);
        gsl_matrix_set(_dm,2,2,dq2);

        toReturn = dtds - deriv_mult_aTm_1b(_da,_dm,_db);
    } else if (i == ChargeIdx::Q) {
        gsl_vector_set(_da,0,dbdq);
        gsl_vector_set(_da,1,dsdq);
        gsl_vector_set(_da,2,dq2);

        gsl_vector_set(_db,0,dt2);
        gsl_vector_set(_db,1,dtdb);
        gsl_vector_set(_db,2,dtds);

        gsl_matrix_set(_dm,0,0,dtdb);
        gsl_matrix_set(_dm,0,1,dtds);
        gsl_matrix_set(_dm,0,2,dtdq);
        gsl_matrix_set(_dm,1,0,db2);
        gsl_matrix_set(_dm,1,1,dbds);
        gsl_matrix_set(_dm,1,2,dbdq);
        gsl_matrix_set(_dm,2,0,dbds);
        gsl_matrix_set(_dm,2,1,ds2);
        gsl_matrix_set(_dm,2,2,dsdq);

        toReturn = dtdq - deriv_mult_aTm_1b(_da,_dm,_db);
    } else {
        std::cout << "Error calculating derivative term 3" << std::endl;
    }

    return toReturn;
}

double EquationOfState::calc_term_4(ChargeIdx j, ChargeIdx i) {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__
								 << ": j, i = "
								<< static_cast<int>(j) << "   " << static_cast<int>(i) << std::endl;
    if (T_only_mode) return 0.0;

    if (baryon_only_mode)
    {
      if (i == ChargeIdx::B && j == ChargeIdx::B)
      {
        const double den = dt2 + 1e-30;
        return db2 - (dtdb*dtdb) / den;
      }
      return 0.0;
    }
    double toReturn = 0;

    if (i == ChargeIdx::B) {
        if(j == ChargeIdx::B) {
            gsl_vector_set(_da,0,dtdb);
            gsl_vector_set(_da,1,dbds);
            gsl_vector_set(_da,2,dbdq);

            gsl_vector_set(_db,0,dtdb);
            gsl_vector_set(_db,1,dbds);
            gsl_vector_set(_db,2,dbdq);

            gsl_matrix_set(_dm,0,0,dt2);
            gsl_matrix_set(_dm,0,1,dtds);
            gsl_matrix_set(_dm,0,2,dtdq);
            gsl_matrix_set(_dm,1,0,dtds);
            gsl_matrix_set(_dm,1,1,ds2);
            gsl_matrix_set(_dm,1,2,dsdq);
            gsl_matrix_set(_dm,2,0,dtdq);
            gsl_matrix_set(_dm,2,1,dsdq);
            gsl_matrix_set(_dm,2,2,dq2);

            toReturn = db2 - deriv_mult_aTm_1b(_da,_dm,_db);
        } else if (j == ChargeIdx::S) {
            gsl_vector_set(_da,0,dtds);
            gsl_vector_set(_da,1,ds2);
            gsl_vector_set(_da,2,dsdq);

            gsl_vector_set(_db,0,dtdb);
            gsl_vector_set(_db,1,db2);
            gsl_vector_set(_db,2,dbdq);

            gsl_matrix_set(_dm,0,0,dt2);
            gsl_matrix_set(_dm,0,1,dtds);
            gsl_matrix_set(_dm,0,2,dtdq);
            gsl_matrix_set(_dm,1,0,dtdb);
            gsl_matrix_set(_dm,1,1,dbds);
            gsl_matrix_set(_dm,1,2,dbdq);
            gsl_matrix_set(_dm,2,0,dtdq);
            gsl_matrix_set(_dm,2,1,dsdq);
            gsl_matrix_set(_dm,2,2,dq2);

            toReturn = dbds - deriv_mult_aTm_1b(_da,_dm,_db);
        } else if (j == ChargeIdx::Q) {
            gsl_vector_set(_da,0,dtdq);
            gsl_vector_set(_da,1,dsdq);
            gsl_vector_set(_da,2,dq2);

            gsl_vector_set(_db,0,dtdb);
            gsl_vector_set(_db,1,db2);
            gsl_vector_set(_db,2,dbds);

            gsl_matrix_set(_dm,0,0,dt2);
            gsl_matrix_set(_dm,0,1,dtds);
            gsl_matrix_set(_dm,0,2,dtdq);
            gsl_matrix_set(_dm,1,0,dtdb);
            gsl_matrix_set(_dm,1,1,dbds);
            gsl_matrix_set(_dm,1,2,dbdq);
            gsl_matrix_set(_dm,2,0,dtds);
            gsl_matrix_set(_dm,2,1,ds2);
            gsl_matrix_set(_dm,2,2,dsdq);

            toReturn = dbdq - deriv_mult_aTm_1b(_da,_dm,_db);
        } else {
            std::cout << "Error calculating derivative term 4" << std::endl;
        }
    } else if (i == ChargeIdx::S) {
        if(j == ChargeIdx::B) {
            gsl_vector_set(_da,0,dtdb);
            gsl_vector_set(_da,1,db2);
            gsl_vector_set(_da,2,dbdq);

            gsl_vector_set(_db,0,dtds);
            gsl_vector_set(_db,1,ds2);
            gsl_vector_set(_db,2,dsdq);

            gsl_matrix_set(_dm,0,0,dt2);
            gsl_matrix_set(_dm,0,1,dtdb);
            gsl_matrix_set(_dm,0,2,dtdq);
            gsl_matrix_set(_dm,1,0,dtds);
            gsl_matrix_set(_dm,1,1,dbds);
            gsl_matrix_set(_dm,1,2,dsdq);
            gsl_matrix_set(_dm,2,0,dtdq);
            gsl_matrix_set(_dm,2,1,dbdq);
            gsl_matrix_set(_dm,2,2,dq2);

            toReturn = dbds - deriv_mult_aTm_1b(_da,_dm,_db);
        } else if (j == ChargeIdx::S) {
            gsl_vector_set(_da,0,dtds);
            gsl_vector_set(_da,1,dbds);
            gsl_vector_set(_da,2,dsdq);

            gsl_vector_set(_db,0,dtds);
            gsl_vector_set(_db,1,dbds);
            gsl_vector_set(_db,2,dsdq);

            gsl_matrix_set(_dm,0,0,dt2);
            gsl_matrix_set(_dm,0,1,dtdb);
            gsl_matrix_set(_dm,0,2,dtdq);
            gsl_matrix_set(_dm,1,0,dtdb);
            gsl_matrix_set(_dm,1,1,db2);
            gsl_matrix_set(_dm,1,2,dbdq);
            gsl_matrix_set(_dm,2,0,dtdq);
            gsl_matrix_set(_dm,2,1,dbdq);
            gsl_matrix_set(_dm,2,2,dq2);

            toReturn = ds2 - deriv_mult_aTm_1b(_da,_dm,_db);
        } else if (j == ChargeIdx::Q) {
            gsl_vector_set(_da,0,dtdq);
            gsl_vector_set(_da,1,dbdq);
            gsl_vector_set(_da,2,dq2);

            gsl_vector_set(_db,0,dtds);
            gsl_vector_set(_db,1,dbds);
            gsl_vector_set(_db,2,ds2);

            gsl_matrix_set(_dm,0,0,dt2);
            gsl_matrix_set(_dm,0,1,dtdb);
            gsl_matrix_set(_dm,0,2,dtdq);
            gsl_matrix_set(_dm,1,0,dtdb);
            gsl_matrix_set(_dm,1,1,db2);
            gsl_matrix_set(_dm,1,2,dbdq);
            gsl_matrix_set(_dm,2,0,dtds);
            gsl_matrix_set(_dm,2,1,dbds);
            gsl_matrix_set(_dm,2,2,dsdq);

            toReturn = dsdq - deriv_mult_aTm_1b(_da,_dm,_db);
        } else {
            std::cout << "Error calculating derivative term 4" << std::endl;
        }
    } else if (i == ChargeIdx::Q) {
        if(j == ChargeIdx::B) {
            gsl_vector_set(_da,0,dtdb);
            gsl_vector_set(_da,1,db2);
            gsl_vector_set(_da,2,dbds);

            gsl_vector_set(_db,0,dtdq);
            gsl_vector_set(_db,1,dsdq);
            gsl_vector_set(_db,2,dq2);

            gsl_matrix_set(_dm,0,0,dt2);
            gsl_matrix_set(_dm,0,1,dtdb);
            gsl_matrix_set(_dm,0,2,dtds);
            gsl_matrix_set(_dm,1,0,dtds);
            gsl_matrix_set(_dm,1,1,dbds);
            gsl_matrix_set(_dm,1,2,ds2);
            gsl_matrix_set(_dm,2,0,dtdq);
            gsl_matrix_set(_dm,2,1,dbdq);
            gsl_matrix_set(_dm,2,2,dsdq);

            toReturn = dbdq - deriv_mult_aTm_1b(_da,_dm,_db);
        } else if (j == ChargeIdx::S) {
            gsl_vector_set(_da,0,dtds);
            gsl_vector_set(_da,1,dbds);
            gsl_vector_set(_da,2,ds2);

            gsl_vector_set(_db,0,dtdq);
            gsl_vector_set(_db,1,dbdq);
            gsl_vector_set(_db,2,dq2);

            gsl_matrix_set(_dm,0,0,dt2);
            gsl_matrix_set(_dm,0,1,dtdb);
            gsl_matrix_set(_dm,0,2,dtds);
            gsl_matrix_set(_dm,1,0,dtdb);
            gsl_matrix_set(_dm,1,1,db2);
            gsl_matrix_set(_dm,1,2,dbds);
            gsl_matrix_set(_dm,2,0,dtdq);
            gsl_matrix_set(_dm,2,1,dbdq);
            gsl_matrix_set(_dm,2,2,dsdq);

            toReturn = dsdq - deriv_mult_aTm_1b(_da,_dm,_db);
        } else if (j == ChargeIdx::Q) {
            gsl_vector_set(_da,0,dtdq);
            gsl_vector_set(_da,1,dbdq);
            gsl_vector_set(_da,2,dsdq);

            gsl_vector_set(_db,0,dtdq);
            gsl_vector_set(_db,1,dbdq);
            gsl_vector_set(_db,2,dsdq);

            gsl_matrix_set(_dm,0,0,dt2);
            gsl_matrix_set(_dm,0,1,dtdb);
            gsl_matrix_set(_dm,0,2,dtds);
            gsl_matrix_set(_dm,1,0,dtdb);
            gsl_matrix_set(_dm,1,1,db2);
            gsl_matrix_set(_dm,1,2,dbds);
            gsl_matrix_set(_dm,2,0,dtds);
            gsl_matrix_set(_dm,2,1,dbds);
            gsl_matrix_set(_dm,2,2,ds2);

            toReturn = dq2 - deriv_mult_aTm_1b(_da,_dm,_db);
        } else {
            std::cout << "Error calculating derivative term 4" << std::endl;
        }
    } else {
        std::cout << "Error calculating derivative term 4" << std::endl;
    }

    return toReturn;
}

double EquationOfState::deriv_mult_aTm_1b(gsl_vector* a, gsl_matrix* m, gsl_vector* b)
{
  // baryon-only safety: if S/Q sector is zero, LU inversion can become singular.
  // Add a small epsilon to diagonal entries that are too small.
  const double EPS = 1e-3; // your requested ~0.001
  const double m11 = gsl_matrix_get(m, 1, 1);
  const double m22 = gsl_matrix_get(m, 2, 2);
  if (std::abs(m11) < EPS) gsl_matrix_set(m, 1, 1, (m11 >= 0.0 ? EPS : -EPS));
  if (std::abs(m22) < EPS) gsl_matrix_set(m, 2, 2, (m22 >= 0.0 ? EPS : -EPS));

  gsl_permutation_init(_dp);
  int s = -1;

  // Compute the LU decomposition of this matrix
  gsl_linalg_LU_decomp(m, _dp, &s);

	gsl_set_error_handler_off();

  // Compute the inverse of the LU decomposition
  gsl_matrix_set_zero(_dminv);
  int inversion_status = gsl_linalg_LU_invert(m, _dp, _dminv);

	if ( inversion_status )	// if an error occurred
	{
		cout << "Current TBQS location: "
				<< hbarc_MeVfm*T() << "   " << hbarc_MeVfm*muB() << "   "
				<< hbarc_MeVfm*muQ() << "   " << hbarc_MeVfm*muS() << endl << endl;

		cout << "Current EoS data:" << endl;
		cout << "pVal = " << pVal << endl
			 << "BVal = " << BVal << endl
			 << "SVal = " << SVal << endl
			 << "QVal = " << QVal << endl
			 << "eVal = " << eVal << endl
			 << "cs2Val = " << cs2Val << endl
			 << "db2 = " << db2 << endl
			 << "ds2 = " << ds2 << endl
			 << "dq2 = " << dq2 << endl
			 << "dt2 = " << dt2 << endl
			 << "dbdq = " << dbdq << endl
			 << "dbds = " << dbds << endl
			 << "dsdq = " << dsdq << endl
			 << "dtdb = " << dtdb << endl
			 << "dtds = " << dtds << endl
			 << "dtdq = " << dtdq << endl
			 << "entrVal = " << entrVal << endl << endl;


		cout << "a=" << endl;
		for (int ii = 0; ii < 3; ii++)
			cout << gsl_vector_get(a, ii) << "   ";
		cout << endl;

		cout << endl << "b=" << endl;
		for (int ii = 0; ii < 3; ii++)
			cout << gsl_vector_get(b, ii) << "   ";
		cout << endl;

		cout << endl << "m=" << endl;
		for (int ii = 0; ii < 3; ii++)
		{
			for (int jj = 0; jj < 3; jj++)
				cout << gsl_matrix_get(m, ii, jj) << "   ";
			cout << endl;
		}
		cout << endl << "minv=" << endl;
		for (int ii = 0; ii < 3; ii++)
		{
			for (int jj = 0; jj < 3; jj++)
				cout << gsl_matrix_get(_dminv, ii, jj) << "   ";
			cout << endl;
		}
		cout << endl;
		cout << "Exiting!" << endl;
		exit (-1);
	}
	gsl_set_error_handler (NULL);


    // Compute _dy = m^-1 @ b
    //gsl_blas_dgemv(CblasNoTrans,1,m,b,0,_dy);
    gsl_blas_dgemv(CblasNoTrans,1,_dminv,b,0,_dy);


	/*cout << endl << endl;
	cout << "=============================================" << endl;
	cout << "y=" << endl;
	for (int ii = 0; ii < 3; ii++)
		cout << gsl_vector_get(_dy, ii) << "   ";
	cout << endl;
	cout << "=============================================" << endl;
	cout << endl << endl;*/


    double toReturn = 0;
    //compute toReturn = aT @ _dy
    gsl_blas_ddot(a,_dy,&toReturn);

    return toReturn;
}