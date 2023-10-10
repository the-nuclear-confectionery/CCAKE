#include "eom_default.h"

using namespace ccake;

/// @brief Calculates the gamma factor = u^0.
/// @details This is the default implementation of the gamma factor calculation.
/// It assumes that the last component of the velocity vector is the longitudinal
/// velocity, and that the metric is Milne.
/// @param u The velocity vector (space components only).
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return 
template<unsigned int D>
double EoM_default<D>::gamma_calc(double *u, const double &time_squared) {
    return sqrt(1.0+dot(u,u,time_squared));
}

/// @brief Computes the inner product of two vectors.
/// @details This considers only the space components of the vectors. This is the 
/// implementation for the special case D=2.
/// @param u The first vector
/// @param v The second vector
/// @param time_squared The square of the time where the gamma factor will be computed
/// @return u^i v^i
template<>
double EoM_default<2>::dot(double *v, double *u, const double &time_squared) {
  double s = 0;
  for (unsigned int i=0; i<2; i++)
    s+= u[i]*v[i];
  return s;
}

/// @brief Transforms a scalar from the lab frame to the LRF.
/// @tparam D 
/// @param lab The quantity in the lab (computational) frame.
/// @param gamma Lorentz contraction factor. 
/// @param time_squared The square of the time where the gamma factor will be computed.
/// @return The quantity in the fluid LRF (local rest frame).
template<unsigned int D>
double EoM_default<D>::get_LRF(const double &lab, const double &gamma, 
                               const double &time_squared) {
                                return lab/gamma/time_squared;
}

/// @brief Computes the inner product of two vectors.
/// @details This considers only the space components of the vectors. This is the 
/// general case. It assumes that the last component of the velocity vector is the
/// longitudinal velocity.
/// @param u The first vector
/// @param v The second vector
/// @param time_squared The square of the time where the gamma factor will be computed
/// @return u^i v^i
template<unsigned int D>
double EoM_default<D>::dot(double *v, double *u, const double &time_squared) {
  double s = 0;
  for (unsigned int i=0; i<D-1; i++)
    s+= u[i]*v[i];
  s += u[D-1]*v[D-1]*time_squared;
  return s;
}

/// @brief Contracts a vector with a vector (LHS contraction).
/// @details We are assuming a contraction in the form x^j = v_i T^{i j}. Bear in
/// mind that we are considering only the space components.
/// @param v The vector
/// @param T The tensor
/// @param time_squared The square of the time where the gamma factor will be computed
/// @param x Vector where the result will be stored
template<unsigned int D>
void EoM_default<D>::dot(double (*v)[D],double (*T)[D][D], const double &time_squared, double (*x)[D]) {
  for (unsigned int j=0; j<D; j++){
    (*x)[j] = 0;
    for (unsigned int i=0; i<D-1; i++){
      (*x)[j]+= (*v)[i]*(*T)[i][j];
    }
    (*x)[j] += (*v)[D-1]*(*T)[D-1][j]*time_squared;
  }
}

/// @brief Contracts a vector with a vector (LHS contraction).
/// @details We are assuming a contraction in the form x^j = v_i T^{i j}. Bear in
/// mind that we are considering only the space components.
/// @param v The vector
/// @param T The tensor
/// @param time_squared The square of the time where the gamma factor will be computed
/// @param x Vector where the result will be stored
template<>
void EoM_default<2>::dot(double (*v)[2],double (*T)[2][2], const double &time_squared, double (*x)[2]) {
  for (unsigned int j=0; j<2; j++){
    (*x)[j] = 0;
    for (unsigned int i=0; i<2; i++){
      (*x)[j]+= (*v)[i]*(*T)[i][j];
    }
  }
}

/// @brief Ensures that the shear tensor will be traceless.
/// @details This is the default implementation of the traceless condition. The last 
/// component is assumed to be the longitudinal one and is modified as to ensure the tensor
/// is traceless.
/// @param pi_diag A D+1 dimensional array containing the diagonal elements of the shear tensor.
/// @return pi^{DD} = pi^{00} - \sum_{i=1}^{D-1} pi^{ii}/tau^2
template<unsigned int D>
double EoM_default<D>::get_shvDD(double* pi_diag, const double &time_squared){
    double s = pi_diag[0];
    for (unsigned int i=1; i<D; i++)
        s -= pi_diag[i];
    s /= time_squared;
    return s;
}

/// @brief Ensures that the shear tensor will be traceless.
/// @details This is the implementation for the special case (2+1)D. It assumes the existence of
/// a longitudinal component. We are using Milne coordinates.
/// @param pi_diag A D+1 dimensional array containing the diagonal elements of the shear tensor.
/// @return pi^{33} = pi^{00} - \sum_{i=1}^{D} pi^{ii}/\tau^2
/// \todo I may be wrong about this implementation. It is worth to double check.
template<>
double EoM_default<2>::get_shvDD(double* pi_diag, const double &time_squared){
    double s = pi_diag[0];
    for (unsigned int i=1; i<3; i++)
        s -= pi_diag[i];
    s /= time_squared;
    return s;
}

/// @brief Get the distance between two particles.
/// @details This is the default implementation of the distance calculation between two
/// particles. This is done in Milne coordinates.
/// @param x1 The position of the first particle.
/// @param x2 The position of the second particle.
/// @param time_squared The square of the time where the distance will be computed.
/// @return The distance between the two particles.
template<unsigned int D>
double EoM_default<D>::get_distance(double* x1, double* x2, const double &time_squared){
    double s = 0;
    for (unsigned int i=0; i<D-1; i++)
        s += (x1[i]-x2[i])*(x1[i]-x2[i]);
    s += (x1[D-1]-x2[D-1])*(x1[D-1]-x2[D-1])*time_squared;
    return sqrt(s);
}

template<unsigned int D>
void EoM_default<D>::evaluate_time_derivatives( Cabana::AoSoA<CabanaParticle, DeviceType, VECTOR_LENGTH> &particles )
{
  CREATE_VIEW( device_, particles );
  auto single_particle_update = KOKKOS_LAMBDA(int const iparticle){

      //1) Cache locally the particle data
      double sigma = device_hydro_scalar(iparticle, ccake::hydro_info::sigma);
      double gamma = device_hydro_scalar(iparticle, ccake::hydro_info::gamma);
      double t = device_hydro_scalar(iparticle, ccake::hydro_info::t);
      double t_squared = t*t;
      double Temperature = device_thermo(iparticle, ccake::thermo_info::T);
      double dwds = device_thermo(iparticle, ccake::thermo_info::dwds);
      double Bulk = device_hydro_scalar(iparticle, ccake::hydro_info::Bulk);
      double w = device_thermo(iparticle, ccake::thermo_info::w);
      double setas = device_hydro_scalar(iparticle, ccake::hydro_info::setas);
      double stauRelax = device_hydro_scalar(iparticle, ccake::hydro_info::stauRelax);
      double s = device_thermo(iparticle, ccake::thermo_info::s);
      double zeta = device_hydro_scalar(iparticle, ccake::hydro_info::zeta);
      double tauRelax = device_hydro_scalar(iparticle, ccake::hydro_info::tauRelax);
      double dwdB = device_thermo(iparticle, ccake::thermo_info::dwdB);
      double rhoB = device_thermo(iparticle, ccake::thermo_info::rhoB);
      double dwdS = device_thermo(iparticle, ccake::thermo_info::dwdS);
      double rhoS = device_thermo(iparticle, ccake::thermo_info::rhoS);
      double dwdQ = device_thermo(iparticle, ccake::thermo_info::dwdQ);
      double rhoQ = device_thermo(iparticle, ccake::thermo_info::rhoQ);
      double shv33 = device_hydro_scalar(iparticle, ccake::hydro_info::shv33);


      double gradV[D][D];
      double v[D];
      double shv[D+1][D+1];
      double u[D];
      double minshv[D];
      double gradshear[D];
      double gradP[D];
      double gradBulk[D];
      double divshear[D];
      double du_dt[D];

      for (unsigned int idir=0; idir<D; idir++){
        v[idir] = device_hydro_vector(iparticle, ccake::hydro_info::v, idir);
        u[idir] = device_hydro_vector(iparticle, ccake::hydro_info::u, idir);
        minshv[idir] = device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, 0, idir+1);
        gradshear[idir] = device_hydro_vector(iparticle, ccake::hydro_info::gradshear, idir);
        gradP[idir] = device_hydro_vector(iparticle, ccake::hydro_info::gradP, idir);
        gradBulk[idir] = device_hydro_vector(iparticle, ccake::hydro_info::gradBulk, idir);
        divshear[idir] = device_hydro_vector(iparticle, ccake::hydro_info::divshear, idir);
        du_dt[idir] = 0.;
        for (unsigned int jdir=0; jdir<D; jdir++)
          gradV[idir][jdir] = device_hydro_space_matrix(iparticle, ccake::hydro_info::gradV, idir, jdir);
      }
      for (unsigned int idir=0; idir<D+1; idir++)
      for (unsigned int jdir=0; jdir<D+1; jdir++)
          shv[idir][jdir] = device_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, idir, jdir);
      
      //Will be filled later
      double aux_vector[D];
      double F[D];
      double gradU[D][D];
      double partU[D][D];
      double uu[D][D];
      double pi_u[D][D];
      double pimin[D][D];
      double piutot[D][D];
      double M[D][D];
      double M_inverse[D][D];
      double sub[D][D];
      double dshv_dt[D][D];


      //2) Perform computations
      double dsigma_dt = -(gradV[0][0]+gradV[1][1])*sigma;
      double gamma_squared = gamma*gamma;
      double gamma_cube = gamma*gamma_squared;
      double gamma_tau = gamma*t;
      double dwdsT    = dwds/Temperature;
      double dwdsT1 = 1.0 - dwdsT;
      double sigl = dsigma_dt/sigma - 1/t;
      double bigPi = Bulk*sigma/gamma_tau;
      double C = w + bigPi;
      

      dot(&v, &gradV, t_squared, &aux_vector);
      for (int idir=0; idir<D; idir++)
      for (int jdir=0; jdir<D; jdir++){
        uu[idir][jdir] = u[idir]*u[jdir];
        gradU[idir][jdir] = gamma*gradV[idir][jdir] + gamma_cube*aux_vector[idir]*v[jdir];
        pi_u[idir][jdir] = shv[0][idir+1]*u[jdir];
        pimin[idir][jdir] = shv[idir+1][jdir+1];
      }
      for (int idir=0; idir<D; idir++)
      for (int jdir=0; jdir<D; jdir++){
        piutot[idir][jdir] = pi_u[idir][jdir] + pi_u[jdir][idir];	
        partU[idir][jdir] = gradU[idir][jdir] + gradU[jdir][idir];
      }
      
      double bsub = 0;
      for (int idir=0; idir<D; idir++)
      for (int jdir=0; jdir<D; jdir++)
        bsub += gradU[idir][jdir]*( pimin[idir][jdir] 
                                    + uu[idir][jdir]*shv[0][0]/gamma_squared 
                                    - piutot[idir][jdir]/gamma );
      
      
      //eta_o_tau
      ///TODO: Disable if shear is not used - Equivalent o setting shear to 0?
      double eta_o_tau = setas/stauRelax;

      double Agam = w - dwds*(s + bigPi/Temperature) - zeta/tauRelax
                      - dwdB*rhoB - dwdS*rhoS - dwdQ*rhoQ; 
      double Agam2 = ( Agam - eta_o_tau/3.0 - dwdsT1*shv[0][0] ) / gamma;
      double Ctot = C + eta_o_tau*(1.0/gamma_squared-1.0);

      double Btot = ( Agam*gamma + 2.0*eta_o_tau/3.0*gamma )*sigl ///TODO: 2/3 or 1/3. See Jaki's (Eq. 274)?
                      + bigPi/tauRelax
                      + dwdsT*( gamma_tau*shv33 + bsub );

      double gamt = 0.0, pre = 0.0, p1 = 0.0;
      ///TODO: Wrap next 3 lines in if statement
      //if ( this->settingsPtr->using_shear )
      //{
        gamt = 1.0/gamma/stauRelax;
        pre  = eta_o_tau/gamma;
        p1   = gamt - 4.0/3.0/sigma*dsigma_dt + 1.0/t/3.0;
      //}

      // set the Mass and the Force
      for (int idir=0; idir<D; idir++){
        for (int jdir=0; jdir<D; jdir++){
          double Imat = idir == jdir ? 1.0 : 0.0;
          M[idir][jdir] = -Agam2*uu[idir][jdir] + Ctot*Imat
                          -(1+4./3./gamma_squared)*pi_u[idir][jdir]
                          + dwdsT1*pi_u[jdir][idir] + gamma*pimin[idir][jdir];
        }
        F[idir] = Btot*u[idir] + gradshear[idir] - ( gradP[idir] + gradBulk[idir] + divshear[idir] );
      }

      //===============
      // shear contribution
      ///TODO: Wrap on if for shear
      //if ( this->settingsPtr->using_shear )
      dot(&v, &partU, t_squared, &aux_vector);
      for (int idir=0; idir<D; idir++)
        F[idir] += pre*aux_vector[idir] + p1*minshv[idir];

      
      Utilities::inverse<D>(&M,&M_inverse);

      //===============
      // compute acceleration
      for (int idir=0; idir<D; idir++){
        for (int jdir=0; jdir<D; jdir++){
          du_dt[idir] += F[jdir]*M_inverse[idir][jdir];
        }
      }

      //===============
      // "coordinate" divergence
      double div_u = (1./ gamma)*dot(u, du_dt, t_squared)
                      - ( gamma/ sigma ) * dsigma_dt;

      //===============
      // "covariant" divergence
      double bigtheta = div_u*t + gamma;

      //===============
      for (int idir=0; idir<D; idir++){
        aux_vector[idir] = minshv[idir] - shv[0][0]*v[idir];
        for (int jdir=0; jdir<D; jdir++)
          sub[idir][jdir] = pimin[idir][jdir] + (shv[0][0]/gamma_squared)*uu[idir][jdir]
                            - piutot[idir][jdir]/gamma;
          
      }
      double aux = 0.;
      for (int idir=0; idir<D; idir++)
        aux += sub[idir][idir]*gradU[idir][idir];
      double inside = t*dot( aux_vector, du_dt, t_squared ) - aux - gamma*t*shv33;

      double d_dt_specific_s = 1./sigma/Temperature*( -bigPi*bigtheta + inside );
      

      //formulating simple setup for Beta_Bulk derivative  //WMS: I don't know what is this
      //hi.finite_diff_cs2 = (ti.cs2 - hi.prev_cs2)/0.05; // Asadek
	    //hi.finite_diff_T   = (ti.T - hi.prev_T)/0.05;     // Asadek
	    //hi.finite_diff_w   = (ti.w - hi.prev_w)/0.05;     // Asadek
	    //hi.dBeta_dt        = 0.5*((-hi.finite_diff_T/(ti.T*ti.T))*(1/ti.w)*(1/((1/3-ti.cs2)*(1/3-ti.cs2))))
	    //                   + 0.5*((-hi.finite_diff_w/(ti.w*ti.w))*(1/ti.T)*(1/((1/3-ti.cs2)*(1/3-ti.cs2))))
			//	          	     + 0.5*((4*ti.cs2*hi.finite_diff_cs2*(1/((1/3-ti.cs2)*(1/3-ti.cs2)*(1/3-ti.cs2))))*(1/ti.T)*(1/ti.w));//Asadek 

      double dBulk_dt = -(zeta*bigtheta/sigma + Bulk/gamma)/tauRelax;

      //===============
      // define auxiliary variables
      double vduk = dot(v, du_dt, t_squared);
      double ulpi[D][D], Ipi[D][D], ududt[D][D], vsub[D][D];
      for (int idir=0;idir<D;++idir)
      for (int jdir=0;jdir<D;++jdir){
        double Imat = idir == jdir ? 1 : 0;
        ulpi[idir][jdir] = u[idir]*shv[0][jdir];
        Ipi[idir][jdir] = - 2.*eta_o_tau*(Imat + uu[idir][jdir] )/3. 
                          + 4.*pimin[idir][jdir]/3.;
        ududt[idir][jdir] = u[idir]*du_dt[jdir];
        vsub[idir][jdir] = 0;
        for (int kdir=0;kdir<D;++kdir)
          vsub[idir][jdir] += (u[idir]*pimin[jdir][kdir]
                           + u[jdir]*pimin[idir][kdir])*du_dt[kdir];
      }

      //===============
      ///TODO: - ADD READABLE TERM NAMES
      ///TODO: Wrap on if for shear
      //if ( this->settingsPtr->using_shear )
      for (int idir=0;idir<D;++idir)
      for (int jdir=0;jdir<D;++jdir)
        dshv_dt[idir][jdir] = - gamt*( pimin[idir][jdir] + setas*partU[idir][jdir] )
                              - eta_o_tau*( ududt[idir][jdir] + ududt[jdir][idir] )
                              + vsub[idir][jdir] + sigl*Ipi[idir][jdir]
                              - vduk*( ulpi[idir][jdir] + ulpi[jdir][idir] + Ipi[idir][jdir]/gamma );
      
    //OK UNTIL ABOVE HERE
    
    //3) Update the particle data
    ///TODO: Do I really need to update all of these?
    device_hydro_scalar(iparticle, ccake::hydro_info::dsigma_dt) = dsigma_dt;
    device_hydro_scalar(iparticle, ccake::hydro_info::gamma_squared) = gamma_squared;
    device_hydro_scalar(iparticle, ccake::hydro_info::gamma_cube) = gamma_cube;
    device_hydro_scalar(iparticle, ccake::hydro_info::gamma_tau) = gamma_tau;
    device_hydro_scalar(iparticle, ccake::hydro_info::dwdsT1) = dwdsT1;
    device_hydro_scalar(iparticle, ccake::hydro_info::sigl) = sigl;
    device_hydro_scalar(iparticle, ccake::hydro_info::bigPI) = bigPi;
    device_hydro_scalar(iparticle, ccake::hydro_info::C) = C;
    device_hydro_scalar(iparticle, ccake::hydro_info::eta_o_tau) = eta_o_tau;
    device_hydro_scalar(iparticle, ccake::hydro_info::Agam) = Agam;
    device_hydro_scalar(iparticle, ccake::hydro_info::Agam2) = Agam2;
    device_hydro_scalar(iparticle, ccake::hydro_info::Ctot) = Ctot;
    device_hydro_scalar(iparticle, ccake::hydro_info::Btot) = Btot;
    device_hydro_scalar(iparticle, ccake::hydro_info::div_u) = div_u;
    device_hydro_scalar(iparticle, ccake::hydro_info::bigtheta) = bigtheta;
    device_hydro_scalar(iparticle, ccake::hydro_info::inside) = inside;
    device_hydro_scalar(iparticle, ccake::hydro_info::dBulk_dt) = dBulk_dt;
    device_d_dt_spec(iparticle, ccake::densities_info::s) = d_dt_specific_s;



    for (int idir=0; idir<D; idir++){
      device_hydro_vector(iparticle, ccake::hydro_info::du_dt, idir) = du_dt[idir];
      for (int jdir=0; jdir<D; jdir++){
        device_hydro_space_matrix(iparticle, ccake::hydro_info::gradU, idir, jdir) = gradU[idir][jdir];
        device_hydro_space_matrix(iparticle, ccake::hydro_info::dshv_dt, idir, jdir) = dshv_dt[idir][jdir];  
      }
    }
  };

  Kokkos::parallel_for( "single_particle_update", particles.size(),
                        single_particle_update );
}
