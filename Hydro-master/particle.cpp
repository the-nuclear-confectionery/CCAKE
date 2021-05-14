#include "vector.h"
#include "matrix.h"
#include "eos.h"
#include "particle.h"
#include <stdio.h>
#include <math.h>

template <int D>
eos Particle<D>::EOS;

// Functions added by Christopher Plumberg
template <int D>
double Particle<D>::EOST() { return particle_T; }

template <int D>
double Particle<D>::EOSmuB() { return particle_muB; }

template <int D>
double Particle<D>::EOSmuS() { return particle_muS; }

template <int D>
double Particle<D>::EOSmuQ() { return particle_muQ; }



template <int D>
double Particle<D>::EOSp()
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.p();
}


// Functions added by Christopher Plumberg
template <int D>
double Particle<D>::EOSs()
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.s();
}



template <int D>
double Particle<D>::EOSe()
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.e();
}



template <int D>
double Particle<D>::EOSw()
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.w();
}



template <int D>
double Particle<D>::EOSA()
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.A();
}



template <int D>
double Particle<D>::EOSdwds()
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.dwds();
}



template <int D>
double Particle<D>::EOSs_terms_T(double Tin)
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.s_terms_T(Tin);
}



template <int D>
double Particle<D>::EOSs_out(double e_In, double rhoB_In, double rhoS_In, double rhoQ_In)
{
	double sVal = EOS.s_out( e_In, rhoB_In, rhoS_In, rhoQ_In );
	particle_T = EOS.T();
	particle_muB = EOS.muB();
	particle_muS = EOS.muS();
	particle_muQ = EOS.muQ();

	return sVal;
}


template <int D>
double Particle<D>::EOSs_out(double e_In)
{
	double sVal = EOS.s_out( e_In, 0.0, 0.0, 0.0 );
	particle_T = EOS.T();
	particle_muB = EOS.muB();
	particle_muS = EOS.muS();
	particle_muQ = EOS.muQ();

	return sVal;
}


template <int D>
void Particle<D>::EOSupdate_s(double s_In, double rhoB_In, double rhoS_In, double rhoQ_In)
{
	EOS.update_s( s_In, rhoB_In, rhoS_In, rhoQ_In );
	particle_T = EOS.T();
	particle_muB = EOS.muB();
	particle_muS = EOS.muS();
	particle_muQ = EOS.muQ();

	return;
}


template <int D>
void Particle<D>::EOSupdate_s(double s_In)
{
	EOS.update_s( s_In, 0.0, 0.0, 0.0 );
	particle_T = EOS.T();
	particle_muB = EOS.muB();
	particle_muS = EOS.muS();
	particle_muQ = EOS.muQ();

	return;
}


