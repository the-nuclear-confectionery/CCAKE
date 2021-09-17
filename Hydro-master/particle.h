#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "vector.h"
#include "matrix.h"
#include "eos.h"
#include <stdio.h>
#include <math.h>

//extern double freezeoutT;

// need to switch to dimension instead!
template <int D>
class Particle
{
private:


public:
////////////////////////////////////////////////////////////////////////////////
//                            Particle Variables                              //
////////////////////////////////////////////////////////////////////////////////


    int btrack;
    static eos EOS;	//use one copy of EOS for all particles
    double Agam, Agam2;
    double sigmaweight;        // specific volume per particle (times s_an)
	double rho_weight;		   // specific volume per particle (without s_an)
	double transverse_area;	   // dx * dy
    Vector<double,D> r;                   // position
    Vector<double,D> v;                   // velocity
    Vector<double,D> u;                   // relativistic velocity
    Vector<double,D> qmom;
    Vector<double,D> du_dt,gradsig;               // relativistic velocity derivative
    Vector<double,D> k1,k2,k3,k4,ksub;
    Vector<double,D> r1,r2,r3,r4,rsub;
    double ets1,ets2,ets3,ets4;
    double b1,b2,b3,b4;
    double bn1,bn2,bn3,bn4;
    Matrix<double,D,D> shv1,shv2,shv3,svh4;
    double shv33;
//    double vmag,vang;

    struct FRZ
    {
        Vector<double,D> r,u,gradP;
        double t,s,T,theta,bulk,sigma,shear33,inside;
        Matrix<double,D+1,D+1> shear;
    };

    FRZ frz1,frz2,fback,fback2,fback3,fback4;

//    double es1,es2,es3,es4;
    double div_u;              // four-divergence of relativistic velocity
    double gamma;               // Lorentz factor
    double gamcalc();
    //double gamcalcbsq();
    double s_sub,e_sub,s_an,s_rat,sigsub;
    double eta_sigma;          // Ratio entropy/especific volume
    double detasigma_dt;
    double Bulk;               // Bulk Viscosity weight
    double bigPI;        // total bulk viscosity
    double C;
    double tauRelax;           // Bulk Relaxation time
    double stauRelax;        // Shear Relxation time;
    double dBulk_dt;           // derivative Bulk Viscosity
    double zeta;               // bulk coefficient
    double setas;
    int count;
    int Freeze;
    double Ctot, Btot;
    Matrix<double,D+1,D+1> shv;
    Matrix<double,D,D> shv0,dshv_dt;
    Vector<double,D>  divshear,gradshear;

    double sv_eta,taupi;

////////////////////////////////////////////////////////////////////////////////
//                           Fluid Variables                                  //
////////////////////////////////////////////////////////////////////////////////
    //double sigmastar;
    double sigma;              // especific volume
    double dsigma_dt;          // derivative of especific volume
    //double div_J;              // divergence of the flux
    //Vector<double,D> flux;               // flux = sigma * velocity
    Vector<double,D> gradP;              // Gradient of Pressure
    Matrix<double,D,D> gradV,gradU;          // Gradient of velocity needed for shear

    Vector<double,D> gradBulk,gradrhoB,gradrhoS,gradrhoQ;           // Gradient of Bulk Viscosity
    Vector<double,D> gradsigma;          // Gradient of especific volume

    double dw_ds;              // derivative of the enthalpy on entropy
    double eta;                   // entropy density
    double rhoB_an, rhoB_sub;                   // Baryon density
    double rhoS_an, rhoS_sub;                   // strange density
    double rhoQ_an, rhoQ_sub;                   // electric charge density
//    double P;                   // pressure
//    double epsilon;               // energy density
//    double s;                   // entropy density
//    double T;                   // temperature
    double eden;

    double bigtheta;

	double B, S, Q;						// baryon, strange, electric charge
    double drhoB_dt,drhoS_dt,drhoQ_dt;

    Particle();
    void start(string enter, eos & EOS_in, bool setup_EOS = true);
    void calc( double tin );
    void calcbsq(double tin);
    void sigset(double tin);
    void vsigset(double tin);
    void svsigset(double tin, int i);
    void bsqsvsigset(double tin, int i);
    void returnv_A();
    void return_sv_A();
    void return_bsqsv_A();
    void returnA() ;
    void setvisc(int etaconst,double bvf, double svf, double zTc, double sTc, double sig, int type);
    void frzcheck(double tin,int &count, int N);
    double g2,g3,gt;
    double eta_o_tau,dwdsT1,sigl;
    Matrix <double,D,D> uu,pimin,piu,piutot;
    //void ucalc();
    double Bsub()  ;
    Matrix<double, D,D> Msub(int i),Imat  ;
    Matrix<double, D,D> dpidtsub();
    void sets(double tin2) ;
    double bcheck;
    double check;
    void setvar();
    double saves;
    double inside;

	double EOST();
	double EOSmuB();
	double EOSmuS();
	double EOSmuQ();

	double EOSp();
	double EOSs();
	double EOSe();
	double EOSB();
	double EOSS();
	double EOSQ();
	double EOSw();
	double EOSs_terms_T(double Tin);
	//double EOScs2out();
	//double EOSwfz();
	//double EOScs2out();
	//double EOSwfz();
	double EOSA();
	double EOSdwds();
	double EOSdwdB();
	double EOSdwdS();
	double EOSdwdQ();
	double EOSs_out(double e_In);
	double EOSs_out(double e_In, double rhoB_In, double rhoS_In, double rhoQ_In);
	void EOSupdate_s(double s_In);
	void EOSupdate_s(double s_In, double rhoB_In, double rhoS_In, double rhoQ_In);

	//double particle_T, particle_muB, particle_muS, particle_muQ;
	struct particle_thermo
    {
		double T, muB, muS, muQ;
        double p, s, B, S, Q, e, cs2, w, A;
		double dwds, dwdB, dwdS, dwdQ;
		//double db2, dq2, ds2, dbdq, dbds, dsdq, dtdb, dtdq, dtds, dt2;
    };
	particle_thermo SPH_cell;

};




template <int D>
Particle<D>::Particle()
{
    div_u = 0.0;
    gamma = 0.0;
    s_sub = 0.0;
    e_sub = 0.0;
    s_an = 0.0;
    s_rat = 0.0;
    sigsub = 0.0;
    eta_sigma = 0.0;
    detasigma_dt = 0.0;
    Bulk = 0.0;
    bigPI = 0.0;
    tauRelax = 0.0;
    stauRelax = 0.0;
    dBulk_dt = 0.0;
    zeta = 0.0;
    setas = 0.0;
    sv_eta = 0.0;
    taupi = 0.0;

    sigma = 0.0; 
    dsigma_dt = 0.0;

    dw_ds = 0.0; 
    eta = 0.0;
    rhoB_an = 0.0;
    rhoB_sub = 0.0;
    rhoS_an = 0.0;
    rhoS_sub = 0.0;
    rhoQ_an = 0.0;
    rhoQ_sub = 0.0;

    bigtheta = 0.0;

	B = 0.0;
    S = 0.0;
    Q = 0.0;
    drhoB_dt = 0.0;
    drhoS_dt = 0.0;
    drhoQ_dt = 0.0;

}

template <int D>
void Particle<D>::start(string enter, eos & EOS_in, bool setup_EOS)
{
	if (setup_EOS)
	{
		EOS = EOS_in;
	    EOS.eosin(enter);
	}
    Imat.identity();
	return;
}

template <int D>
double Particle<D>::gamcalc()
{
    return sqrt( Norm2(u) + 1 );

}


/*template <int D>
double Particle<D>::gamcalcbsq()
{
    return sqrt( Norm2(u) + 1 );

}*/

template <int D>
void Particle<D>::frzcheck(double tin,int &count, int N)
{
    if( Freeze==0)
    {
        if (EOST()<=freezeoutT)
        {
            Freeze=1;
            frz2.t=tin;

        }
    }
    else if( Freeze==1)
    {


        if (btrack==-1) {
            count +=1;
            Freeze=3;
            frz1.t=tin;
        }
        else if (EOST()>frz1.T) {
            Freeze=1;
            frz2.t=tin;
        }
        else if(EOST()<=freezeoutT)
        {
            count +=1;
            Freeze=3;
            frz1.t=tin;
        }
        else
        {
            Freeze=0;
        }
    }


}



//  Computes gamma and velocity
template <int D>
void Particle<D>::calc(double tin)
{

    gamma=gamcalc();
    v =(1/gamma)*u;
    double s_in2= eta/gamma/tin;
    qmom=((EOSe()+ EOSp())*gamma/sigma)*u;
    EOSupdate_s(s_in2);    // single-argument version
}

//  Computes gamma and velocity
template <int D>
void Particle<D>::calcbsq(double tin)
{
    gamma=gamcalc();
    v =(1/gamma)*u;
    double s_in2= eta/gamma/tin;
    qmom=((EOSe()+ EOSp())*gamma/sigma)*u;
	double rhoB_in2 = B*sigma/sigmaweight;
	double rhoS_in2 = S*sigma/sigmaweight;
	double rhoQ_in2 = Q*sigma/sigmaweight;
	//double rhoB_in2 = B*sigma/sigmaweight/gamma/tin;
	//double rhoS_in2 = S*sigma/sigmaweight/gamma/tin;
	//double rhoQ_in2 = Q*sigma/sigmaweight/gamma/tin;
	rhoB_an = rhoB_in2;
	rhoS_an = rhoS_in2;
	rhoQ_an = rhoQ_in2;
cout << "\t - finding EoS solution for sBSQ: " << s_in2 << "   "
		<< rhoB_in2 << "   " << rhoS_in2 << "   " << rhoQ_in2 << endl;
	EOSupdate_s(s_in2, rhoB_in2, rhoS_in2, rhoQ_in2);
}

template <int D>
void Particle<D>::returnA()
{
    Agam=EOSA()/gamma;
}

template <int D>
void Particle<D>::returnv_A()
{
	//if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    Agam=EOS.w()-EOSdwds()*(EOSs()+ bigPI/EOST() )- zeta/tauRelax ;
    Agam/=gamma;
}

template <int D>
void Particle<D>::return_sv_A()
{
	//if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    eta_o_tau=setas/stauRelax;

    Agam=EOSw()-EOSdwds()*(EOSs()+ bigPI/EOST() )- zeta/tauRelax ;

    Agam2=(Agam-eta_o_tau*(0.5-1/3.) -dwdsT1*shv.x[0][0])/gamma;
    Ctot=C+eta_o_tau*(1/g2-1)/2.;

}

template <int D>
void Particle<D>::return_bsqsv_A()
{
	//if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    eta_o_tau=setas/stauRelax;

    Agam=EOSw()-EOSdwds()*(EOSs()+ bigPI/EOST() )- zeta/tauRelax
			- EOS.dwdB() * rhoB_an
			- EOS.dwdS() * rhoS_an
			- EOS.dwdQ() * rhoQ_an;

    Agam2=(Agam-eta_o_tau*(0.5-1/3.) -dwdsT1*shv.x[0][0])/gamma;
    Ctot=C+eta_o_tau*(1/g2-1)/2.;
//	cout << "   " << eta_o_tau << "   " << EOSs() << "   " << Agam2
//			<< "   " << Ctot << "   " << EOS.dwdB() * rhoB_an
//			<< "   " << EOS.dwdS() * rhoS_an << "   " << EOS.dwdQ() * rhoQ_an;

}

template <int D>
void Particle<D>::sigset(double tin)
{
    Agam2=Agam*gamma*gamma*(dsigma_dt/sigma -1/tin);

}

template <int D>
void Particle<D>::vsigset(double tin)
{

    bigPI = Bulk*sigma/gamma/tin ;

    C=EOSw()+bigPI;

    returnv_A();
    Agam2=Agam*gamma*gamma*(dsigma_dt/sigma -1./tin)+ bigPI/tauRelax;

}

template <int D>
void Particle<D>::svsigset(double tin,int i)
{
    g2=gamma*gamma;
    g3=gamma*g2;
    gt=gamma*tin;
    double dwdsT=EOSdwds()/EOST();
    dwdsT1=1- EOSdwds()/EOST();
    sigl=dsigma_dt/sigma -1/tin;


    gradU=gamma*gradV+g3*(v*(v*gradV));

    bigPI= Bulk*sigma/gt ;

    C=EOSw()+ bigPI;


    return_sv_A();



    Btot=(Agam*gamma+eta_o_tau/3*gamma)*sigl+ bigPI/tauRelax + dwdsT* (gt*shv33+  Bsub());



    check=sigl;
}

template <int D>
double Particle<D>::Bsub()
{

    mini(pimin,shv);
    uu=u*u;
    piu=rowp1(0,shv)*u;
    piutot=piu+transpose(piu);
    double bsub=0.;
    double pig=shv.x[0][0]/g2;



    for (int i=0; i<=1; i++) {
        for (int j=0; j<=1; j++) {
            bsub+=(pimin.x[i][j]+pig*uu.x[j][i]-(piu.x[i][j] +piu.x[j][i])/gamma)*gradU.x[i][j];


        }
    }

    return bsub;
}

template <int D>
Matrix<double, D,D> Particle<D>::Msub(int i)
{


    piu=rowp1(0,shv)*u;

    Matrix<double,D,D> msub=Agam2*uu+Ctot*gamma*Imat-(1+4/3./g2)*piu+dwdsT1*transpose(piu)+gamma*pimin ;

cout << "CHECK Msub: " << Agam2 << "   "
		<< uu << "   "
		<< Ctot << "   "
		<< gamma << "   "
		<< Imat << "   "
		<< (1+4/3./g2) << "   "
		<< piu << "   "
		<< dwdsT1 << "   "
		<< transpose(piu) << "   "
		<< pimin << endl;

    //if (i==11000) cout << "A=" << Agam2 << " C=" << Ctot << endl;


    return msub;
}

template <int D>
Matrix<double, D,D> Particle<D>::dpidtsub()
{


    Matrix<double,D,D> vsub;


    for (int i=0; i<=1; i++) {
        for (int j=0; j<=1; j++) {
            for (int k=0; k<=1; k++) {
                vsub.x[i][j]+=(u.x[i]*pimin.x[j][k]+u.x[j]*pimin.x[i][k])*du_dt.x[k];
            }

        }
    }

    return vsub;
}

template <int D>
void Particle<D>::bsqsvsigset(double tin,int i)
{
	// from svsigset
    g2=gamma*gamma;
    g3=gamma*g2;
    gt=gamma*tin;
    double dwdsT=EOSdwds()/EOST();
    dwdsT1=1- EOSdwds()/EOST();
    sigl=dsigma_dt/sigma -1/tin;
    gradU=gamma*gradV+g3*(v*(v*gradV));
    bigPI= Bulk*sigma/gt ;
    C=EOSw()+ bigPI;
//	cout << "Check thermo1: " << i << "   " << tin << "   " << gamma
//			<< "   " << EOST()*197.327 << "   " << EOSdwds() << "   " << dwdsT1
//			<< "   " << sigma << "   " << dsigma_dt << "   " << sigl
//			<< "   " << gradU << "   " << gradV << "   " << v;
    return_bsqsv_A();
    Btot=(Agam*gamma+eta_o_tau/3*gamma)*sigl+ bigPI/tauRelax + dwdsT* (gt*shv33+  Bsub());
    check=sigl;
//	cout << "   " << EOSw() << "   " << bigPI << "   " << Agam
//			<< "   " << shv33 << "   " << Bsub() << "   " << Btot << endl;
	cout << "Check thermo: " << i << "   "
			<< tin << "   " << r << "   "
			<< SPH_cell.T << "   "
			<< SPH_cell.muB << "   "
			<< SPH_cell.muS << "   "
			<< SPH_cell.muQ << "   "
			<< SPH_cell.p << "   "
			<< SPH_cell.s << "   "
			<< SPH_cell.B << "   "
			<< SPH_cell.S << "   "
			<< SPH_cell.Q << "   "
			<< SPH_cell.e << "   "
			<< SPH_cell.cs2 << "   "
			<< SPH_cell.w << "   "
			<< SPH_cell.A << "   "
			<< SPH_cell.dwds << "   "
			<< SPH_cell.dwdB << "   "
			<< SPH_cell.dwdS << "   "
			<< SPH_cell.dwdQ << endl;

}

template <int D>
void Particle<D>::setvisc(int etaconst,double bvf, double svf, double zTc, double sTc, double sig, int type)
{
    //check variables!!!!

    if (type==1) // bulk viscosity
    {


        //double cscheck=EOScs2out(EOST());
        //zeta = bvf*0.5*zconst*(1/3.-cscheck);
        //zeta=-0.0295874+4.22749/sqrt(522.496+pow(-182.807+temp,2.))+5.09207833442629*pow(10,-7)*temp*temp;

        //double cscheck=EOScs2out(EOST());

        //double sig=5;
        //double mf=2.0768;

        double temp=EOST()*197.3;
        zeta =bvf/(sig*sqrt(2*PI))*exp(-pow(temp-zTc,2)/(2.*sig*sig));
        zeta*=EOSs();
        if (zeta<0.001) zeta=0.001;
        tauRelax = 9*zeta/(EOSe()-3*EOSp());
        if (tauRelax <0.1) tauRelax=0.1 ;

    }
    else if (type==2) // shear viscosity
    {

        setas =svf*0.08;

        stauRelax=5*setas/EOSw();

    }
    else if (type==3) // shear+bulk viscosity
    {

        if ((etaconst==1)||(etaconst==3)||(etaconst==4))
        {   // const eta/s
            setas=2*EOSs()*svf;  // svf defines eta/s const (the two is needed for the definition in the code, don't remove!
        }

        //    for TECHQM/Gubser set svf=0.08
    }
    else if (type==4) // BSQ+shear+bulk viscosity
    {

        if ((etaconst==1)||(etaconst==3)||(etaconst==4)) {// const eta/s
            setas=2*EOSs()*svf;  // svf defines eta/s const (the two is needed for the definition in the code, don't remove!



            //    for TECHQM/Gubser set svf=0.08
        }
        else { // eta/s (T)



//               if( temp>TC){
//                setas = svf*EOSs()*(-0.289 + 0.288*temp + 0.0818*temp*temp);
//            }
//            else {
//                setas = svf*EOSs()*(0.681 - 0.0594*temp - 0.544 *temp*temp);
//            }
//        if (setas<0.01) setas=0.001;


            if (etaconst==5) {
                double TC=173.9; // 173.9/197.3
                double temp=EOST()*197.3/TC;
                double TC0=149.4/TC; // 149.4/197.3
                if( temp<TC0) {
                    setas = EOSs()*(8.0191 - 16.4659* temp  +  8.60918* temp *temp  );
                }
                else if( temp>1) {
                    //setas = EOSs()*(0.48 - 0.36*temp );
                    setas = EOSs()*(0.007407515123054544 +  0.06314680923610914* temp + 0.08624567564083624* temp *temp  );
                }
                else {
                    //setas = EOSs()*(-0.107143 + 0.227143*temp);
                    setas = EOSs()*(0.397807 + 0.0776319* temp - 0.321513* temp *temp  );
                }
            }
            if (etaconst==6) {
                double TC=155; // 173.9/197.3
                double temp=EOST()*197.3/TC;
                double z=pow(0.66*temp,2);
                double alpha=33./(12.*PI)*(z-1)/(z*log(z));

                setas = EOSs()*(0.0416762/pow(alpha,1.6)+ 0.0388977/pow(temp,5.1) );

            }
            else {
                double TC=sTc/197.3;
                double temp=EOST()/TC;
                if( temp>TC) {
                    setas = EOSs()*(0.3153036437246963 + 0.051740890688259315* temp  -  0.24704453441295576* temp *temp  );
                }
                else {
                    setas = EOSs()*(0.0054395278010882795 + 0.08078575671572835*temp + 0.033774715483183566* temp *temp );
                }
            }








        }
        stauRelax=5*setas/EOSw();
        if (stauRelax <0.005) stauRelax=0.005 ;

        /// defining bulk viscosity

        if (bvf==0)
        {
            zeta =0;
            tauRelax = 1;
        }
        else
        {

            double temp=EOST()*197.3;

            if ((etaconst==2)||(etaconst==3)) {

                double t2=temp/zTc;
//        zeta=0.01162/sqrt(pow((t2-1.104),2)+ 0.0569777  )+  -0.1081/(t2*t2+23.7169);
//        tauRelax =5*(-0.0506*sTc/(temp*temp) + 10.453/(temp*sqrt(0.156658+pow((temp/sTc-1.131),2) )) );
                double min1=t2-1;
                if (t2>1.05) zeta=0.9*exp(-min1/0.025)+0.25*exp(-min1/0.13)+0.001;
                else if (t2<0.995) zeta=0.9*exp(min1/0.0025)+0.22*exp(min1/0.022)+0.03;
                else zeta=-13.77*t2*t2+27.55*t2-13.45;


                tauRelax =5.*zeta/(pow((1-EOS.cs2out(EOST())),2)*(EOSe()+EOSp()));    // single-argument version of cs2out

            }
            else if (etaconst==4) {
                double t2=temp/zTc;
                zeta=0.01162/sqrt(pow((t2-1.104),2)+ 0.0569777  )+  -0.1081/(t2*t2+23.7169);
                tauRelax =5*(-0.0506*sTc/(temp*temp) + 10.453/(temp*sqrt(0.156658+pow((temp/sTc-1.131),2) )) );


            }
            else
            {
                //bvf=5 and sig=2.1 for thin peak

//        zeta =bvf/(sig*sqrt(2*PI))*exp(-pow(temp-zTc,2)/(2.*sig*sig));
//
//        tauRelax =9*zeta/(EOSe()-3*EOSp());

                double t2=temp/zTc;
//        zeta=0.01162/sqrt(pow((t2-1.104),2)+ 0.0569777  )+  -0.1081/(t2*t2+23.7169);
//        tauRelax =5*(-0.0506*sTc/(temp*temp) + 10.453/(temp*sqrt(0.156658+pow((temp/sTc-1.131),2) )) );
                double min1=t2-1;
                if (t2>1.05) zeta=0.9*exp(-min1/0.025)+0.25*exp(-min1/0.13)+0.001;
                else if (t2<0.995) zeta=0.9*exp(min1/0.0025)+0.22*exp(min1/0.022)+0.03;
                else zeta=-13.77*t2*t2+27.55*t2-13.45;

                tauRelax =5.*zeta/(pow((1-EOS.cs2out(EOST())),2)*(EOSe()+EOSp()));    // single-argument version of cs2out

            }



            zeta*=EOSs();
            if (zeta<0.001) zeta=0.001;
            if (tauRelax <0.2) tauRelax=0.2 ;
        }




    }
}

template <int D>
void Particle<D>::sets(double tin2)
{

    gamma=gamcalc();
    shv.x[2][1]=shv.x[1][2];
    shv.x[0][1]=1./gamma*inner(u,colp1(1,shv));
    shv.x[0][2]=1./gamma*inner(u,colp1(2,shv));
    shv.x[1][0]=shv.x[0][1];
    shv.x[2][0]=shv.x[0][2];

    setvar();
    shv.x[0][0]=1./gamma/gamma*con(uu,pimin);
    shv33=(shv.x[0][0]-shv.x[1][1]-shv.x[2][2])/tin2;



}

template <int D>
void Particle<D>::setvar()
{
    mini(pimin,shv);
    uu=u*u;
}













//==================================================================
// Functions added by Christopher Plumberg

// THESE ARE THE NEW VERSIONS

template <int D>
double Particle<D>::EOST() { return SPH_cell.T; }
template <int D>
double Particle<D>::EOSmuB() { return SPH_cell.muB; }
template <int D>
double Particle<D>::EOSmuS() { return SPH_cell.muS; }
template <int D>
double Particle<D>::EOSmuQ() { return SPH_cell.muQ; }

template <int D>
double Particle<D>::EOSp() { return SPH_cell.p; }
template <int D>
double Particle<D>::EOSs() { return SPH_cell.s; }
template <int D>
double Particle<D>::EOSe() { return SPH_cell.e; }
template <int D>
double Particle<D>::EOSB() { return SPH_cell.B; }
template <int D>
double Particle<D>::EOSS() { return SPH_cell.S; }
template <int D>
double Particle<D>::EOSQ() { return SPH_cell.Q; }
template <int D>
double Particle<D>::EOSw() { return SPH_cell.w; }
template <int D>
double Particle<D>::EOSA() { return SPH_cell.A; }
template <int D>
double Particle<D>::EOSdwds() { return SPH_cell.dwds; }
template <int D>
double Particle<D>::EOSdwdB() { return SPH_cell.dwdB; }
template <int D>
double Particle<D>::EOSdwdS() { return SPH_cell.dwdS; }
template <int D>
double Particle<D>::EOSdwdQ() { return SPH_cell.dwdQ; }
template <int D>
double Particle<D>::EOSs_terms_T(double Tin) { return EOS.s_terms_T(Tin); }

template <int D>
double Particle<D>::EOSs_out(double e_In, double rhoB_In, double rhoS_In, double rhoQ_In)
{
	double sVal   = EOS.s_out( e_In, rhoB_In, rhoS_In, rhoQ_In );
	SPH_cell.T    = EOS.T();
	SPH_cell.muB  = EOS.muB();
	SPH_cell.muS  = EOS.muS();
	SPH_cell.muQ  = EOS.muQ();

	SPH_cell.p    = EOS.p();
	SPH_cell.s    = EOS.s();
	SPH_cell.B    = EOS.B();
	SPH_cell.S    = EOS.S();
	SPH_cell.Q    = EOS.Q();
	SPH_cell.e    = EOS.e();
	SPH_cell.cs2  = EOS.cs2();
	SPH_cell.w    = EOS.w();
	SPH_cell.A    = EOS.A();
	SPH_cell.dwds = EOS.dwds();
	SPH_cell.dwdB = EOS.dwdB();
	SPH_cell.dwdS = EOS.dwdS();
	SPH_cell.dwdQ = EOS.dwdQ();

	//cout << "Check entropy: " << sVal << " =?= " << EOS.s() << endl;

	return sVal;
}

template <int D>
double Particle<D>::EOSs_out(double e_In)
{
	double sVal = EOS.s_out( e_In, 0.0, 0.0, 0.0 );
	SPH_cell.T    = EOS.T();
	SPH_cell.muB  = EOS.muB();
	SPH_cell.muS  = EOS.muS();
	SPH_cell.muQ  = EOS.muQ();

	SPH_cell.p    = EOS.p();
	SPH_cell.s    = EOS.s();
	SPH_cell.B    = EOS.B();
	SPH_cell.S    = EOS.S();
	SPH_cell.Q    = EOS.Q();
	SPH_cell.e    = EOS.e();
	SPH_cell.cs2  = EOS.cs2();
	SPH_cell.w    = EOS.w();
	SPH_cell.A    = EOS.A();
	SPH_cell.dwds = EOS.dwds();
	SPH_cell.dwdB = EOS.dwdB();
	SPH_cell.dwdS = EOS.dwdS();
	SPH_cell.dwdQ = EOS.dwdQ();

	return sVal;
}

template <int D>
void Particle<D>::EOSupdate_s(double s_In, double rhoB_In, double rhoS_In, double rhoQ_In)
{
	bool update_s_success = EOS.update_s( s_In, rhoB_In, rhoS_In, rhoQ_In );
	SPH_cell.T    = EOS.T();
	SPH_cell.muB  = EOS.muB();
	SPH_cell.muS  = EOS.muS();
	SPH_cell.muQ  = EOS.muQ();

	SPH_cell.p    = EOS.p();
	SPH_cell.s    = EOS.s();
	SPH_cell.B    = EOS.B();
	SPH_cell.S    = EOS.S();
	SPH_cell.Q    = EOS.Q();
	SPH_cell.e    = EOS.e();
	SPH_cell.cs2  = EOS.cs2();
	SPH_cell.w    = EOS.w();
	SPH_cell.A    = EOS.A();
	SPH_cell.dwds = EOS.dwds();
	SPH_cell.dwdB = EOS.dwdB();
	SPH_cell.dwdS = EOS.dwdS();
	SPH_cell.dwdQ = EOS.dwdQ();

	return;
}

template <int D>
void Particle<D>::EOSupdate_s(double s_In)
{
	EOS.update_s( s_In, 0.0, 0.0, 0.0 );
	SPH_cell.T    = EOS.T();
	SPH_cell.muB  = EOS.muB();
	SPH_cell.muS  = EOS.muS();
	SPH_cell.muQ  = EOS.muQ();

	SPH_cell.p    = EOS.p();
	SPH_cell.s    = EOS.s();
	SPH_cell.B    = EOS.B();
	SPH_cell.S    = EOS.S();
	SPH_cell.Q    = EOS.Q();
	SPH_cell.e    = EOS.e();
	SPH_cell.cs2  = EOS.cs2();
	SPH_cell.w    = EOS.w();
	SPH_cell.A    = EOS.A();
	SPH_cell.dwds = EOS.dwds();
	SPH_cell.dwdB = EOS.dwdB();
	SPH_cell.dwdS = EOS.dwdS();
	SPH_cell.dwdQ = EOS.dwdQ();

	return;
}



/*
//THESE ARE THE OLD VERSIONS


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
double Particle<D>::EOSB()
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.B();
}


template <int D>
double Particle<D>::EOSS()
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.S();
}


template <int D>
double Particle<D>::EOSQ()
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	return EOS.Q();
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
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
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
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
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
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	bool update_s_success = EOS.update_s( s_In, rhoB_In, rhoS_In, rhoQ_In );
	particle_T = EOS.T();
	particle_muB = EOS.muB();
	particle_muS = EOS.muS();
	particle_muQ = EOS.muQ();

	return;
}




template <int D>
void Particle<D>::EOSupdate_s(double s_In)
{
	EOS.tbqs( particle_T, particle_muB, particle_muQ, particle_muS );
	EOS.update_s( s_In, 0.0, 0.0, 0.0 );
	particle_T = EOS.T();
	particle_muB = EOS.muB();
	particle_muS = EOS.muS();
	particle_muQ = EOS.muQ();

	return;
}
*/


#endif
