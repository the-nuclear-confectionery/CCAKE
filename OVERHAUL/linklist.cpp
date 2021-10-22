#ifndef _LINKLIST_H_
#define _LINKLIST_H_

#include "vector.h"
#include "matrix.h"
#include "particle.h"
#include "mathdef.h"
#include "random.h"
#include "eos.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
#include "Stopwatch.h"

// this contains functions for calling EoS directly
#include <lib.h>


LinkList::LinkList()
{
    range = 2; //number of boxes on the sides
    for (int i=0; i<2; i++) uni.x[i]=1.0;


}

void LinkList::destroy()
{
    delete [] lead;

}

void LinkList::setup(double it0, int ntot,double h, Particle<D> *_pin,double dtsave,int & numpart)
{
    t0=it0;
    _h=h;
    _n=ntot;
    _p=_pin;
    knorm=10/7./PI/(_h*_h);
    knorm2=knorm*0.25;
    kgrad=-10/7./PI/pow(_h,3)*3/4.;	//FIX MISSING MINUS SIGN!!!!!!  CONFIRM WITH JAKI
    kgrad2=10/7./PI/pow(_h,3)/_h;
    link=new int[_n];
    dael=new Vector<int,D>[_n];
    steps=100*(floor(tend-t0)+1);

//    ep1=new double [steps];
//    ep2=new double [steps];
//    ec1=new double [steps];
//    ec2=new double [steps];
//    tsave=new double [steps];

    dt=dtsave;
    number_part=numpart;
    avgetasig=0.;
}




void LinkList::endEV()
{
    delete [] link;
    delete [] dael;
}

LinkList::~LinkList(){}

int LinkList::triToSum(Vector<int,2>  dael, Vector<int,2>  size)
{
    return dael.x[0] + dael.x[1]*size.x[0];

}



void LinkList::bsqsvconservation()
{

//    conservation_entropy();
    bsqsvconservation_E();
    Etot=E+Ez;
    Eloss= (E0-Etot)/E0*100;
    rk2=0;
//   cout << Eloss << "% of Energy loss at time t=" << t << endl;
//    cout << (S0-S)/S0 << "% of Entropy loss at time t=" << t << endl;

}

void LinkList::conservation_entropy()
{

    S=0.;

    for (int i=0; i<_n; i++) {
        S+= _p[i].eta_sigma*_p[i].sigmaweight;
if (i==0)
		cout << "\t\t --> " << i << "   " << _p[i].eta_sigma << "   "
				<< _p[i].sigmaweight << "   " << S << endl;
    }

    if (first==1)
    {

        S0=S;
    }
}


// =============================================
// function to check conservation of B, S, and Q
void LinkList::conservation_BSQ()
{
    Btotal = 0.0;
    Stotal = 0.0;
    Qtotal = 0.0;

    for (int i=0; i<_n; i++)
	{
        //Btotal += _p[i].B;
        //Stotal += _p[i].S;
        //Qtotal += _p[i].Q;
        Btotal += _p[i].rhoB_sub*_p[i].rho_weight;
        Stotal += _p[i].rhoS_sub*_p[i].rho_weight;
        Qtotal += _p[i].rhoQ_sub*_p[i].rho_weight;
    }

    if (first==1)
    {
        Btotal0 = Btotal;
        Stotal0 = Stotal;
        Qtotal0 = Qtotal;
    }
	return;
}




void LinkList::bsqsvconservation_E()
{

    E=0.;
    for (int i=0; i<_n; i++)
    {
        E += ( _p[i].C* _p[i].g2 - _p[i].EOSp() - _p[i].bigPI + _p[i].shv.x[0][0] )
              / _p[i].sigma*_p[i].sigmaweight*t;
if (i==0)
cout << "E: " << i << "   " << t
		<< "   " << _p[i].EOST()
		<< "   " << _p[i].EOSe()
		<< "   " << _p[i].C
		<< "   " << _p[i].g2
		<< "   " << _p[i].EOSp()
		<< "   " << _p[i].bigPI
		<< "   " << _p[i].shv.x[0][0]
		<< "   " << _p[i].sigma
		<< "   " << _p[i].sigmaweight << endl;    }

    if (first==1)
    {
        first=0;
        E0=E;
    }
}




void LinkList::etas_set()
{

    for (int i=0; i<_n; i++) {

        _p[i].eta_sigma=_p[i].sigmaweight;
        _p[i].sigmaweight=1;
    }

}



void LinkList::bsqsv_set()
{

    for (int i=0; i<_n; i++) {
        double gg=_p[i].gamcalc();
        _p[i].g2=gg*gg;
        _p[i].shv33=0;
    }

}



void LinkList::bsqsvconservation_Ez()
{

    dEz=0.;

    double t2=t*t;
    for (int i=0; i<_n; i++) {
        dEz+=( _p[i].EOSp()+_p[i].bigPI+_p[i].shv33*t2)/_p[i].sigma*_p[i].sigmaweight;
if (false)
cout << "dEz: " << i << "   " << t
		<< "   " << _p[i].EOSp()
		<< "   " << _p[i].bigPI
		<< "   " << _p[i].shv33*t2
		<< "   " << _p[i].sigma
		<< "   " << _p[i].sigmaweight << endl;
    }

}




void LinkList::setshear()
{
    double t2=t*t;
    for (int i=0; i<_n; i++) _p[i].sets(t2);
}


void LinkList::initialize()
{
// check what happens with particle separates by itself?  Where in fortran code?
    //find system boundaries

    max=min=_p[0].r;

    for(int i=1; i<_n; i++)
    {
        for(int j=0; j<D; j++) {

            if(_p[i].r.x[j]>max.x[j]) max.x[j]=_p[i].r.x[j];

            if(_p[i].r.x[j]<min.x[j]) min.x[j]=_p[i].r.x[j];

        }
    }

    //evaluate system size

    //2*range puts extra boxes on sides of grid
    double sub=1./_h;

    size=sub*(max-min)+(2.*range+1.)*uni;

    //Size is the volume
    Size=1;

    Vector<double,D> dsub;

    // finds total volume of the system
    for(int i=0; i<D; i++) Size*=size.x[i];


    //dael: relates every particle with its linklist cube
    // also convert particle position to an integer

    for(int j=0; j<_n; j++) dael[j]=sub*(_p[j].r-min)+(1.*range)*uni;

    //lead: relates every linklist cube with one of the particles (leader) in it
    //link: links the leader particle of one cube with the others of the same cube
    // if only one particle in cube then it is the lead

    lead = new int[Size];

    for(int j=0; j<Size; j++) lead[j]=-1;


    for(int k=_n-1; k>=0; k--)
	{
        int tt=triToSum(dael[k],size); // need to understand still... seems like it needs another coordinate
        link[k]=lead[tt];
        lead[tt]=k;
    }

    return;
}





//if we include the SPH over rhoB, rhoS, rhoQ
void LinkList::bsqsvoptimization(int a, bool init_mode /*== false*/)
{
    _p[a].sigma = 0;
    _p[a].eta = 0;
    _p[a].rhoB_sub = 0.0;
    _p[a].rhoS_sub = 0.0;
    _p[a].rhoQ_sub = 0.0;
	int neighbor_count = 0;
    Vector<int,D> i;
    for(i.x[0]=-2; i.x[0]<=2; i.x[0]++)
        for(i.x[1]=-2; i.x[1]<=2; i.x[1]++)
        {
            int b = lead[triToSum(dael[a]+i, size)];
            while( b!=-1 )
            {
                double kern  = kernel(_p[a].r-_p[b].r);
				if (kern>0.0) neighbor_count++;
                _p[a].sigma += _p[b].sigmaweight*kern;
                _p[a].eta   += _p[b].sigmaweight*_p[b].eta_sigma*kern;
                _p[a].rhoB_sub  += _p[b].rho_weight*_p[b].rhoB_an*kern;    //confirm with Jaki
                _p[a].rhoS_sub  += _p[b].rho_weight*_p[b].rhoS_an*kern;    //confirm with Jaki
                _p[a].rhoQ_sub  += _p[b].rho_weight*_p[b].rhoQ_an*kern;    //confirm with Jaki

if (false)
std::cout << "bsqsvoptimization(SPH particle == " << a << " ): "
			<< b << "   " << _p[a].r
			<< "   " << _p[a].sigma
			<< "   " << _p[a].eta
			<< "   " << _p[b].r
			<< "   " << _p[b].sigmaweight
			<< "   " << _p[b].eta_sigma
			<< "   " << _p[b].rhoB_an
			<< "   " << _p[a].rhoB_sub
			<< "   " << _p[b].rhoS_an
			<< "   " << _p[a].rhoS_sub
			<< "   " << _p[b].rhoQ_an
			<< "   " << _p[a].rhoQ_sub
			<< "   " << kern << std::endl;

                b=link[b];
            }
        }

	//cout << "Check neighbor count: " << a << "   " << neighbor_count << endl;

	// reset total B, S, and Q charge of each SPH particle to reflect
	// smoothing from kernel function (ONLY ON FIRST TIME STEP)
	//cout << "-----------------------------------------------------------------" << endl;
	if ( init_mode )
	{
		//cout << "BEFORE: " << a << "   " << _p[a].B << "   " << _p[a].S << "   " << _p[a].Q << endl;
		//cout << _p[a].rho_weight << "   " << _p[a].rhoB_an << "   " << _p[a].rhoS_an << "   " << _p[a].rhoQ_an << endl;
		_p[a].B = _p[a].rhoB_sub * _p[a].rho_weight;
		_p[a].S = _p[a].rhoS_sub * _p[a].rho_weight;
		_p[a].Q = _p[a].rhoQ_sub * _p[a].rho_weight;
		//cout << "AFTER: " << a << "   " << _p[a].B << "   " << _p[a].S << "   " << _p[a].Q << endl;
		//cout << _p[a].rho_weight << "   " << _p[a].rhoB_sub << "   " << _p[a].rhoS_sub << "   " << _p[a].rhoQ_sub << endl;
		//cout << "-----------------------------------------------------------------" << endl;
	}

    return;
}




void LinkList::bsqsvoptimization2(int a,double tin,int & count)
{

    _p[a].gradP=0.;
    _p[a].gradBulk = 0.;
    //_p[a].gradrhoB = 0.;
    //_p[a].gradrhoS = 0.;
    //_p[a].gradrhoQ = 0.;
    _p[a].gradV = 0.;
    _p[a].gradshear=0.;
    _p[a].divshear=0.;

    Vector<int,D> i;

    if (_p[a].btrack!=-1) _p[a].btrack=0;

    double rdis=0;

    for(i.x[0]=-2; i.x[0]<=2; i.x[0]++)
    {
        for(i.x[1]=-2; i.x[1]<=2; i.x[1]++)
        {

            int b=lead[triToSum(dael[a]+i, size)];


            while(b!=-1 )
            {

                Vector<double,D> gradK=gradKernel(_p[a].r-_p[b].r, static_cast<bool>(a == 30 && b == 43));
                Vector<double,D> va=rowp1(0,_p[a].shv);
                Vector<double,D> vb=rowp1(0,_p[b].shv);
                Matrix<double,D,D> vminia,vminib;
                mini(vminia,_p[a].shv);
                mini(vminib,_p[b].shv);
                double sigsqra=1/(_p[a].sigma*_p[a].sigma);
                double sigsqrb=1/(_p[b].sigma*_p[b].sigma);
                Vector<double,D> sigsigK=_p[b].sigmaweight*_p[a].sigma*gradK;

/*if (a == 30 && b == 43)
	cout << "CHECK PARTICLE: " << a << "   " << b << "   " << t << "   "
		<< _p[a].r << "   " << _p[b].r << "   "
		<< _p[b].sigmaweight << "   " << _p[a].sigma << "   "
		<< gradK << "   " << sigsigK << endl;*/

                _p[a].gradP +=( sigsqrb*_p[b].EOSp()+ sigsqra*_p[a].EOSp() )*sigsigK;

                if (((Norm(_p[a].r-_p[b].r)/_h)<=2)&&(a!=b)) {
                    if (_p[a].btrack!=-1) _p[a].btrack++;
                    if (_p[a].btrack==1) rdis=Norm(_p[a].r-_p[b].r)/_h;
                }

                _p[a].gradBulk += ( _p[b].Bulk/_p[b].sigma/_p[b].gamma + _p[a].Bulk/_p[a].sigma/_p[a].gamma)/tin*sigsigK;
                //_p[a].gradrhoB += ( _p[b].rhoB/_p[b].sigma/_p[b].gamma + _p[a].rhoB/_p[a].sigma/_p[a].gamma)/tin*sigsigK;
                //_p[a].gradrhoS += ( _p[b].rhoS/_p[b].sigma/_p[b].gamma + _p[a].rhoS/_p[a].sigma/_p[a].gamma)/tin*sigsigK;
                //_p[a].gradrhoQ += ( _p[b].rhoQ/_p[b].sigma/_p[b].gamma + _p[a].rhoQ/_p[a].sigma/_p[a].gamma)/tin*sigsigK;
                _p[a].gradV+=(_p[b].sigmaweight/_p[a].sigma*( _p[b].v -  _p[a].v ))*gradK;



                _p[a].gradshear+=(inner(sigsigK,_p[a].v))*(sigsqrb*vb+sigsqra*va);
                _p[a].divshear+=(sigsqrb*(sigsigK*transpose(vminib))+sigsqra*(sigsigK*transpose(vminia)));



                if(isnan(_p[a].gradP.x[0])) {
                    cout << "gradP stopped working" << endl;
                    cout << t <<" "  << _p[a].gradP <<  " " << a << " " << b << endl;
                    cout << _p[b].sigmaweight << " " <<     _p[a].sigma << " " <<  _p[b].EOSp() << " " << endl;
                    cout <<   Size << " " <<  _p[b].EOSs() << " " << _p[a].EOSs() <<endl;

                    cout << _p[a].r << endl;
                    cout << _p[b].r << endl;
                    cout << kernel(_p[a].r-_p[b].r) << endl;

                }
                else if (isnan(_p[a].gradP.x[1])) {
                    cout << "1 " << gradPressure_weight(a,b) <<  " " << a << " " << b << endl;

                }
                else if (isnan(_p[a].gradP.x[2])) {
                    cout << "2 " << gradPressure_weight(a,b) <<  " " << a << " " << b << endl;

                }



                b=link[b];
            }
        }
    }

    if ((_p[a].btrack==1)&&((_p[a].EOST()*197.3)>=150)) {
        _p[a].frz2.t=tin;
    }
    else if ((_p[a].btrack==0)&&((_p[a].EOST()*197.3)>=150)&&(_p[a].Freeze<4)) {
        cout <<"Missed " << a << " "<< tin << "  " << _p[a].EOST()*197.3 << " " << rdis << " " << cfon <<  endl;
    }


}




#endif
