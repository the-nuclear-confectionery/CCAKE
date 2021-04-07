#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <string.h>
#include <vector>


using namespace std;

string convertInt(int number);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

double ecc(double  & psi,vector <double> rx,vector <double> ry,vector <double> r2, vector <double> phi,vector <double> e, int m,int n,double & rout);

extern double g2;
double g2;

int main (int argc, char *argv[]) //execute with ./a.out start_event# end_event#
{

	string cgc="cgc";
	string nex="nexus";
	string glb="glauber";
	string glas="glasma";



	string infile,type;
	int ev1,ev2;
	double grid;
	if (!argv[1]||!argv[2])
	{
		cout << "Please include folder name and event numbers" << endl;
		exit(1);
	}
	else
	{
		infile=argv[1];
		
		stringstream s;
		s << argv[2];
		s >> ev1;
		
		if (argv[3]){
		stringstream s1;
		s1 << argv[3];
		s1 >> ev2;
		}
		else ev2=ev1;
		
		if (argv[4]){
		stringstream s3;
		s3 << argv[4];
		s3 >> grid;
		}
		
		if (argv[5]) type=argv[5];
		
	}
	g2=grid*grid;
	cout << "grid=" << g2 << " type=" << type << endl;
	

	int evtot=ev2+1;
	
	vector< double > ec2,rad,psi2,ec3,psi3,ec4,psi4,ec5,psi5,ec6,psi6,ec15,psi15,ec33,psi33,ec34,psi34,ec35,psi35,ec36,psi36,ec11,psi11,ec12,psi12,ec44,psi44,ec13,psi13;
	ec2.resize(evtot);
	rad.resize(evtot);
	psi2.resize(evtot);
	ec3.resize(evtot);
	psi3.resize(evtot);
	ec4.resize(evtot);
	psi4.resize(evtot);
	ec5.resize(evtot);
	psi5.resize(evtot);
	ec6.resize(evtot);
	psi6.resize(evtot);
	ec15.resize(evtot);
	psi15.resize(evtot);
	ec33.resize(evtot);
	psi33.resize(evtot);
	ec34.resize(evtot);
	psi34.resize(evtot);
	ec35.resize(evtot);
	psi35.resize(evtot);
	ec36.resize(evtot);
	psi36.resize(evtot);
	ec11.resize(evtot);
	psi11.resize(evtot);
	ec12.resize(evtot);
	psi12.resize(evtot);
	ec44.resize(evtot);
	psi44.resize(evtot);
	ec13.resize(evtot);
	psi13.resize(evtot);
	
	
	vector<double> npart;
	npart.resize(evtot);
	
	for(int ev=ev1;ev<=ev2;ev++)
	{
	
	string name=infile+"/ic"+convertInt(ev)+".dat";
//	cout << name << endl;
	ifstream input(name.c_str());
	if (!input.is_open())
 	{
 	cout << "Can't open " << name << endl;
 	exit(1);
 	}

	
	string line;
	
	if ((type==cgc)||(type==glb)) {
	getline(input,line);
	
//	std::vector<std::string> pp = split(line, ' ');
//	stringstream ss;
//	ss << pp[3];
//	ss >> npart[ev];
	
	}
	
	int keep=0;
	
	vector< double > rx,ry,e,r2,phi;
	while (getline(input,line)){
	
	std::vector<double> y (3,0) ;
	std::vector<std::string> x = split(line, ' ');
	
	if (x.size()<3) { 
	x.clear();
	std::vector<std::string> g = split(line, '\t');
	x=g;
	}
	
	
	for(int j=0;j<3;j++)
	{
	stringstream s;
	s << x[j];
	s >> y[j];
	}
	
	e.push_back(y[2]);
	rx.push_back(y[0]);
	ry.push_back(y[1]);
//	r2.push_back(y[0]*y[0]+y[1]*y[1]);
//	phi.push_back(atan2(y[1],y[0]) );
	
	
	
	
	
	
	}
	input.close();
	
	int max2=e.size();
	r2.resize(max2);
	phi.resize(max2);
	
	//cout << e.size() << endl;
	
	
	int two=2,three=3;
	double psisb,sub;
	
	double rads;
	
	sub=ecc(psisb,rx,ry,r2,phi,e,two,two,rads);
	
	ec2[ev]=sub;
	psi2[ev]=psisb;
	rad[ev]=rads;
	
	
	ec3[ev]=ecc(psisb,rx,ry,r2,phi,e,two,three,rads);
	psi3[ev]=psisb;
	
	ec4[ev]=ecc(psisb,rx,ry,r2,phi,e,2,4,rads);
	psi4[ev]=psisb;

	ec5[ev]=ecc(psisb,rx,ry,r2,phi,e,2,5,rads);
	psi5[ev]=psisb;
	
	ec6[ev]=ecc(psisb,rx,ry,r2,phi,e,2,6,rads);
	psi6[ev]=psisb;
	
	ec15[ev]=ecc(psisb,rx,ry,r2,phi,e,1,5,rads);
	psi15[ev]=psisb;
	
	ec33[ev]=ecc(psisb,rx,ry,r2,phi,e,3,3,rads);
	psi33[ev]=psisb;
	
	
	ec34[ev]=ecc(psisb,rx,ry,r2,phi,e,3,4,rads);
	psi34[ev]=psisb;
	
	ec35[ev]=ecc(psisb,rx,ry,r2,phi,e,3,5,rads);
	psi35[ev]=psisb;
	
	ec36[ev]=ecc(psisb,rx,ry,r2,phi,e,3,6,rads);
	psi36[ev]=psisb;
	
	
	ec12[ev]=ecc(psisb,rx,ry,r2,phi,e,1,2,rads);
	psi12[ev]=psisb;
	
	ec44[ev]=ecc(psisb,rx,ry,r2,phi,e,4,4,rads);
	psi44[ev]=psisb;
	
	ec13[ev]=ecc(psisb,rx,ry,r2,phi,e,1,3,rads);
	psi13[ev]=psisb;
	 
	 
	 
	} // finish running over all events
	
	
		
	
	
	
	
	
	
	
	
	
	string name8;
	name8=infile+"/eccCM.dat";
	
	
	
	
	ofstream OUT2(name8.c_str());
	if (!OUT2.is_open())
 	{
 	cout << "Can't open " << name8 << endl;
 	exit(1);
 	}
	
	for(int ev=ev1;ev<=ev2;ev++){
		
		
		OUT2 <<  ec2[ev] << " " << psi2[ev] <<  " " <<   ec3[ev] << " " << psi3[ev] <<  " " <<  ec4[ev] << " " << psi4[ev] <<  " " <<  ec5[ev] << " " << psi5[ev] <<  " " <<  ec6[ev] << " " << psi6[ev] <<  " " <<  ec15[ev] << " " << psi15[ev] <<  " " <<  ec33[ev] << " " << psi33[ev] <<  " " <<  ec34[ev] << " " << psi34[ev] <<  " " <<  ec35[ev] << " " << psi35[ev] <<  " " <<  ec36[ev] << " " << psi36[ev] <<    " " <<   ec12[ev] << " " << psi12[ev] <<  " " << ec44[ev] << " " << psi44[ev] << " " << ec13[ev] << " " << psi13[ev] <<   " " << rad[ev] <<endl;
	}
	
	OUT2.close();
	
	




}

string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {    	
        if (!item.empty()) elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


double ecc(double  & psi,vector <double> rx,vector <double> ry,vector <double> r2, vector <double> phi,vector <double> e, int m,int n,double & rout){

	

	int max=e.size();
	
	double xcm=0,ycm=0,etot=0;
	for (int s=0;s<max;s++){
	xcm+=rx[s]*e[s];
	ycm+=ry[s]*e[s];
	etot+=e[s];
	}
	xcm/=etot;
	ycm/=etot;
	
	
	double psit=0,psib=0,rb=0;
	for (int s=0;s<max;s++){
	double xsub=(rx[s]-xcm);
	double ysub=(ry[s]-ycm);
	r2[s]=xsub*xsub+ysub*ysub;
	phi[s]=atan2(ysub,xsub);
	double rv=e[s]*pow(r2[s],(m/2.));
	psit+=rv*sin(1.0*n*phi[s]);
	psib+=rv*cos(1.0*n*phi[s]);
	
	rb+=rv;
	}
	
	
	psit/=max;
	psib/=max;
	
	
	
	
	psi=1./(1.0*n)*atan2(psit,psib);
	//if (n==3&&m==3) cout << 1./(1.0*n) << " " << atan2(psit,psib) << " " <<  psi << endl;
	
	
	double ec=0;
	for (int s=0;s<max;s++) ec+=e[s]*pow(r2[s],m/2.)*cos(n*(phi[s]-psi));

	ec/=rb;

	rout=rb/etot;
	
	return ec;

}

