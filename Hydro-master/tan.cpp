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

struct var{

	double tau,T;

};

struct chr{

	double cen,npart;
	double Tmax,derv;

};

double tdis( vector <var> list);
string convertInt(int number);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);


int main (int argc, char *argv[]) //execute with ./a.out start_event# end_event#
{

	string infol,npfile;
	int ev1,ev2;
	if (!argv[1]||!argv[2])
	{
		cout << "Please include folder name and event numbers" << endl;
		exit(1);
	}
	else
	{
		infol=argv[1];
		
		stringstream s;
		s << argv[2];
		s >> ev1;
		
		if (argv[3]){
		stringstream s1;
		s1 << argv[3];
		s1 >> ev2;
		}
		else ev2=ev1;
		
		
		npfile=argv[4];
		
	}
	
	std::vector<std::string> minfol = split(infol, '/');

        int afol=minfol.size();
        string sfol;
        sfol=minfol[0];
        for (int j=1;j<(afol-1);j++){
          sfol=sfol+"/"+minfol[j];
        }
        cout << sfol << endl;
        string npfol=minfol[0]+"/"+minfol[1];
        string outfol=minfol[0]+"/"+minfol[1]+"/"+minfol[2]+"/"+minfol[3];
	
	
	
	
	
	vector<chr> all;
	for (int t=ev1;t<=ev2;t++){
	vector<double> templist,taulist;
	for (int ev=0;ev<=999;ev++){
	string name="outputfiles/"+infol+"/"+convertInt(t)+"/ev"+convertInt(ev)+".dat";
	ifstream input(name.c_str());
	if (!input.is_open())
 	{
 		chr tsub;
		tsub.Tmax=0;
		all.push_back(tsub);
 		continue;
 	}

	
	string line;
	vector< var > list;
	while (getline(input,line)){
	
	std::vector<double> y (2,0) ;
	std::vector<std::string> x = split(line, ' ');

	
	
	for(int j=0;j<2;j++)
	{
	stringstream s;
	s << x[j];
	s >> y[j];
	}
	
	var sub;
	sub.tau=y[0];
	sub.T=y[1];
	list.push_back(sub);
	}
	input.close();


	chr tsub;
	tsub.Tmax=list[0].T;
	tsub.derv=tdis(list);
	
	all.push_back(tsub);	
	}}
	
	
	//read in npart                                                                                                                                                  

        string nameb="df/out/"+npfol+"/npart.dat";
        ifstream inputb(nameb.c_str());
        if (!inputb.is_open())
        {
        cout << "Can't open " << nameb << endl;
        exit(1);
        }

        string line;

        int first=0;
        int pl=0;
	int cou=0;
	int tot=all.size();
	while (getline(inputb,line)&&(cou<tot)) {
	std::vector<std::string> x = split(line, ' ');
        if (x.size()<2) std::vector<std::string> x = split(line, '\t');
                                                                                                             
	double np;
	stringstream s;
        s << x[pl];
        s >> np;
        
        if (np==1&&first==0){
          pl=1;
          stringstream ss;
          ss << x[pl];
          ss >> np;          
	}

	first=1;                                                                                                                                
	all[cou].npart=np;
	cou++;
	}


        inputb.close();

	
	
    // read in centralities                                                                                                                                          
        vector<double> nph,npl;
        vector<double> cmax,cmin;
        string named="df/out/sorts/"+npfile+".dat";
        ifstream inputd(named.c_str());
        if (!inputd.is_open())
        {
        cout << "Can't open " << named << endl;
        exit(1);
        }
        while (getline(inputd,line)) {

	std::vector<double> y (2,0) ;
	std::vector<double> y2 (4,0) ;
        std::vector<std::string> x = split(line, ' ');
	if (x.size()<4) std::vector<std::string> x = split(line, '\t');



        for(int j=0;j<2;j++)
	{
        stringstream s;
	s << x[j];
	s >> y[j];
	}
        for(int j=2;j<4;j++)
        {
	stringstream s;
	s << x[j];
	s >> y2[j];
	}

	nph.push_back(y2[2]);
        npl.push_back(y2[3]);
	cmax.push_back(y[1]);
	cmin.push_back(y[0]);                                                                              
        }
        inputd.close();

	
	
	
	int cm=cmax.size();
	for(int c=0;c<cm;cm++)
	{
	

	string name2="outputfiles/"+outfol+"/cen"+convertInt(c)+".dat";
	ofstream OUT2(name2.c_str());
	if (!OUT2.is_open())
 	{
 	cout << "Can't open " << name2 << endl;
 	}
	
	
	for(int ev=0;ev<tot;ev++){
	
		if ((all[ev].npart<nph[c]) &&(all[ev].npart>npl[c])&&(all[ev].Tmax>0)) OUT2 <<  all[ev].Tmax << " " << all[ev].derv << endl;
	}
	
	OUT2.close();
	}

}

double tdis(vector <var> list){

	

	int max=list.size();
	
	vector<double> d;
	vector<double> tau;
	double davg=0;
	for (int s=1;s<max;s++){
	d.push_back((list[s].T-list[s-1].T)/(list[s].tau-list[s-1].tau));
	tau.push_back(0.5*(list[s].tau+list[s-1].tau));	
	davg+=d[s-1];
	}
	davg/=d.size();

	return davg;

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
