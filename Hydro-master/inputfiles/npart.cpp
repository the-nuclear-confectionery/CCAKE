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

int main (int argc, char *argv[]) //execute with ./a.out start_event# end_event#
{





	string infile,type;
	int ev1,ev2;
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
		
		
	}


	int evtot=ev2+1;
	
	vector<double> npart;
	npart.resize(evtot);
	
	for(int ev=ev1;ev<=ev2;ev++)
	{
	
	string name=infile+"/ic"+convertInt(ev)+".dat";
	ifstream input(name.c_str());
	if (!input.is_open())
 	{
 	cout << "Can't open " << name << endl;
 	exit(1);
 	}

	
	string line;
	getline(input,line);
	
	std::vector<std::string> pp = split(line, ' ');
	stringstream ss;
	ss << pp[3];
	ss >> npart[ev];
	input.close();
	}
	
	
	
	string name8;
	name8=infile+"/npart.dat";
	
	
	
	
	ofstream OUT2(name8.c_str());
	if (!OUT2.is_open())
 	{
 	cout << "Can't open " << name8 << endl;
 	exit(1);
 	}
	
	for(int ev=ev1;ev<=ev2;ev++){
		
		
		OUT2 <<  ev << " " << npart[ev] <<endl;
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



