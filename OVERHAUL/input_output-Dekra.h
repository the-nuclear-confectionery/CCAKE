#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>


#include "mathdef.h"
#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "LinkList.h"
#include "eostables.h"
#include "eos.h"


//public

void readICs_iccing(string &firstry,  int &_Ntable3,Particle<2> *&_p,double factor,double const& sfcheck, int & numpart, eos EOS);//iccing (energy density+ conserved charges) 
void read_in_iccing_IC();
void read_in_settings_file();
void print_system_state();

//private



#endif




















