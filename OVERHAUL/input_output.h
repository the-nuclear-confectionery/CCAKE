#ifndef InputOutput_H
#define InputOutput_H

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
#include "particle.h"
#include "linklist.h"
#include "system_state.h"


class InputOutput
{
public:

InputOutput();
~InputOutput();


void load_settings_file( string path_to_settings_file ); // load setting
// paramters for simulation

void set_EoS_type(); // load in table for interpolation

void set_results_directory( string path_to_results_directory ); // sets up
// output directory, will update the outfile as time goes on

void read_in_initial_conditions(const SystemState &system); // talks to
// system state so that it can set initial system state

void print_system_state(); //at every time step, will write to output file


private:

string input_directory;
string output_directory;
void read_in_initial_conditions();

};



#endif




















