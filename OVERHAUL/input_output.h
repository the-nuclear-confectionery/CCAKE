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
#include "system.h"


class input_output
{

public:

  input_output();
  ~input_output();

  void load_settings_file( string path_to_settings_file );
  void set_results_directory( string path_to_results_directory );
  void read_in_initial_conditions();
  void print_results();

private:

  string input_directory
  string output_directory;

}



#endif




















