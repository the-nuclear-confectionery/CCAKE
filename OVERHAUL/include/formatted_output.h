#ifndef FORMATTED_OUTPUT_H
#define FORMATTED_OUTPUT_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using std::string;
using std::vector;

namespace formatted_output
{
  void announce(string message);
  void report(string message);
  void update(string message);
  void comment(string message);
  void detail(string message);
//  void exposit(string message);
};

#endif