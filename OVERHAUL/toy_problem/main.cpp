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

#include "animal.h"
#include "animal.cpp"


using std::cout;
using std::endl;
using std::string;
using std::vector;

int main (int argc, char *argv[])
{
 vector<string> types = ({"dog","cat","komono dragom","duck"});
 vector<int> legs = ({4,4,4,2});
 vector<animals> vecOfAnimals;
 for (int i = 0; i<types.size;i++)
 {
    vecOfAnimals.push_back(animal(types[i],legs[i]))
 }

 




    return 0;
}