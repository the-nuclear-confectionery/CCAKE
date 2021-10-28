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
 vector<string> types = {"dog","cat","komono dragom","duck"};
 vector<int> legs = {4,4,4,2};
 vector<animal> vecOfAnimals;
 string OWNER = "Jimbo";
 for (int i = 0; i<types.size();i++)
 {
    vecOfAnimals.push_back(animal(types[i],legs[i]));
 }

 vector<animal> tempAnimals;
  for (int i = 0; i<types.size();i++)
 {
     vecOfAnimals[i].setOwnerPtr(*OWNER);
     cout << vecOfAnimals[i] << endl;
     if(vecOfAnimals[i].count_legs()==4)
     {
         tempAnimals.push_back(vecOfAnimals[i]);
     }
 }

 vecOfAnimals = tempAnimals;

   for (int i = 0; i<vecOfAnimals.size();i++)
 {
     cout << vecOfAnimals[i] << endl;
 }




    return 0;
}