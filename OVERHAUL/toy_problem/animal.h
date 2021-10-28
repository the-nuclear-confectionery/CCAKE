#ifndef animal_h
#define animal_h

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

using std::cout;
using std::endl;
using std::string;
using std::vector;

// #include "owner.h"
// #include "owner.cpp"

class animal
{
public:
    animal();
    animal(string Sound, string Type, int Legs);
    ~animal();


    string animal_sound();
    int count_legs();
    string animal_type();

    void setOwnerPtr(string *ownerPtrIn);
    string returnOwner();
private:
    int legs;
    string sound,type,name;
    string *ownerPtr;

    // void setNameByOwner(Owner &owner);
    
};

#endif