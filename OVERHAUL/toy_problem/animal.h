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
    animal();
    animal(string Sound, string Type, int Legs);
    ~animal();
};

public:
    string animal_sound();
    int count_legs();
    string animal_type();

    // setOwnerPtr(Owner *ownerPtrIn);
private:
    int legs,
    string sound,type,name;

    // void setNameByOwner(Owner &owner)
    