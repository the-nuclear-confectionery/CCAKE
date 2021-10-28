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
// #include "owner.h"
// #include "owner.cpp"


using std::cout;
using std::endl;
using std::string;
using std::vector;

animal::animal(){}

animal::animal(string Type, int Legs)
{
    // sound = Sound;
    type = Type;
    legs = Legs;
}

animal::~animal(){}



string animal::animal_sound()
{
    return sound;
}

int animal::count_legs()
{
    return legs;
}

string animal::animal_type()
{
    return type;
}

void animal::setOwnerPtr(string* ownerPtrIn)
{
    ownerPtr = ownerPtrIn;
}

// void animal::setNameByOwner(Ownder *owner)
// {
//     name = owner.getName();
// }

string animal::returnOwner()
{
    return *ownerPtr;
}