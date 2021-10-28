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
#include "owner.h"
#include "owner.cpp"

animal::animal(){}

animal::animal(string Sound, String Type, int Legs)
{
    sound = Sound;
    type = Type;
    legs = Legs
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

string animal::type()
{
    return type;
}

void animal::setNameByOwner(Ownder *owner)
{
    name = owner.getName();
}