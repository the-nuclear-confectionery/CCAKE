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

#include "owner.h"

owner::owner(){}

owner::owner(string name)
{
    petName = name;
}

owner::~owner(){}

void getName()
{
    return petName;
}