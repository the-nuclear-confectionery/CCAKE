#include <iostream>

#include "unit_tests.h"

using namespace std;

int main()
{
  Vector<double, 2> a, b;
  a(0) = 1; a(1) = 2;
  b(0) = 1; b(1) = 2;

  if ( a == b )
    cout << "They are equal!" << endl;
  else
    cout << "They are not equal!" << endl;

  return 0;
}