#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
 
#include "kdtree.h"

int main() {
    try {
        test_wikipedia();
        std::cout << '\n';
        test_random(1000);
        std::cout << '\n';
        test_random(1000000);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
    return 0;
}
