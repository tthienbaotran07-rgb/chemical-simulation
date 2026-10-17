#include <iostream>
#include <string>
#include "config.h"
//ngon lành
int main() {
    int z;std::cin >> z;
    std::string s;std::cin >> s;
    std::cout << get_symbol(z) << "\n";
    std::cout << get_atomic_number(s);
    return 0;
}