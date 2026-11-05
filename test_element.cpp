#include <iostream>

#include "element_obj.h"

int main() {
    int z;
    std::cout << "Enter atomic number (1-118): ";
    if (!(std::cin >> z)) {
        std::cerr << "Invalid input.\n";
        return 1;
    }

    Element element(z);
    const bool* index = element.get_index();
    std::cout << "Atomic number: " << element.get_atomic_n() << '\n'
              << "Symbol: " << (element.get_sb().empty() ? "(invalid)" : element.get_sb()) << '\n'
              << "Outer-shell electrons: " << element.get_outer_e() << '\n'
              << std::boolalpha
              << "Can yield: " << index[0] << '\n'
              << "Can accept: " << index[1] << '\n';
    return 0;
}