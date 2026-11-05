#include <iostream>
#include <stdexcept>
#include <string>
#include "reaction.h"

int main() {
    std::string equation;
    std::cout << "Enter an equation (for example: H2 + O2 -> H2O): ";
    std::getline(std::cin, equation);

    try {
        std::cout << balance_equation(equation) << '\n';
    } catch (const std::invalid_argument& error) {
        std::cerr << "Cannot balance equation: " << error.what() << '\n';
        return 1;
    }
    return 0;
}