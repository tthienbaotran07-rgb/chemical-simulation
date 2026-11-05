#include "reaction.h"

#include <cassert>
#include <stdexcept>
#include <string>

int main() {
    assert(balance_equation("H2 + O2 -> H2O") == "2H2 + O2 -> 2H2O");
    assert(balance_equation("Fe + O2 -> Fe2O3") == "4Fe + 3O2 -> 2Fe2O3");
    assert(balance_equation(
               "Al2(SO4)3 + Ca(OH)2 -> Al(OH)3 + CaSO4") ==
           "Al2(SO4)3 + 3Ca(OH)2 -> 2Al(OH)3 + 3CaSO4");
    assert(balance_equation("2 H2 + O2 -> H2O") == "2H2 + O2 -> 2H2O");

    bool rejected_unknown_element = false;
    try {
        balance_equation("Xx + O2 -> XxO2");
    } catch (const std::invalid_argument&) {
        rejected_unknown_element = true;
    }
    assert(rejected_unknown_element);

    bool rejected_ambiguous_equation = false;
    try {
        balance_equation("C + O2 -> CO + CO2");
    } catch (const std::invalid_argument&) {
        rejected_ambiguous_equation = true;
    }
    assert(rejected_ambiguous_equation);
}
