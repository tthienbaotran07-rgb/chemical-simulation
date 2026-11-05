#include <string>
#include "config.h"
//dài vãi
const std::string symbol[118] = {
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};
std::string get_symbol(int z) {
    if (z < 1 || z > 118) {
        return "";
    }
    return symbol[z - 1];
}
int get_atomic(std::string x) {
    for (int i = 0; i < 118; i++) {
        if (symbol[i] == x) return (i + 1);
    }
    return 0;
}
int chemotherapy(int z) {
    if (z <= 0 || z > 118) {
        return 0;
    }
    const int subshell_n[] = {
        1, 2, 2, 3, 3, 4, 3, 4, 5, 4, 5, 6,
        4, 5, 6, 7, 5, 6, 7
    };
    const int subshell_capacity[] = {
        2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2,
        14, 10, 6, 2, 14, 10, 6
    };

    int remaining = z;
    int outer_shell = 0;
    int outer_electrons = 0;
    for (int i = 0; i < 19 && remaining > 0; i++) {
        int filled = remaining < subshell_capacity[i]
            ? remaining
            : subshell_capacity[i];
        if (subshell_n[i] > outer_shell) {
            outer_shell = subshell_n[i];
            outer_electrons = filled;
        } else if (subshell_n[i] == outer_shell) {
            outer_electrons += filled;
        }
        remaining -= filled;
    }
    return outer_electrons;
}
bool yield(int outer, int z) {
    if (outer == 1 || outer == 2 || outer == 3 || outer == 4) {
        return !(z == 1 || z == 2);
    }
    return false;
}
bool accept(int outer, int z) {
    if (z == 1) return true;
    if (outer == 5 || outer == 6 || outer == 7) {
        return !(z == 83 || z == 115 || z == 116);
    }
    return false;
}