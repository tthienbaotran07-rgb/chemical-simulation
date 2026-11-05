#include "element_obj.h"
#include "config.h"

Element::Element() : Element(1) {}

Element::Element(int an) : atomic_n(an), sb(get_symbol(an)), outer_e(chemotherapy(an)), index(new bool[2]) {
    if (an < 1 || an > 118) {
        atomic_n = 1;
        sb = get_symbol(atomic_n);
        outer_e = chemotherapy(atomic_n);
    }

    index[0] = yield(outer_e, atomic_n);
    index[1] = accept(outer_e, atomic_n);
}

int Element::get_atomic_n() const { return atomic_n; }

std::string Element::get_sb() const { return sb; }

int Element::get_outer_e() const { return outer_e; }

bool* Element::get_index() { return index; }

const bool* Element::get_index() const { return index; }

Element::~Element() {
    delete[] index;
    index = nullptr;
}