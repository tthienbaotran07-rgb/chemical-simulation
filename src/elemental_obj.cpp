#include "config.h"
class Element {
public:
    int atomic_n;
    Element(int an) : atomic_n(an) {}
    std::string sb = get_symbol(atomic_n);
    int chemo = chemotherapy(atomic_n);
};