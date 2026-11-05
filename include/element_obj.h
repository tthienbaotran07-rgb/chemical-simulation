#ifndef ELEMENT_OBJ_H
#define ELEMENT_OBJ_H

#include <string>

class Element {
private:
    int atomic_n;
    std::string sb;
    int outer_e;
    bool* index;

public:
    Element();
    explicit Element(int an);

    int get_atomic_n() const;
    std::string get_sb() const;
    int get_outer_e() const;
    bool* get_index();
    const bool* get_index() const;

    ~Element();
};

#endif