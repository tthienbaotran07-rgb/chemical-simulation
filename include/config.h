#ifndef CONFIG_H
#define CONFIG_H
#include <string>
extern const std::string symbol[118];
std::string get_symbol(int z);
int get_atomic(std::string x);
int chemotherapy(int z);
bool yield(int outer, int z);
bool accept(int outer, int z);
#endif