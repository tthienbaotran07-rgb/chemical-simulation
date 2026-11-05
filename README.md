# chemical-simulation

A small C++ program for balancing chemical equations with known reactants and
products.

## Build and run

Compile from the repository root with a C++17 compiler:

```sh
g++ -std=c++17 -Iinclude test.cpp src/reaction.cpp src/config.cpp -o chemical-simulation
```

Run the program and enter an unbalanced equation, for example:

```text
H2 + O2 -> H2O
```

It prints the balanced equation:

```text
2H2 + O2 -> 2H2O
```

Element symbols and nested parentheses or square brackets are supported, for
example `Al2(SO4)3 + Ca(OH)2 -> Al(OH)3 + CaSO4`. A leading coefficient is
accepted but recalculated. This version balances neutral formulas only; it does
not predict products, handle ionic charges, or parse hydrate notation. Equations
whose coefficients are not uniquely determined are reported instead of guessed.
