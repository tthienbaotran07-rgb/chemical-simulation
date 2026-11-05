#include "reaction.h"

#include "config.h"

#include <algorithm>
#include <cctype>
#include <numeric>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>

namespace {

using AtomCounts = std::map<std::string, long long>;

std::string trim(const std::string& value) {
    const std::string whitespace = " \t\r\n";
    const std::size_t first = value.find_first_not_of(whitespace);
    if (first == std::string::npos) {
        return "";
    }
    const std::size_t last = value.find_last_not_of(whitespace);
    return value.substr(first, last - first + 1);
}

class FormulaParser {
public:
    explicit FormulaParser(const std::string& formula) : formula_(formula) {}

    AtomCounts parse() {
        AtomCounts counts = parse_group('\0');
        if (position_ != formula_.size()) {
            throw std::invalid_argument("Unexpected character in formula: " +
                                        formula_.substr(position_));
        }
        if (counts.empty()) {
            throw std::invalid_argument("A chemical formula cannot be empty.");
        }
        return counts;
    }

private:
    AtomCounts parse_group(char closing) {
        AtomCounts counts;
        while (position_ < formula_.size()) {
            const char current = formula_[position_];
            if (current == ')' || current == ']') {
                if (closing == '\0' || current != closing) {
                    throw std::invalid_argument("Mismatched parentheses in formula: " + formula_);
                }
                ++position_;
                return counts;
            }

            if (current == '(' || current == '[') {
                ++position_;
                const char expected = current == '(' ? ')' : ']';
                AtomCounts group = parse_group(expected);
                if (group.empty()) {
                    throw std::invalid_argument("Empty group in formula: " + formula_);
                }
                const long long multiplier = parse_number();
                for (const auto& atom : group) {
                    counts[atom.first] += atom.second * multiplier;
                }
                continue;
            }

            if (current < 'A' || current > 'Z') {
                throw std::invalid_argument("Expected an element symbol in formula: " + formula_);
            }

            std::string element(1, current);
            ++position_;
            while (position_ < formula_.size() &&
                   formula_[position_] >= 'a' && formula_[position_] <= 'z') {
                element += formula_[position_++];
            }
            if (get_atomic(element) == 0) {
                throw std::invalid_argument("Unknown element symbol: " + element);
            }
            counts[element] += parse_number();
        }

        if (closing != '\0') {
            throw std::invalid_argument("Unclosed parentheses in formula: " + formula_);
        }
        return counts;
    }

    long long parse_number() {
        if (position_ >= formula_.size() ||
            formula_[position_] < '0' || formula_[position_] > '9') {
            return 1;
        }

        long long value = 0;
        while (position_ < formula_.size() &&
               formula_[position_] >= '0' && formula_[position_] <= '9') {
            const int digit = formula_[position_++] - '0';
            if (value > (std::numeric_limits<long long>::max() - digit) / 10) {
                throw std::invalid_argument("Number is too large in formula: " + formula_);
            }
            value = value * 10 + digit;
        }
        if (value == 0) {
            throw std::invalid_argument("Atom counts and multipliers must be positive.");
        }
        return value;
    }

    const std::string& formula_;
    std::size_t position_ = 0;
};

struct Species {
    std::string formula;
    AtomCounts atoms;
};

Species parse_species(std::string text) {
    text = trim(text);
    if (text.empty()) {
        throw std::invalid_argument("Each side of the equation must contain a substance.");
    }

    std::size_t position = 0;
    while (position < text.size() && text[position] >= '0' && text[position] <= '9') {
        ++position;
    }
    if (position > 0) {
        const std::string coefficient = text.substr(0, position);
        if (coefficient.find_first_not_of('0') == std::string::npos) {
            throw std::invalid_argument("A substance coefficient must be positive.");
        }
        while (position < text.size() &&
               std::isspace(static_cast<unsigned char>(text[position]))) {
            ++position;
        }
        text = text.substr(position);
    }

    FormulaParser parser(text);
    return Species{text, parser.parse()};
}

std::vector<Species> parse_side(const std::string& side) {
    std::vector<Species> species;
    std::size_t start = 0;
    while (true) {
        const std::size_t separator = side.find('+', start);
        species.push_back(parse_species(side.substr(
            start, separator == std::string::npos ? std::string::npos : separator - start)));
        if (separator == std::string::npos) {
            break;
        }
        start = separator + 1;
    }
    return species;
}

class Fraction {
public:
    Fraction(long long numerator = 0, long long denominator = 1)
        : numerator_(numerator), denominator_(denominator) {
        normalize();
    }

    long long numerator() const { return numerator_; }
    long long denominator() const { return denominator_; }
    bool is_zero() const { return numerator_ == 0; }

    Fraction operator-() const { return Fraction(-numerator_, denominator_); }

    Fraction& operator+=(const Fraction& other) {
        numerator_ = numerator_ * other.denominator_ + other.numerator_ * denominator_;
        denominator_ *= other.denominator_;
        normalize();
        return *this;
    }

    Fraction& operator-=(const Fraction& other) {
        return *this += -other;
    }

    Fraction& operator*=(const Fraction& other) {
        numerator_ *= other.numerator_;
        denominator_ *= other.denominator_;
        normalize();
        return *this;
    }

    Fraction& operator/=(const Fraction& other) {
        if (other.is_zero()) {
            throw std::logic_error("Cannot divide by zero while balancing equation.");
        }
        numerator_ *= other.denominator_;
        denominator_ *= other.numerator_;
        normalize();
        return *this;
    }

private:
    void normalize() {
        if (denominator_ == 0) {
            throw std::logic_error("Invalid fraction while balancing equation.");
        }
        if (denominator_ < 0) {
            numerator_ = -numerator_;
            denominator_ = -denominator_;
        }
        const long long divisor = std::gcd(numerator_, denominator_);
        numerator_ /= divisor;
        denominator_ /= divisor;
    }

    long long numerator_;
    long long denominator_;
};

Fraction operator*(Fraction left, const Fraction& right) {
    left *= right;
    return left;
}

Fraction operator/(Fraction left, const Fraction& right) {
    left /= right;
    return left;
}

std::vector<long long> solve_coefficients(
    const std::vector<Species>& reactants,
    const std::vector<Species>& products) {
    std::vector<std::string> elements;
    for (const Species& species : reactants) {
        for (const auto& atom : species.atoms) {
            elements.push_back(atom.first);
        }
    }
    for (const Species& species : products) {
        for (const auto& atom : species.atoms) {
            elements.push_back(atom.first);
        }
    }
    std::sort(elements.begin(), elements.end());
    elements.erase(std::unique(elements.begin(), elements.end()), elements.end());

    const std::size_t columns = reactants.size() + products.size();
    std::vector<std::vector<Fraction>> matrix(
        elements.size(), std::vector<Fraction>(columns));
    for (std::size_t row = 0; row < elements.size(); ++row) {
        for (std::size_t col = 0; col < reactants.size(); ++col) {
            auto atom = reactants[col].atoms.find(elements[row]);
            if (atom != reactants[col].atoms.end()) {
                matrix[row][col] = Fraction(atom->second);
            }
        }
        for (std::size_t col = 0; col < products.size(); ++col) {
            auto atom = products[col].atoms.find(elements[row]);
            if (atom != products[col].atoms.end()) {
                matrix[row][reactants.size() + col] = Fraction(-atom->second);
            }
        }
    }

    std::vector<std::size_t> pivot_columns;
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < columns && pivot_row < matrix.size(); ++col) {
        std::size_t selected = pivot_row;
        while (selected < matrix.size() && matrix[selected][col].is_zero()) {
            ++selected;
        }
        if (selected == matrix.size()) {
            continue;
        }

        std::swap(matrix[pivot_row], matrix[selected]);
        const Fraction pivot = matrix[pivot_row][col];
        for (std::size_t j = col; j < columns; ++j) {
            matrix[pivot_row][j] = matrix[pivot_row][j] / pivot;
        }
        for (std::size_t row = 0; row < matrix.size(); ++row) {
            if (row == pivot_row || matrix[row][col].is_zero()) {
                continue;
            }
            const Fraction factor = matrix[row][col];
            for (std::size_t j = col; j < columns; ++j) {
                matrix[row][j] -= factor * matrix[pivot_row][j];
            }
        }
        pivot_columns.push_back(col);
        ++pivot_row;
    }

    if (columns - pivot_columns.size() != 1) {
        throw std::invalid_argument(
            "The equation does not have a unique balancing ratio. "
            "Check the substances on both sides.");
    }

    std::size_t free_column = 0;
    while (std::find(pivot_columns.begin(), pivot_columns.end(), free_column) !=
           pivot_columns.end()) {
        ++free_column;
    }

    std::vector<Fraction> solution(columns);
    solution[free_column] = Fraction(1);
    for (std::size_t row = 0; row < pivot_columns.size(); ++row) {
        solution[pivot_columns[row]] = -matrix[row][free_column];
    }

    long long common_denominator = 1;
    for (const Fraction& value : solution) {
        common_denominator =
            std::lcm(common_denominator, value.denominator());
    }

    std::vector<long long> coefficients(columns);
    for (std::size_t i = 0; i < columns; ++i) {
        coefficients[i] = solution[i].numerator() *
                          (common_denominator / solution[i].denominator());
    }

    if (coefficients.front() < 0) {
        for (long long& coefficient : coefficients) {
            coefficient = -coefficient;
        }
    }
    if (std::any_of(coefficients.begin(), coefficients.end(),
                    [](long long value) { return value <= 0; })) {
        throw std::invalid_argument(
            "No positive balancing coefficients exist for these substances.");
    }

    long long divisor = coefficients.front();
    for (long long coefficient : coefficients) {
        divisor = std::gcd(divisor, coefficient);
    }
    for (long long& coefficient : coefficients) {
        coefficient /= divisor;
    }
    return coefficients;
}

void append_side(
    std::ostringstream& output,
    const std::vector<Species>& species,
    const std::vector<long long>& coefficients,
    std::size_t offset) {
    for (std::size_t i = 0; i < species.size(); ++i) {
        if (i > 0) {
            output << " + ";
        }
        if (coefficients[offset + i] != 1) {
            output << coefficients[offset + i];
        }
        output << species[i].formula;
    }
}

} // namespace

std::string balance_equation(const std::string& equation) {
    const std::size_t arrow = equation.find("->");
    if (arrow == std::string::npos ||
        equation.find("->", arrow + 2) != std::string::npos) {
        throw std::invalid_argument("Enter one reaction arrow using '->'.");
    }

    const std::vector<Species> reactants = parse_side(equation.substr(0, arrow));
    const std::vector<Species> products = parse_side(equation.substr(arrow + 2));
    const std::vector<long long> coefficients =
        solve_coefficients(reactants, products);

    std::ostringstream result;
    append_side(result, reactants, coefficients, 0);
    result << " -> ";
    append_side(result, products, coefficients, reactants.size());
    return result.str();
}
