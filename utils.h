#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <sstream>
#include <iomanip>

namespace utils {

size_t get_number_length(long long int n) {
    n = std::abs(n);
    size_t number_of_digits = 0;
    do {
        ++number_of_digits; 
        n /= 10;
    } while (n);
    return number_of_digits;
}

template <class T> size_t get_fixed_length(T t, size_t precision) {
    std::stringstream stream;
    stream << std::setprecision(precision + get_number_length(t)) << t;
    return stream.str().length();
}

}

#endif