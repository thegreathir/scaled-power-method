#ifndef UTILS_H_
#define UTILS_H_

#include <iomanip>
#include <iostream>
#include <sstream>

namespace scaled_power_method {

namespace utils {

template <class T>
inline size_t GetFixedLength(T t, size_t precision) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << t;
    return stream.str().length();
}

}  // namespace utils

}  // namespace scaled_power_method
#endif