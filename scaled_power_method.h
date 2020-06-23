#ifndef SCALED_POWER_METHOD_H_
#define SCALED_POWER_METHOD_H_

#include "matrix.h"

namespace scaled_power_method {

std::pair<double, matrix::Matrix<double>> CalculateScaledPowerMethod(const matrix::Matrix<double>& A,
                                                                     const matrix::Matrix<double>& x0, size_t n);

}

#endif
