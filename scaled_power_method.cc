#include "scaled_power_method.h"

namespace scaled_power_method {

using scaled_power_method::matrix::Matrix;

std::pair<double, Matrix<double>> CalculateScaledPowerMethod(const matrix::Matrix<double>& A,
                                                             const matrix::Matrix<double>& x0, size_t n) {
    if (x0.n != 1)
        throw std::invalid_argument("x0 is not vector");

    if (A.n != A.m)
        throw std::invalid_argument("A is not square");

    if (x0.m != A.m)
        throw std::invalid_argument("x0 and A should have same dimansion");

    auto normalizedX0 = x0 / x0.Norm2();
    auto xK = std::move(normalizedX0);

    for (size_t k = 0; k < n; ++k) {
        auto yK = A * xK;
        xK = yK / yK.Norm2();
    }

    auto lambdaK = (xK.T() * A * xK).body[0][0];
    return std::pair(lambdaK, xK);
}

}  // namespace scaled_power_method