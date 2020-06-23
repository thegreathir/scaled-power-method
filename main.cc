#include "matrix.h"
#include "scaled_power_method.h"

using scaled_power_method::matrix::Matrix;
using scaled_power_method::CalculateScaledPowerMethod;

int main() {
    auto A = Matrix<double>({
        {1, 1, 1},
        {-10, -1, 6},
        {10, -2, -9}
    });

    auto x0 = Matrix<double>({
        {1, 0, 0}
    }).T();

    std::cout << "A:" << std::endl;
    std::cout << A << std::endl;
    
    std::cout << "x0:" << std::endl;
    std::cout << x0 << std::endl;

    size_t n;

    std::cout << "Input N: ";
    std::cin >> n;

    auto [lambda, x] = CalculateScaledPowerMethod(A, x0, n);

    std::cout << "lambda: " << lambda << std::endl;
    std::cout << "x: " << std::endl << x << std::endl;

    return 0;
}