#include "matrix.h"

using matrix::Matrix;

int main() {

    auto matrix1 = Matrix<float>({
        {1, 3, 12, 12},
        {2, 5, 12, 12},
        {8, 10, 12, 12},
        {-10, 5, 12, 12}
    });

    auto matrix2 = Matrix<float>({
        {1, 3, 21, 17},
        {19, -23, 22, 17},
        {-5, -10, 23, 34},
        {1, -100, -1, 22},
    });

    matrix2 *= matrix1;
    matrix2 *= 2;

    std::cout << matrix2 << std::endl;

    return 0;
}