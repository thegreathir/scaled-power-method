#ifndef MATRIX_H_
#define MATRIX_H_

#include "utils.h"

#include <ios>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

namespace matrix
{

constexpr size_t kPrecision = 2;

template<class V>
class Matrix {
public:

    using DataType = std::vector<std::vector<V>>;

    DataType body;
    int m;
    int n;

    Matrix(const DataType& data): body(data), m(0), n(0) {
        m = data.size();
        if (m != 0) {
            n = data[0].size();
        }
    }

    Matrix(): m(0), n(0) {
    }

    template<class U = V>
    static Matrix<U> get_empty(size_t m, size_t n) {
        auto res = Matrix<U>();
        res.m = m;
        res.n = n;
        res.body.resize(res.m);
        for (size_t i = 0; i < res.m; ++i) {
            res.body[i].resize(res.n);
        }

        return res;
    }

    template<class U = V>
    Matrix<U> T() {
        auto res = get_empty(n, m);
        
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                res.body[j][i] = body[i][j];
            }
        }
        return res;
    }

    double norm2() {
        if (n != 1 && m != 1)
            throw std::invalid_argument("Norm2 not implemented for matrix, just vector");

        double res = 0;
        for (const auto& row : body) {
            for (const auto& cell : row) {
                res += std::pow(cell, 2);
            }
        }
        return std::sqrt(res);
    }

    template<class U = V>
    Matrix<U> operator * (const Matrix<U>& other) {
        if (n != other.m)
            throw std::invalid_argument("Can not multiply, dimension error");
        
        auto res = get_empty(m, other.n);

        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < other.n; ++j){
                U cell = 0;
                for (size_t k = 0; k < n; k++) {
                    cell += body[i][k] * other.body[k][j];
                }
                res.body[i][j] = cell;
            }
        }

        return res;

    }

    template<class U = V, class W>
    Matrix<U> operator * (W w) {
        auto res = get_empty(m, n);
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                res.body[i][j] = body[i][j] * w;
            }
        }

        return res;
    }

    template<class W>
    void operator *= (W w) {
        *this = *this * w;
    }

    template<class U = V, class W>
    void operator *= (const Matrix<U>& other)
    {
        *this = *this * other;
    }

};

}

namespace std
{

template<class T>
std::ostream& operator << (std::ostream& os, const matrix::Matrix<T>& matrix) {

    os <<  matrix.m << "x" << matrix.n << std::endl;
    if (matrix.m == 0 || matrix.n == 0) {
        os << "[]";
        return os;
    }
    
    auto max = std::vector<int>(matrix.n, -std::numeric_limits<int>::infinity());
    for (const auto& row : matrix.body) {
        size_t col = 0;
        for (const auto& cell : row) {
          max[col] =
              std::max(max[col], static_cast<int>(utils::get_fixed_length(cell, matrix::kPrecision)));
          ++col;
        }
    }

    size_t i = 0;
    for (const auto& row : matrix.body) {
        size_t col = 0;
        for (const auto& cell : row) {
            os << std::setprecision(matrix::kPrecision + utils::get_number_length(cell))
                    << std::setw(max[col] + 1) << cell;
            ++col;
        }

        if (i != matrix.body.size() - 1)
            os << std::endl;
        ++i;
    }

    return os;
}

}

#endif