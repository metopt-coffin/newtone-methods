#include "util/VectorOps.h"

#include <algorithm>
#include <functional>
#include <numeric>

namespace util {

VectorT negate(VectorT vec)
{
    std::transform(vec.begin(), vec.end(), vec.begin(), std::negate<double>());
    return std::move(vec);
}

VectorT plus(VectorT lhs, VectorT rhs)
{
    std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), std::plus<double>{});
    return std::move(lhs);
}

MatrixT plus(MatrixT lhs, MatrixT rhs)
{
    std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
        [](auto & l_line, auto & r_line) { return plus(std::move(l_line), std::move(r_line)); });
    return std::move(lhs);
}

VectorT minus(VectorT lhs, VectorT rhs)
{
    std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), std::minus<double>{});
    return std::move(lhs);
}

VectorT normalize(VectorT vec)
{
    return mul(std::move(vec), 1. / length(vec));
}


VectorT mul(VectorT vec, double scalar)
{
    std::transform(vec.begin(), vec.end(), vec.begin(), [scalar](double x_i) { return x_i * scalar; });
    return std::move(vec);
}

VectorT mul(const MatrixT & matrix, const VectorT & vec)
{
    VectorT res(vec.size());

    std::transform(matrix.begin(), matrix.end(), res.begin(), [&vec](const auto & line) {
        return scalar(line, vec);
    });
    return res;
}

MatrixT mul(const VectorT & lhs, const VectorT & rhs)
{
    MatrixT res(lhs.size());
    std::transform(lhs.begin(), lhs.end(), res.begin(), [&](double l_el) { return mul(rhs, l_el); });
    return res;
}

MatrixT mul(MatrixT matrix, double scalar)
{
    for (auto & line : matrix) {
        for (auto & el : line) {
            el *= scalar;
        }
    }
    return std::move(matrix);
}

MatrixT mul(MatrixT lhs, MatrixT rhs)
{
    VectorT line(lhs.size());
    rhs = util::trans(std::move(rhs));

    for (unsigned i = 0; i < lhs.size(); ++i) {
        for (unsigned j = 0; j < lhs.size(); ++j) {
            line[j] = scalar(lhs[i], rhs[j]);
        }
        lhs[i] = line;
    }
    return std::move(lhs);
}


MatrixT trans(MatrixT mtx)
{
    for (unsigned i = 0; i < mtx.size(); ++i) {
        for (unsigned j = 0; j < i; ++j) {
            std::swap(mtx[i][j], mtx[j][i]);
        }
    }
    return std::move(mtx);
}


double scalar(const VectorT & lhs, const VectorT & rhs)
{
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), 0.);
}

double length(const VectorT & vec)
{
    return scalar(vec, vec);
}

} // namespace util
