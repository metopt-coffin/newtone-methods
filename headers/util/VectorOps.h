#pragma once

#include <vector>

namespace util {

using VectorT = std::vector<double>;
using MatrixT = std::vector<VectorT>;

VectorT negate(VectorT vec);
VectorT plus(VectorT lhs, VectorT rhs);
MatrixT plus(MatrixT lhs, MatrixT rhs);
VectorT minus(VectorT lhs, VectorT rhs);
VectorT normalize(VectorT vec);

VectorT mul(VectorT vec, double scalar);
VectorT mul(const MatrixT & matrix, const VectorT & vec);
MatrixT mul(const VectorT & lhs, const VectorT & rhs);
MatrixT mul(MatrixT matrix, double scalar);
MatrixT mul(MatrixT lhs, MatrixT rhs);

MatrixT trans(MatrixT mtx);

double scalar(const VectorT & lhs, const VectorT & rhs);
double length(const VectorT & vec);

} // namespace util
