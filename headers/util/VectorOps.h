/*
 *  Functions for matrix and vector operations
 */
#pragma once

#include <vector>

namespace util {

using VectorT = std::vector<double>;
using MatrixT = std::vector<VectorT>;

VectorT neg(VectorT vec);
VectorT add(VectorT lhs, VectorT rhs);
MatrixT add(MatrixT lhs, MatrixT rhs);
VectorT sub(VectorT lhs, VectorT rhs);
MatrixT sub(MatrixT lhs, MatrixT rhs);
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
