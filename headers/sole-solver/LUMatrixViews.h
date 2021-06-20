#pragma once

#include "sole-solver/Matrix.h"

#include <stdexcept>

namespace detail {

struct MatrixView : public Matrix
{
    MatrixView(const Matrix & from)
        : matrix(from)
    {}

    void set(id_t row, id_t col, value_t val) override {}

    id_t row_cnt() const override { return matrix.col_cnt(); }
    id_t col_cnt() const override { return matrix.row_cnt(); }

protected:
    const Matrix & matrix;
};

} // namespace detail

class LMatrixView : public detail::MatrixView
{
    using Super = detail::MatrixView;
    using Super::Super;

public:
    value_t get(id_t row, id_t col) const override
    {
        if (row > col) {
            return matrix.get(row, col);
        } else if (row == col) {
            return value_t{1.};
        } else {
            return value_t{};
        }
    }

    void set(id_t row, id_t col, value_t val) override { throw std::runtime_error{"LMatrixView should always be const"}; }
};

class UMatrixView : public detail::MatrixView
{
    using Super = detail::MatrixView;
    using Super::Super;

public:
    value_t get(id_t row, id_t col) const override
    {
        if (row <= col) {
            return matrix.get(row, col);
        } else {
            return value_t{};
        }
    }
};
