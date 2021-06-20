#include "sole-solver/ProfileMatrix.h"
#include "sole-solver/Matrix.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>

ProfileMatrix::ProfileMatrix(const Matrix & matrix)
{
    assert(matrix.row_cnt() * matrix.col_cnt() > 0 && "Matrix must contain at least single row");
    assert(matrix.row_cnt() == matrix.col_cnt()
        && "Matrix must be square");

    id_t dim = matrix.row_cnt();
    diag.resize(dim);
    for (id_t i = 0; i < dim; i++)
    {
        diag[i] = matrix.get(i, i);
    }
    prof.resize(dim + 1);
    prof[0] = prof[1] = 1;
    for (id_t i = 1; i < dim; i++)
    {
        id_t pos = 0;
        for (; pos < i && matrix.get(i, pos) == 0; pos++);
        prof[i + 1] = prof[i] + i - pos;
        for (id_t j = pos; j < i; j++)
        {
            a_low.push_back(matrix.get(i, j));
            a_up.push_back(matrix.get(j, i));
        }
    }
}

ProfileMatrix::ProfileMatrix(
    value_vec diag,
    value_vec low,
    value_vec up,
    id_vec prof)
    : diag(std::move(diag))
    , a_low(std::move(low))
    , a_up(std::move(up))
    , prof(std::move(prof))
{}

/*static*/ ProfileMatrix ProfileMatrix::lu_decompose(Matrix && matrix)
{
    auto dims = matrix.row_cnt();

    assert(dims == matrix.col_cnt() && "Square matrix expected for LU-decomposition");

    auto subtract_from_cell = [&matrix](id_t row, id_t col, id_t count, value_t cell) {
        for (id_t k = 0; k < count; ++k) {
            cell -= matrix.get(row, k) * matrix.get(k, col);
        }
        return cell;
    };

    value_t new_val;
    for (id_t i = 0; i < dims; ++i) {
        for (id_t j = 0; j < dims; ++j) {
            new_val = matrix.get(i, j);
            new_val = subtract_from_cell(i, j, std::min(i, j), new_val);
            if (i > j) {
                new_val /= matrix.get(j, j);
            }
            matrix.set(i, j, new_val);
        }
    }

    return ProfileMatrix(matrix);
}

auto ProfileMatrix::get(id_t row, id_t col) const -> value_t /*override*/
{
    auto res = get_ptr(row, col);
    return res ? *res : 0;
}

void ProfileMatrix::set(id_t row, id_t col, value_t val) /*override*/
{
    auto changed = get_ptr(row, col);
    if (changed) {
        *changed = val;
    }
}

auto ProfileMatrix::get_ptr(id_t row, id_t col) const -> const value_t *
{
    assert(row >= 0 && col >= 0
        && row < diag.size() && col < diag.size()
        && "Row and column indexes must be in bounds to get");

    if (row == col)
    {
        return &diag[row];
    }

    bool is_lower = col < row;
    if (!is_lower)
    {
        std::swap(col, row);
    }

    id_t prof_val = prof[row + 1] - prof[row];
    if (prof_val < row - col)
    {
        return nullptr;
    }
    id_t pos = (prof[row] - 1) + prof_val - (row - col);
    return is_lower ? &a_low[pos] : &a_up[pos];
}

auto ProfileMatrix::get_ptr(id_t row, id_t col) -> value_t*
{
    return const_cast<value_t*>(const_cast<const ProfileMatrix*>(this)->get_ptr(row, col));
}

template<class T>
void print_vector(std::ostream& os, const std::vector<T> & vec)
{
    for (std::size_t i = 0; i < vec.size(); i++)
    {
        os << vec[i] << '\t';
    }
}
std::ostream& operator<<(std::ostream& os, const ProfileMatrix& pm)
{
    os << std::setprecision(20);
    print_vector(os, pm.diag);
    os << '\n';
    print_vector(os, pm.a_low);
    os << '\n';
    print_vector(os, pm.a_up);
    os << '\n';
    print_vector(os, pm.prof);
    os << '\n';
    return os;
}
