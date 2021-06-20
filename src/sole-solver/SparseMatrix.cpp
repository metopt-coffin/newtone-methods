#include "sole-solver/SparseMatrix.h"

#include <cassert>
#include <iomanip>

SparseMatrix::SparseMatrix(const Matrix & matrix)
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
    i_prof.resize(dim + 1);
    i_prof[0] = i_prof[1] = 1;
    for (id_t i = 1; i < dim; i++)
    {
        id_t cnt = 0;
        for (int pos = 0; pos < i; pos++)
        {
            if (matrix.get(i, pos) != 0)
            {
                cnt++;
                a_low.push_back(matrix.get(i, pos));
                a_up.push_back(matrix.get(pos, i));
                j_prof.push_back(pos);
            }
        }
        i_prof[i + 1] = i_prof[i] + cnt;
    }
}

SparseMatrix::SparseMatrix(
    value_vec diag,
    value_vec low,
    value_vec up,
    id_vec j_profile,
    id_vec i_profile)
    : diag(std::move(diag))
    , a_low(std::move(low))
    , a_up(std::move(up))
    , j_prof(std::move(j_profile))
    , i_prof(std::move(i_profile))
{}

auto SparseMatrix::get(id_t row, id_t col) const -> value_t /*override*/
{
    auto res = get_ptr(row, col);
    return res ? *res : 0;
}

void SparseMatrix::set(id_t row, id_t col, value_t val) /*override*/
{
    auto changed = get_ptr(row, col);
    if (changed)
    {
        *changed = val;
    }
}

auto SparseMatrix::get_ptr(id_t row, id_t col) const -> const value_t*
{
    assert(row >= 0 && col >= 0
           && row < diag.size() && col < diag.size()
           && "Row and column indexes must be in bounds to get");

    if (row == col)
    { return &diag[row]; }

    bool is_lower = col < row;
    if (!is_lower)
    {
        std::swap(col, row);
    }

    id_t prof_val = i_prof[row + 1] - i_prof[row];
    id_t start_pos = i_prof[row] - 1;
    for (id_t i = start_pos; i < start_pos + prof_val; i++)
    {
        if (j_prof[i] == col)
        {
            return is_lower ? &a_low[i] : &a_up[i];
        }
    }
    return nullptr;
}

auto SparseMatrix::get_ptr(id_t row, id_t col) -> value_t*
{
    return const_cast<value_t*>(const_cast<const SparseMatrix*>(this)->get_ptr(row, col));
}

template<class T>
void print_vector(std::ostream& os, const std::vector<T>& vec)
{
    for (std::size_t i = 0; i < vec.size(); i++)
    {
        os << vec[i] << '\t';
    }
}
std::ostream& operator<<(std::ostream& os, const SparseMatrix& sm)
{
    os << std::setprecision(20);
    print_vector(os, sm.diag);
    os << '\n';
    print_vector(os, sm.a_low);
    os << '\n';
    print_vector(os, sm.a_up);
    os << '\n';
    print_vector(os, sm.j_prof);
    os << '\n';
    print_vector(os, sm.i_prof);
    os << '\n';
    return os;
}
