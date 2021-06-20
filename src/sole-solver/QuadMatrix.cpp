#include "sole-solver/QuadMatrix.h"

#include <cassert>
#include <iomanip>

QuadMatrix::QuadMatrix(std::vector<value_vec> other_matrix)
    : matrix(std::move(other_matrix))
{}

QuadMatrix::QuadMatrix(id_t dim1, id_t dim2)
    : matrix(dim1, value_vec(dim2))
{}

auto QuadMatrix::get_matrix() const -> const std::vector<value_vec> &
{
    return matrix;
}

auto QuadMatrix::get(id_t row, id_t col) const -> value_t /*override*/
{
    assert(row >= 0 && col >= 0 && row < matrix.size() && (!matrix.size() || col < matrix.front().size())
        && "Indexes out of bounds in QuadMatrix get");
    return matrix[row][col];
}

void QuadMatrix::set(id_t row, id_t col, value_t val) /*override*/
{
    assert(row >= 0 && col >= 0 && row < matrix.size() && (!matrix.size() || col < matrix.front().size())
        && "Indexes out of bounds in QuadMatrix set");
    matrix[row][col] = val;
}

std::ostream& operator<<(std::ostream& os, const QuadMatrix& qm)
{
    os << std::setprecision(20);
    int size = qm.matrix.size();
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            os << qm.get(i, j) << '\t';
        }
        os << '\n';
    }
    return os;
}