#include <cassert>

#include "sole-solver/Matrix.h"

auto Matrix::operator*(const std::vector<value_t>& vec) const -> std::vector<value_t>
{
    assert(col_cnt() == vec.size() && "Dimension mismatch in matrix * vector");
    std::vector<value_t> answer(row_cnt(), 0.);
    for (id_t row = 0; row < answer.size(); row++)
    {
        for (id_t column = 0; column < vec.size(); column++)
        {
            answer[row] += get(row, column) * vec[column];
        }
    }
    return answer;
}