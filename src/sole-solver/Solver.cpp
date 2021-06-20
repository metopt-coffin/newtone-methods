#include "sole-solver/Solver.h"
#include "sole-solver/LUMatrixViews.h"
#include "sole-solver/ProfileMatrix.h"

#include <cmath>
#include <numeric>

/*static*/ auto Solver::solve_gauss(Matrix&& a, std::vector<value_t>&& b) -> Result
{
    id_t actions_cnt = 0;                                           // counter for mult and div operations
    /*
     * We don't know what kind of matrix we received.
     * Therefore, we can't swap rows in-place and have to remember all permutations we made.
     */
    std::vector<id_t> permutations(b.size());
    std::iota(permutations.begin(), permutations.end(), 0);         // permutations[i] = i;

    // straight
    for (int i = 0; i < b.size(); i++)
    {
        /*
         * Starting from here, 
         *  "virtual" - indexes in matrix with swapped rows
         *  "physical" - indexes in real matrix.
         * Usually, we don't have to map column indexes as they remain the same. 
         */
        id_t i_row = permutations[i];                               // get physical index of current row
        value_t pivot = a.get(i_row, i);                            // and current max pivot element
        id_t pivot_row = i;                                         // remember the virtual index of max pivot row

        /*
         * Find the max pivot element.
         * Look up in the next rows, because previous rows are already processed
         */
        for (int j = i + 1; j < b.size(); j++)
        {
            id_t new_row = permutations[j];                         // get physical index of row we look at now
            value_t new_val = a.get(new_row, i);                    // get next value from physical matrix
            if (std::abs(new_val) > std::abs(pivot))                // relaxate maximum
            {
                pivot = new_val;
                pivot_row = j;
            }
        }
        if (std::abs(pivot) < 1e-20)                                // if pivot element is zero, exit. Hope that 10^-20 is enough for epsilon.
        { return { b, Result::FAILED }; }                           //    (just a thought: probably, should add epsilon to the function's signature?)

        std::swap(permutations[i], permutations[pivot_row]);        // "swap" rows in virtual matrix
        i_row = permutations[i];                                    // get physical index of current row (after swap)
        /*
         * Perform row operations and process the rest of the matrix
         */
        for (int j = i + 1; j < b.size(); j++)
        {
            id_t j_row = permutations[j];                           // get physical index of row we change now
            value_t factor = a.get(j_row, i) / a.get(i_row, i);     // get factor that will be used to change values in this row
            b[j_row] -= factor * b[i_row];                          // change the vector b
            actions_cnt += 2;
            /*
             * Process the row using formula: 
             *  a[i][j] = a[i][j] - factor * a[k][j], where k - index of row we subtract
             */
            for (int k = i + 1; k < b.size(); k++)
            {
                value_t val = a.get(j_row, k);                      // get previous value
                a.set(j_row, k, val - factor * a.get(i_row, k));    // set new value to the same place
                actions_cnt++;
            }
        }
    }

    
    // reverse
    /*
     * Now we have an upper triangular matrix. 
     * To get the solution we have to transfrom it to diagonal matrix, 
     *  changing the vector b accordingly
     */
    for (int i = b.size() - 1; i >= 0; i--)
    {
        id_t row = permutations[i];                                 // get physical index of current row
        for (int j = b.size() - 1; j > i; j--)
        {
            b[row] -= a.get(row, j) * b[permutations[j]];           // subtract all elements a[i][j] * x[j], where j > i;
            actions_cnt++;
        }                                                           //  we know answer x[j] as it is the already counted value in the vector b
        b[row] /= a.get(row, i);                                    // divide by diagonal element a[i][i];
        actions_cnt++;
    }
    
    /*
     * Write down the result vector in correct order
     */
    std::vector<value_t> answer(b.size());
    for (int i = 0; i < b.size(); i++)
    {
        answer[i] = b[permutations[i]];
    }
    return { answer, actions_cnt };
}

/*static*/ auto Solver::solve_lu(Matrix && a, std::vector<value_t> && b) -> Result
{
    bool untrusted_flag = false;                                    // flag for division by near-epsilon number
    id_t actions_cnt = 0;                                           // counter for mult and div operations
    /*
     * As we have LU-decomposition of matrix solving a system of linear equations is easy.
     * First, we have to solve equation L * y = b, where L - lower triangular, b - given vector of values,
     *      and, thus, find vector y.
     * Second, we solve equation U * x = y, where U - upper triangular. Thus, we find answer - vector x.
     */
    ProfileMatrix pm = ProfileMatrix::lu_decompose(std::move(a));
    LMatrixView l(pm);
    UMatrixView u(pm);

    /*
     * L is a lower triangular matrix.
     * Therefore, all we have to do is to transform it to diagonal form.
     * We can do it by using Gaussian elimination algorithm and only in a straight way.
     * Moreover, we don't have to change the matrix, we need to change only vector b.
     * Note, that we don't have to divide by these elements, as they are equal to one.
     */
    for (int i = 0; i < b.size(); i++)
    {
        for (int j = i + 1; j < b.size(); j++)
        {
            b[j] -= b[i] * l.get(j, i);
            actions_cnt++;
        }
    }

    /*
     * U is an upper triangular matrix.
     * We also have to transform it and we will use Gaussian elimination again,
     *      but in a reversed way now.
     */
    for (int i = b.size() - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < b.size(); j++)
        {
            b[i] -= u.get(i, j) * b[j];
            actions_cnt++;
        }
        b[i] /= u.get(i, i);
        actions_cnt++;
    }

    return { b, actions_cnt };
}
