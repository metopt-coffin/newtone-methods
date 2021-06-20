#pragma once

#include "sole-solver/Matrix.h"

#include <iostream>
#include <vector>

class SparseMatrix : public Matrix
{
    using value_vec = std::vector<value_t>;
    using id_vec = std::vector<id_t>;
public:
    /*
     * Constructs a SparseMatrix from a standard two-dimensional square matrix.
     */
    explicit SparseMatrix(const Matrix& matrix);

    /*
     * Constructs a SparseMatrix directly from respective vectors
     * Note, that vectors are not checked to be valid vectors of profile matrix.
     */
    SparseMatrix(
        value_vec diag,
        value_vec low,
        value_vec up,
        id_vec j_prof,
        id_vec i_prof);

    /*
     * Copy and move constructors and assign operators.
     */
    SparseMatrix(const SparseMatrix& other) = default;
    SparseMatrix(SparseMatrix&& other) = default;
    SparseMatrix& operator=(const SparseMatrix& other) = default;
    SparseMatrix& operator=(SparseMatrix&& other) = default;

    /*
     * Returns value from "row" row and "col" column.
     * If the required cell is contained, than the corresponding value is returned.
     * Otherwise, returns 0.
     */
    value_t get(id_t row, id_t col) const override;

    /*
     * Sets value "val" to "row" row and "col" column.
     * If the required cell is contained, than the corresponding value is set.
     * Otherwise, does nothing.
     */
    void set(id_t row, id_t col, value_t val) override;

    id_t row_cnt() const override { return diag.size(); }
    id_t col_cnt() const override { return diag.size(); }

    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix& pm);
private:
    const value_t* get_ptr(id_t row, id_t col) const;
    value_t* get_ptr(id_t row, id_t col);

private:
    value_vec diag;                                                 // vector, containing diagonal elements
    value_vec a_low;                                                // vector, containing elements from the lower part of matrix
    value_vec a_up;                                                 // vector, containing elements from the upper part of matrix
    id_vec j_prof;                                                  // vector, containing number of columns (rows) of lower (upper) part of the matrix
    id_vec i_prof;                                                  // vector, containing first indices of non-zero elements in row
};