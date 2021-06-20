#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

/*
 * Abstract class representing two-dimensional matrix
 */
class Matrix
{
protected:
    using value_t = double;
    using id_t = std::size_t;
public:
    virtual ~Matrix() = default;
    /*
     * Returns value from "row" row and "col" column.
     */
    virtual value_t get(id_t row, id_t col) const = 0;

    /* 
     * Sets value "val" to "row" row and "col" column.
     */
    virtual void set(id_t row, id_t col, value_t val) = 0;

    virtual id_t row_cnt() const = 0;
    virtual id_t col_cnt() const = 0;

    std::vector<value_t> operator*(const std::vector<value_t>& vec) const;
};