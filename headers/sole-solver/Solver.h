#pragma once

#include "sole-solver/Matrix.h"

#include <vector>

class Solver
{
public:
    using value_t = double;
    using id_t = std::size_t;

    struct Result
    {
        /*
         * Special values for different error cases
         */
        constexpr static id_t FAILED = -1;
        constexpr static id_t UNTRUSTED = -2;

        std::vector<value_t> answer;
        id_t actions;
    };

    /*
     * Solve a system of linear equations using Gauss with LU-decomposition method
     */
    static Result solve_lu(Matrix && a, std::vector<value_t> && b);
    /*
     * Solve a system of linear equations using Gauss with pivot element choice
     */
    static Result solve_gauss(Matrix && a, std::vector<value_t> && b);

};
