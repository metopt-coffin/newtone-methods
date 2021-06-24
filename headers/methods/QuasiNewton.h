#pragma once

#include "methods/Searcher.h"
#include "util/Function.h"
#include "util/VectorOps.h"

#include <vector>

/*
 * Quasinewton methods' implementation
 */
struct QuasiNewton : Searcher
{
    using Searcher::Searcher;

private:
    // enum to trace current method
    enum struct UpdateRule
    {
        BFSh,
        Powell,
    };
    // Broyden-Fletcher-Shanno algorithm
    util::MatrixT bfs(util::MatrixT anti_hessian, util::VectorT w_diff, const util::VectorT & x_diff);
    // Powell algorithm
    util::MatrixT powell(util::MatrixT anti_hessian, util::VectorT w_diff, const util::VectorT & x_diff);

    // Common code for all quasinewton methods
    PointT search_common(UpdateRule rule, const Function & func, PointT init);

public:
    // public wrappers for methods
    PointT search_bfs(const Function & func, PointT init = {}) { return search_common(UpdateRule::BFSh, func, std::move(init)); }
    PointT search_powell(const Function & func, PointT init = {}) { return search_common(UpdateRule::Powell, func, std::move(init)); }
};
