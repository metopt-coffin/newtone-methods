#pragma once

#include "methods/Searcher.h"
#include "util/Function.h"
#include "util/VectorOps.h"

#include <vector>

struct QuasiNewton : Searcher
{
    using Searcher::Searcher;

private:
    enum struct UpdateRule
    {
        BFSh,
        Powell,
    };
    util::MatrixT bfs(util::MatrixT anti_hessian, util::VectorT w_diff, const util::VectorT & x_diff);
    util::MatrixT powell(util::MatrixT anti_hessian, util::VectorT w_diff, const util::VectorT & x_diff);

    PointT search_common(UpdateRule rule, const Function & func, PointT init);

public:
    PointT search_bfs(const Function & func, PointT init = {}) { return search_common(UpdateRule::BFSh, func, std::move(init)); }
    PointT search_powell(const Function & func, PointT init = {}) { return search_common(UpdateRule::Powell, func, std::move(init)); }
};
