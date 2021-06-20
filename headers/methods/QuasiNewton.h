#pragma once

#include "methods/Searcher.h"
#include "util/Function.h"

#include <vector>

struct QuasiNewton : Searcher
{
    using Searcher::Searcher;

    PointT bfs(const Function & func, PointT init = {});

private:
    double m_eps;
};
