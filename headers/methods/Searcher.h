#pragma once

#include "util/Function.h"
#include "util/ReplayData.h"

struct Searcher
{
    using PointT = std::vector<double>;

    Searcher(double eps)
        : m_eps(eps)
    {}


    const Function & last_func() const noexcept { return *m_last_func; }
    const util::ReplayData & replay_data() const noexcept { return m_replay_data; }

protected:
    PointT init_method(const Function & func, PointT init);

    void log_x(unsigned iter_num, const std::vector<double> & x);
    void log_alpha(unsigned iter_num, double alpha);

protected:
    static constexpr unsigned MaxIter = 3000;

    double m_eps;

    const Function * m_last_func = nullptr;
    util::ReplayData m_replay_data;
};
