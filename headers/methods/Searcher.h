#pragma once

#include "sd_methods/Brent.h"
#include "util/Function.h"
#include "util/ReplayData.h"

/*
 * Base class for Newton methods implementations
 */
struct Searcher
{
    using PointT = std::vector<double>;

    Searcher(double eps)
        : m_eps(eps)
        , m_sd_searcher(m_eps)
    {}

    const Function & last_func() const noexcept { return *m_last_func; }
    const util::ReplayData & replay_data() const noexcept { return m_replay_data; }

protected:
    // Initialize values before starting method
    PointT init_method(const Function & func, PointT init);

    // Find coefficient alpha by solving one-dimensional minimization problem
    double find_alpha(const PointT & curr, const std::vector<double> & shift);

    // Log current point
    void log_x(unsigned iter_num, const std::vector<double> & x);
    // Log current alpha coefficient
    void log_alpha(unsigned iter_num, double alpha);

protected:
    static constexpr unsigned MaxIter = 3000;       // Maximum iterations treshold

    double m_eps;                                   // Current precision

    const Function * m_last_func = nullptr;         // Minimum of this function is searched now
    util::ReplayData m_replay_data;                 // Object for recording tracing information
    min1d::Brent m_sd_searcher;                     // One-dimensional minimization problem solver
};
