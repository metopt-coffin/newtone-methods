#pragma once

#include "MinSearcher.h"
#include "sd_methods/Brent.h"
#include "sd_methods/Function.h"
#include "sd_methods/MinSearcher.h"

#include "util/Function.h"
#include "util/Misc.h"
#include "util/VectorOps.h"

namespace min_nd {

struct FastestDescent : MinSearcher
{
    FastestDescent(double eps, double max_step = 10.)
        : m_sd_searcher(eps)
        , m_eps(eps)
        , m_alpha(max_step)
    {}

protected:
    /*
     * Find n-dimensional function's minimum 
     * using fastest descent method.
     */
    util::VectorT find_min_impl() override;
    /*
     * Find n-dimensional function's minimum
     * using fastest descent method.
     * Outputs tracing information.
     */
    util::VectorT find_min_traced_impl(util::VectorT init) override;

protected:
    /*
     * Solve the one dimensional minimization problem.
     */
    double find_sd_min(min1d::Function && func) { return m_sd_searcher.find_min(std::move(func)); }
    /*
     * Returns function for which one dimensional minimization problem was solved.
     */
    const min1d::Function & last_sd_func() const noexcept { return m_sd_searcher.last_func(); }

private:
    min1d::Brent m_sd_searcher;     // Current one dimensional minimization method
    double m_eps;
    double m_alpha;
};

} // namespace min_nd
