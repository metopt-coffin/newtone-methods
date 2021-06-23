#include "nd_methods/FastestDescent.h"

#include "nd_methods/MinSearcher.h"
#include "sd_methods/MinSearcher.h"

#include "util/Function.h"
#include "util/VectorOps.h"
#include "util/VersionedData.h"

namespace min_nd {
/*
 * Idea: after finding gradient of the function do not make a small step in the direction of the antigradient.
 * Instead, move, until the function decreases. 
 * After finding minimum on the chosen direction, find function's gradient again. 
 * Repeat the algorithm.
 */
util::VectorT FastestDescent::find_min_impl()
{
    const double eps_pow2 = m_eps * m_eps;
    auto & func = last_func();

    util::VectorT curr(func.dims()); // Vector of current coordinates
    double f_curr = func(curr);

    auto grad = func.grad();
    util::VectorT shift = grad(curr);
    double sd_min;    // Minimum found on the chosen direction
    uint iter_num = 0;      // To prevent infinite or very long cycles
    while (util::length(shift) >= eps_pow2 && iter_num < MAX_ITER) {
        sd_min = find_sd_min({[&](double x) { return func(util::sub(curr, util::mul(shift, x))); }, {0., m_alpha}});
        curr = util::sub(std::move(curr), util::mul(std::move(shift), sd_min));
        shift = grad(curr);
        iter_num++;
    }

    return curr;
}

/*
 * Version with tracing output.
 */
util::VectorT FastestDescent::find_min_traced_impl(util::VectorT init)
{
    const double eps_pow2 = m_eps * m_eps;
    auto & func = last_func();

    util::VectorT curr(std::move(init));
    double f_curr = func(curr);

    auto grad = func.grad();
    util::VectorT shift = grad(curr);
    double sd_min;

    uint iter_num = 0;
    while (util::length(shift) >= eps_pow2 && iter_num < MAX_ITER) {;
        m_replay_data.emplace_back<util::VdPoint>(iter_num, curr);

        sd_min = find_sd_min({[&](double x) { return func(util::sub(curr, util::mul(shift, x))); }, {0., m_alpha}});

        curr = util::sub(std::move(curr), util::mul(std::move(shift), sd_min));
        shift = grad(curr);

        iter_num++;
    }

    return curr;
}

} // namespace min_nd
