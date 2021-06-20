#include "sd_methods/Brent.h"

#include "util/ReplayData.h"
#include "util/VersionedData.h"

#include <cmath>

namespace min1d {
    /*
     * Helper functions for Brent's method.
     */
namespace {

bool all_different(double val1, double val2, double val3, double eps)
{
    return (std::abs(val1 - val2) > eps) &&
           (std::abs(val1 - val3) > eps) &&
           (std::abs(val2 - val3) > eps);
}

double count_parabole(std::pair<double, double> p1, std::pair<double, double> p2, std::pair<double, double> p3)
{
    double a0 = p1.second;
    double a1 = (p2.second - a0) / (p2.first - p1.first);
    double a2 = ((p3.second - a0) / (p3.first - p1.first) - a1) / (p3.first - p2.first);

    return (p1.first + p2.first - a1 / a2) / 2;
}

double sign(double val)
{
    return std::signbit(val) ? 1 : -1;
}

} // anonymous namespace

double Brent::find_min_impl() noexcept /*override*/
{
    const auto & fn = last_func();
    auto bnds = fn.bounds();
    const uint ITER_MAX = 100;

    /*
     * Start with finding three points x, w, v that will be used in constructing parabola and function values in this points.
     * On each iteration, if points x, w, v are good, try to construct quadratic polynomial as in parabole method.
     * Find minimum of the parabola - point u. If it lies within current segment and it is not very far from current function's minimum x, accept it.
     * Otherwise, use golden ratio method to find value of point u.
     * Choose next segment according to x, u and function values in this point.
     * Break, when the exit condition is met (it guarantess required precision) or when the iteration threshold is met.
     */
    double x, w, v;
    double f_x, f_w, f_v;
    x = w = v = bnds.from + TAU * bnds.length();
    f_x = f_w = f_v = fn(x);

    double step, prev_step;
    step = prev_step = bnds.length();
    for (uint iter = 0; iter < ITER_MAX; iter++) {
        double prev_prev_step = prev_step;
        prev_step = step;
        double to_leave = m_eps * std::abs(x) + m_eps / 10;
        if (std::abs(x - bnds.middle()) + bnds.length() / 2 - 2 * to_leave <= m_eps) {
            break;
        }

        double u;
        bool is_accepted = false;
        if (all_different(x, w, v, m_eps) && all_different(f_x, f_w, f_v, m_eps)) {
            /*
             * Points x, w, v and function values are good. Try to find parabola's minimum. 
             */
            u = count_parabole({x, f_x}, {w, f_w}, {v, f_v});
            if (bnds.from + m_eps <= u && u <= bnds.to - m_eps && std::abs(u - x) < prev_prev_step / 2) {
                /*
                 * Accept point got from parabolic interpolation. 
                 * If u is too close to segment's ends, place it between x and segment's middle.
                 */
                is_accepted = true;
                if (u - bnds.from < 2 * to_leave || bnds.to - u < 2 * to_leave) {
                    u = x - sign(x - bnds.middle()) * to_leave;
                }
            }
        }

        if (!is_accepted) {
            /*
             * Points x, w, v were not good or point got from parabolic interpolation was not accepted.
             * Use golden ratio method.
             */
            if (x < bnds.middle()) {
                u = x + TAU * (bnds.to - x);
                step = bnds.to - x;
            } else {
                u = x - TAU * (x - bnds.from);
                step = x - bnds.from;
            }
        }
        step = std::abs(u - x);
        double f_u = fn(u);
        if (f_u <= f_x) {
            if (u >= x) {
                /*
                 * Minimum of the function is to the right from x.
                 * Proceed with [x, b].
                 */
                bnds.from = x;
            } else { // u < x
                /*
                 * Minimum of the function is to the left from x.
                 * Proceed with [a, x].
                 */
                bnds.to = x;
            }

            /*
             * Save points according to rules.
             */
            v = w;
            w = x;
            x = u;
            f_v = f_w;
            f_w = f_x;
            f_x = f_u;
        } else { // f_u > f_x
            if (u >= x) {
                /*
                 * Minimum of the function is to the left from u.
                 * Proceed with [a, u].
                 */
                bnds.to = u;
            } else { // u < x
                /*
                 * Minimum of the function is to the right from u.
                 * Proceed with [u, b].
                 */
                bnds.from = u;
            }

            /*
             * Save points according to rules.
             */
            if (f_u <= f_w || w == x) {
                v = w;
                w = u;
                f_v = f_w;
                f_w = f_u;
            } else if (f_u <= f_v || v == x || v == w) {
                v = u;
                f_v = f_u;
            }
        }
    }

    return x;
}

} // namespace min1d
