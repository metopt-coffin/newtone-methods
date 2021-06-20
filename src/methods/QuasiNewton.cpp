#include "methods/QuasiNewton.h"
#include "sd_methods/Brent.h"
#include "sole-solver/QuadMatrix.h"
#include "util/VectorOps.h"

namespace {

auto identity_matrix(unsigned dims)
{
    std::vector res(dims, std::vector<double>(dims));
    for (unsigned i = 0; i < dims; ++i) {
        res[i][i] = 1.;
    }
    return res;
}
}

auto QuasiNewton::bfs(const Function &func, PointT init) -> PointT
{
    using VectorT = std::vector<double>;
    using MatrixT = std::vector<VectorT>;

    double eps_2 = m_eps * m_eps;

    PointT curr = init_method(func, std::move(init));
    VectorT curr_diff(func.dims()), prev_w(func.dims());

    auto grad = func.grad();
    auto anti_hessian = identity_matrix(func.dims());

    auto next_anti_hessian = [](MatrixT ah, VectorT w_diff, const VectorT & curr_diff) {
        VectorT ah_wd = util::mul(ah, w_diff);
        double roe = util::scalar(ah_wd, w_diff);
        double wd_cd = util::scalar(w_diff, curr_diff);

        VectorT r = util::minus(
            util::mul(ah_wd, 1. / roe),
            util::mul(curr_diff, 1. / wd_cd));

        MatrixT fst = util::mul(util::mul(curr_diff, curr_diff), -1. / wd_cd);
        MatrixT sec = util::mul(util::mul(util::mul(ah_wd, w_diff), util::trans(ah)), -1. / roe);
        MatrixT thd = util::mul(util::mul(r, roe), std::move(r));

        return util::plus(
            util::plus(std::move(ah), std::move(fst)),
            util::plus(std::move(sec), std::move(thd)));
    };

    do {
        VectorT w = util::negate(grad(curr));
        VectorT p = util::mul(anti_hessian, w);
        double alpha = find_alpha(curr, p);

        curr_diff = util::mul(std::move(p), alpha);
        VectorT w_diff = util::minus(w, std::move(prev_w));
        prev_w = std::move(w);

        anti_hessian = next_anti_hessian(std::move(anti_hessian), std::move(w_diff), curr_diff);
        curr = util::plus(std::move(curr), curr_diff);
    } while(util::length(curr_diff) > eps_2);

    return curr;
}
