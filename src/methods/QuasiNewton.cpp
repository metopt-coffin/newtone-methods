#include "methods/QuasiNewton.h"
#include "sd_methods/Brent.h"
#include "sole-solver/QuadMatrix.h"
#include "util/VectorOps.h"
#include <functional>


using VectorT = util::VectorT;
using MatrixT = util::MatrixT;

namespace {

MatrixT identity_matrix(unsigned dims)
{
    std::vector res(dims, std::vector(dims, 0.));
    for (unsigned i = 0; i < dims; ++i) {
        res[i][i] = 1.;
    }
    return res;
}
}

MatrixT QuasiNewton::bfs(MatrixT ah, VectorT w_diff, const VectorT & curr_diff)
{
    VectorT ah_wd = util::mul(ah, w_diff);
    double roe = util::scalar(ah_wd, w_diff);
    double wd_cd = util::scalar(w_diff, curr_diff);

    VectorT r = util::sub(
        util::mul(ah_wd, 1. / roe),
        util::mul(curr_diff, 1. / wd_cd));

    MatrixT fst = util::mul(util::mul(curr_diff, curr_diff), 1. / wd_cd);
    MatrixT sec = util::mul(util::mul(ah_wd, ah_wd), 1. / roe);
    MatrixT thd = util::mul(util::mul(r, r), roe);

    ah = util::sub(std::move(ah), std::move(fst));
    ah = util::sub(std::move(ah), std::move(sec));
    return util::add(std::move(ah), std::move(thd));
}

util::MatrixT QuasiNewton::powell(util::MatrixT ah, util::VectorT w_diff, const util::VectorT & x_diff)
{
    VectorT x_wave = util::add(util::mul(ah, w_diff), x_diff);

    return util::sub(
        std::move(ah),
        util::mul(util::mul(x_wave, x_wave), 1. / util::scalar(w_diff, x_wave)));
}

auto QuasiNewton::search_common(UpdateRule rule, const Function &func, PointT init) -> PointT
{
    decltype(&QuasiNewton::bfs) next_anti_hessian_method;
    switch (rule) {
        case UpdateRule::BFSh:
            next_anti_hessian_method = &QuasiNewton::bfs;
            break;
        case UpdateRule::Powell:
            next_anti_hessian_method = &QuasiNewton::powell;
            break;
    }

    using namespace std::placeholders;
    auto next_anti_hessian = std::bind(next_anti_hessian_method, this, _1, _2, _3);

    double eps_2 = m_eps * m_eps;

    PointT curr = init_method(func, std::move(init));
    VectorT curr_diff, w;

    unsigned iter_num = 0;
    log_x(iter_num++, curr);

    auto grad = func.grad();
    MatrixT anti_hessian = identity_matrix(func.dims());

    {
        w = util::neg(grad(curr));
        VectorT p = w;
        double alpha = find_alpha(curr, p);
        curr_diff = util::mul(std::move(p), alpha);
        curr = util::add(std::move(curr), curr_diff);
        log_x(iter_num++, curr);
    }

    do {
        VectorT next_w = util::neg(grad(curr));
        VectorT w_diff = util::sub(next_w, std::move(w));
        w = std::move(next_w);

        anti_hessian = next_anti_hessian(std::move(anti_hessian), std::move(w_diff), curr_diff);

        VectorT p = util::mul(anti_hessian, w);
        double alpha = find_alpha(curr, p);

        curr_diff = util::mul(std::move(p), alpha);
        curr = util::add(std::move(curr), curr_diff);
        log_x(iter_num++, curr);
    } while(util::length(curr_diff) > eps_2);

    return curr;
}
