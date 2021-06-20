#include "methods/Newton.h"

#include "sd_methods/Brent.h"
#include "sole-solver/QuadMatrix.h"
#include "sole-solver/Solver.h"
#include "util/VectorOps.h"
#include "util/VersionedData.h"

#include <functional>
#include <type_traits>

std::vector<double> NewtonMethods::classic(const Function & func, std::vector<double> init)
{
    auto curr = init_method(func, std::move(init));
    auto eps_2 = m_eps * m_eps;

    auto grad = func.grad();
    auto hessian = grad.grad();

    for (unsigned iter_num = 0; iter_num < MaxIter; ++iter_num) {
        log_x(iter_num, curr);

        auto shift = Solver::solve_lu(QuadMatrix(hessian(curr)), util::negate(grad(curr)));

        auto shift_len = util::length(shift.answer);
        if (shift_len < eps_2) {
            break;
        } else {
            curr = util::plus(std::move(curr), std::move(shift.answer));
        }
    }
    return curr;
}

std::vector<double> NewtonMethods::with_sd_search(const Function & func, std::vector<double> init)
{
    auto curr = init_method(func, std::move(init));
    auto eps_2 = m_eps * m_eps;

    auto grad = func.grad();
    auto hessian = grad.grad();

    for (unsigned iter_num = 0; iter_num < MaxIter; ++iter_num) {
        auto shift = Solver::solve_lu(QuadMatrix(hessian(curr)), util::negate(grad(curr)));

        auto alpha = find_alpha(curr, shift.answer);

        log_x(iter_num, curr);
        log_alpha(iter_num, alpha);

        shift.answer = util::mul(std::move(shift.answer), alpha);
        if (util::length(shift.answer) < eps_2) {
            break;
        } else {
            curr = util::plus(std::move(curr), std::move(shift.answer));
        }
    }

    return curr;
}

std::vector<double> NewtonMethods::with_desc_dir(const Function & func, std::vector<double> init)
{
    auto curr = init_method(func, std::move(init));
    auto eps_2 = m_eps * m_eps;

    auto grad = func.grad();
    auto hessian = grad.grad();

    for (unsigned iter_num = 0; iter_num < MaxIter; ++iter_num) {
        auto curr_grad = grad(curr);
        auto curr_grad_neg = util::negate(curr_grad);

        auto shift = Solver::solve_lu(QuadMatrix(hessian(curr)), std::vector(curr_grad_neg));

        if (util::scalar(shift.answer, curr_grad) > 0) {
            shift.answer = std::move(curr_grad_neg);
        }

        auto alpha = find_alpha(curr, shift.answer);

        log_x(iter_num, curr);
        log_alpha(iter_num, alpha);

        shift.answer = util::mul(std::move(shift.answer), alpha);
        if (util::length(shift.answer) < eps_2) {
            break;
        } else {
            curr = util::plus(std::move(curr), std::move(shift.answer));
        }
    }

    return curr;
}
