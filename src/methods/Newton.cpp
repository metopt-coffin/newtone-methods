#include "methods/Newton.h"

#include "sd_methods/Brent.h"
#include "sole-solver/QuadMatrix.h"
#include "sole-solver/Solver.h"
#include "util/VectorOps.h"
#include "util/VersionedData.h"

#include <functional>
#include <iostream>
#include <type_traits>

std::vector<double> NewtonMethods::classic(const Function & func, std::vector<double> init)
{
    // Init
    auto curr = init_method(func, std::move(init));
    auto eps_2 = m_eps * m_eps;

    // Count initial hessian and gradient
    auto grad = func.grad();
    auto hessian = grad.grad();

    // End if max iterations treshold is reached
    for (unsigned iter_num = 0; iter_num < MaxIter; ++iter_num) {
        log_x(iter_num, curr);

        // Solve sole to find p_k
        auto shift = Solver::solve_lu(QuadMatrix(hessian(curr)), util::neg(grad(curr)));

        auto shift_len = util::length(shift.answer);
        if (shift_len < eps_2) {
            // end iterating if we got enough precision
            break;
        } else {
            // count next point
            curr = util::add(std::move(curr), std::move(shift.answer));
        }
    }
    return curr;
}

std::vector<double> NewtonMethods::with_sd_search(const Function & func, std::vector<double> init)
{
    // Init
    auto curr = init_method(func, std::move(init));
    auto eps_2 = m_eps * m_eps;

    // Count initial hessian and gradient
    auto grad = func.grad();
    auto hessian = grad.grad();

    // End if max iterations treshold is reached
    for (unsigned iter_num = 0; iter_num < MaxIter; ++iter_num) {
        // Solve sole to find p_k
        auto shift = Solver::solve_lu(QuadMatrix(hessian(curr)), util::neg(grad(curr)));

        // Find coefficient alpha by solving one-dimensional minimization problem
        auto alpha = find_alpha(curr, shift.answer);

        log_x(iter_num, curr);
        log_alpha(iter_num, alpha);

        // Count the next step
        shift.answer = util::mul(std::move(shift.answer), alpha);
        if (util::length(shift.answer) < eps_2) {
            // end iterating if we got enough precision
            break;
        } else {
            // count the next point
            curr = util::add(std::move(curr), std::move(shift.answer));
        }
    }

    return curr;
}

std::vector<double> NewtonMethods::with_desc_dir(const Function & func, std::vector<double> init)
{
    // Init
    auto curr = init_method(func, std::move(init));
    auto eps_2 = m_eps * m_eps;

    // Count initial hessian and gradient
    auto grad = func.grad();
    auto hessian = grad.grad();

    // End if max iterations treshold is reached
    for (unsigned iter_num = 0; iter_num < MaxIter; ++iter_num) {
        // Count current gradient and antigradient
        auto curr_grad = grad(curr);
        auto curr_grad_neg = util::neg(curr_grad);

        // Solve sole to find p_k
        auto shift = Solver::solve_lu(QuadMatrix(hessian(curr)), std::vector(curr_grad_neg));

        // if current direction is not the descent direction, use antigradient instead
        if (util::scalar(shift.answer, curr_grad) > 0) {
            shift.answer = std::move(curr_grad_neg);
        }

        // Find coefficient alpha by solving one-dimensional minimization problem
        auto alpha = find_alpha(curr, shift.answer);

        log_x(iter_num, curr);
        log_alpha(iter_num, alpha);

        // Count the next step
        shift.answer = util::mul(std::move(shift.answer), alpha);
        if (util::length(shift.answer) < eps_2) {
            // end iterating if we got enough precision
            break;
        } else {
            // count the next point
            curr = util::add(std::move(curr), std::move(shift.answer));
        }
    }

    return curr;
}
