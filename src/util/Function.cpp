#include "util/Function.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>


std::unique_ptr<Fn> cns(double value) { return std::make_unique<Const>(value); }
std::unique_ptr<Fn> var(unsigned idx) { return std::make_unique<Variable>(idx); }

std::unique_ptr<Fn> add(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r)
{
    std::optional<double> l_as_const;
    if (l->dims() == 0) {
        l_as_const = l->as_const();

        if (*l_as_const == 0.) {
            return std::move(r);
        }
    }
    if (r->dims() == 0) {
        double cnst = r->as_const();

        if (cnst == 0.) {
            return std::move(l);
        } else if (l_as_const) {
            return cns(*l_as_const + cnst);
        }
    }
    return std::make_unique<Add>(std::move(l), std::move(r));
}

std::unique_ptr<Fn> sub(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r)
{
    std::optional<double> r_as_const;
    if (r->dims() == 0) {
        r_as_const = r->as_const();
        if (*r_as_const == 0.) {
            return std::move(l);
        }
    }

    if (r_as_const && l->dims() == 0) {
        return cns(l->as_const() - *r_as_const);
    }
    return std::make_unique<Sub>(std::move(l), std::move(r));
}

std::unique_ptr<Fn> mul(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r)
{
    std::optional<double> l_as_const;
    if (l->dims() == 0) {
        double cnst = l->as_const();
        if (cnst == 0.) {
            return cns(0.);
        } else if (cnst == 1.) {
            return std::move(r);
        }
        l_as_const = cnst;
    }

    if (r->dims() == 0) {
        double cnst = r->as_const();
        if (cnst == 0.) {
            return cns(0.);
        } else if (cnst == 1.) {
            return std::move(l);
        } else if (l_as_const) {
            return cns(*l_as_const * cnst);
        }
    }

    return std::make_unique<Mul>(std::move(l), std::move(r));
}

std::unique_ptr<Fn> pow(std::unique_ptr<Fn> base, int p)
{
    if (p == 0) {
        return cns(1.);
    } else if (p == 1) {
        return std::move(base);
    }

    if (base->dims() == 0) {
        return cns(std::pow(base->as_const(), p));
    }

    return std::make_unique<Pow>(std::move(base), p);
}


std::unique_ptr<Fn> Const::part_der(unsigned idx) const noexcept
{
    return cns(0.);
}

std::unique_ptr<Fn> Variable::part_der(unsigned idx) const noexcept
{
    if (index == idx) {
        return cns(1.);
    } else {
        return cns(0.);
    }
}

std::unique_ptr<Fn> Add::part_der(unsigned idx) const noexcept
{
    auto l_gr = m_l->part_der(idx);
    auto r_gr = m_r->part_der(idx);

    return add(std::move(l_gr), std::move(r_gr));
}

std::unique_ptr<Fn> Sub::part_der(unsigned idx) const noexcept
{
    auto l_gr = m_l->part_der(idx);
    auto r_gr = m_r->part_der(idx);

    return sub(std::move(l_gr), std::move(r_gr));
}

std::unique_ptr<Fn> Mul::part_der(unsigned idx) const noexcept
{
    auto l_gr = m_l->part_der(idx);
    auto r_gr = m_r->part_der(idx);

    return (m_l->part_der(idx) * m_r->clone()) + (m_l->clone() * m_r->part_der(idx));
}

std::unique_ptr<Fn> Pow::part_der(unsigned idx) const noexcept
{
    return cns(m_pow) * m_base->part_der(idx) * (m_base->clone() ^ (m_pow - 1));
}
