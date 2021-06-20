#include "util/Function.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>

/*static*/ std::vector<int> Term::create_pows(unsigned dims, std::vector<std::pair<unsigned, int>> idx_pows)
{
    std::vector<int> pows(dims);
    for (auto [idx, pow] : idx_pows) {
        assert(idx < dims && "Variable indicies must be less than dims count");
        pows[idx] = pow;
    }
    return pows;
}

std::size_t Term::count_hash() const noexcept
{
    std::size_t res = 0;
    for (unsigned i = 0; i < m_pows.size(); ++i) {
        if (m_pows[i] != 0.) {
            res += i * m_pows[i];
        }
    }
    return res;
}

double Term::operator()(const std::vector<double> x) const noexcept
{
    assert(x.size() == dims() && "x dims must be equal to term's dims");

    return std::inner_product(x.begin(), x.end(), m_pows.begin(),
        m_coeff,
        std::multiplies<double>{},
        [](double x_i, int pow) { return pow == 0 ? 1. : std::pow(x_i, pow); });
}

double Function::operator()(const std::vector<double> &x) const noexcept
{
    assert(x.size() == dims() && "x dims must be equal to func's dims");

    return std::accumulate(m_terms.begin(), m_terms.end(),
        0.,
        [&x](double res, const Term & term) { return res + term(x); });
}
