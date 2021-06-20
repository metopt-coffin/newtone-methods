#include "util/VectorOps.h"

#include <algorithm>
#include <functional>
#include <numeric>

namespace util {

std::vector<double> negate(std::vector<double> vec)
{
    std::transform(vec.begin(), vec.end(), vec.begin(), std::negate<double>());
    return std::move(vec);
}

std::vector<double> plus(std::vector<double> lhs, std::vector<double> rhs)
{
    std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), std::plus<double>{});
    return std::move(lhs);
}

std::vector<double> mul(std::vector<double> vec, double scalar)
{
    std::transform(vec.begin(), vec.end(), vec.begin(), [scalar](double x_i) { return x_i * scalar; });
    return std::move(vec);
}

std::vector<double> normalize(std::vector<double> vec)
{
    return mul(std::move(vec), 1. / length(vec));
}

double scalar(const std::vector<double> & lhs, const std::vector<double> & rhs)
{
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), 0.);
}

double length(const std::vector<double> & vec)
{
    return scalar(vec, vec);
}

} // namespace util
