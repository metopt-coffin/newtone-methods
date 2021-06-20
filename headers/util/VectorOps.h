#pragma once

#include <vector>

namespace util {

std::vector<double> negate(std::vector<double> vec);
std::vector<double> plus(std::vector<double> lhs, std::vector<double> rhs);
std::vector<double> mul(std::vector<double> vec, double scalar);
std::vector<double> normalize(std::vector<double> vec);

double scalar(const std::vector<double> & lhs, const std::vector<double> & rhs);
double length(const std::vector<double> & vec);

} // namespace util
