#pragma once

#include "util/Misc.h"

#include <string>

namespace min1d {

struct Function
{
    struct Bounds
    {
        double from, to;

        double length() noexcept { return to - from; }

        double middle() noexcept { return (from + to) / 2; }
    };

public:
    Function(util::CalculateFunc calculate, Bounds bounds, std::string as_string = "");

    void reset(std::string as_string, util::CalculateFunc calculate);

    /*
     * Calculate function's value in point x
     */
    double operator()(double x) const { ++m_call_count; return m_calculate(x); }

    uint call_count() const noexcept { return m_call_count; }

    Bounds bounds() const noexcept { return m_bounds; }

    const std::string & to_string() const noexcept { return m_as_string; }

    void reset() noexcept { m_call_count = 0; }

    friend std::ostream & operator<<(std::ostream & out, const Function & func);

private:
    std::string m_as_string;
    util::CalculateFunc m_calculate;
    Bounds m_bounds;
    mutable uint m_call_count;
};

} // namespace min1d
