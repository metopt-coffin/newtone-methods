#pragma once

#include "util/ReplayData.h"

#include "sd_methods/Function.h"

#include <ostream>
#include <string_view>

namespace min1d {

struct MinSearcher
{
    virtual ~MinSearcher() = default;

    double find_min(Function func)
    {
        m_last_func.emplace(std::move(func));
        return find_min_impl();
    }

    const Function & last_func() const noexcept { return *m_last_func; }

public:
    virtual std::string_view method_name() const noexcept = 0;

protected:
    virtual double find_min_impl() noexcept = 0;

protected:
    std::optional<Function> m_last_func;
};

} // namespace min1d
