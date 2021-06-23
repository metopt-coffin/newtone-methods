#pragma once

#include "util/Function.h"
#include "util/ReplayData.h"
#include "util/VectorOps.h"

#include <vector>

namespace min_nd {

struct MinSearcher
{
    util::VectorT find_min() { return find_min_impl(); }
    util::VectorT find_min(const Function & func)
    {
        set_func(func);
        return find_min_impl();
    }

    util::VectorT find_min_traced(const Function & func, util::VectorT init = {})
    {
        set_func(func);
        m_replay_data.clear();
        return find_min_traced_impl(init.empty() ? util::VectorT(func.dims()) : std::move(init));
    }

    void set_func(const Function & func) { m_last_func = &func; }
    const Function & last_func() const { return *m_last_func; }

    const util::ReplayData & replay_data() const { return m_replay_data; }

protected:
    static const uint MAX_ITER = 1000;
protected:
    virtual util::VectorT find_min_impl() = 0;
    virtual util::VectorT find_min_traced_impl(util::VectorT init) = 0;

protected:
    util::ReplayData m_replay_data;
    const Function * m_last_func;
};

} // namespace min_nd
