#pragma once

#include "util/Function.h"
#include "util/ReplayData.h"

#include <vector>

struct NewtonMethods
{
    NewtonMethods(double eps)
        : m_eps(eps)
    {}

    std::vector<double> classic(const Function & func, std::vector<double> init = {});
    std::vector<double> with_sd_search(const Function & func, std::vector<double> init = {});
    std::vector<double> with_desc_dir(const Function & func, std::vector<double> init = {});

    const Function & last_func() const noexcept { return *m_last_func; }
    const util::ReplayData & replay_data() const noexcept { return m_replay_data; }
private:
    std::vector<double> init_method(const Function & func, std::vector<double> init);

    void log_x(unsigned iter_num, const std::vector<double> & x);
    void log_alpha(unsigned iter_num, double alpha);

protected:
    static constexpr unsigned MaxIter = 3000;

    double m_eps;

    const Function * m_last_func = nullptr;
    util::ReplayData m_replay_data;
};
