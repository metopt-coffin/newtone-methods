#include "methods/Searcher.h"


auto Searcher::init_method(const Function & func, std::vector<double> init) -> PointT
{
    m_last_func = &func;
    m_replay_data.clear();

    if (init.empty()) {
        return std::vector<double>(func.dims(), 0.);
    } else {
        return std::move(init);
    }
}

void Searcher::log_x(unsigned iter_num, const std::vector<double> & x)
{
    m_replay_data.emplace_back<util::VdComment>(iter_num, "x:");
    m_replay_data.emplace_back<util::VdPoint>(iter_num, x);
}

void Searcher::log_alpha(unsigned iter_num, double alpha)
{
    m_replay_data.emplace_back<util::VdComment>(iter_num, "alpha:");
    m_replay_data.emplace_back<util::VdValue>(iter_num, alpha);
}
