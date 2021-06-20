#pragma once

#include "methods/Searcher.h"
#include "util/Function.h"
#include "util/ReplayData.h"

#include <vector>

struct NewtonMethods : Searcher
{
    using Searcher::Searcher;

    std::vector<double> classic(const Function & func, std::vector<double> init = {});
    std::vector<double> with_sd_search(const Function & func, std::vector<double> init = {});
    std::vector<double> with_desc_dir(const Function & func, std::vector<double> init = {});

};
