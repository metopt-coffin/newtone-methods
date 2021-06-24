#pragma once

#include "methods/Searcher.h"
#include "util/Function.h"
#include "util/ReplayData.h"

#include <vector>

/*
 * Standard newton methods' implementation
 */
struct NewtonMethods : Searcher
{
    using Searcher::Searcher;

    // Classic newton method
    std::vector<double> classic(const Function & func, std::vector<double> init = {});
    // Newton method with using one-dimensional search
    std::vector<double> with_sd_search(const Function & func, std::vector<double> init = {});
    // Newton method with descent direction choice
    std::vector<double> with_desc_dir(const Function & func, std::vector<double> init = {});

};
