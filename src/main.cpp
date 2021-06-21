#include "methods/Newton.h"
#include "methods/QuasiNewton.h"
#include "util/Function.h"
#include "util/Misc.h"
#include "util/ReplayData.h"
#include "util/VectorOps.h"
#include "util/VersionedData.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <ostream>
#include <sstream>

template <class T>
std::ostream & print(std::ostream & out, const std::vector<T> & vec, std::string delim = ", ")
{
    out << "[";
    if (!vec.empty()) {
        out << vec.front();
        std::for_each(++vec.begin(), vec.end(), [&](const T & el) { out << delim << el; });
    }
    return out << "]";
}

std::string as_string(const Function & func)
{
    std::stringstream ss;
    ss << func;
    return std::move(ss).str();
}

std::string format_for_matlab(const util::ReplayData & replay_data, const Function & func)
{
    if (func.dims() != 2) {
        return "";
    }

    std::vector<double> x, y;

    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();

    auto collect_points = util::overload(
        [&](const util::VdPoint & point) {
            min = std::min({min, point.coords[0], point.coords[1]});
            max = std::max({max, point.coords[0], point.coords[1]});

            x.push_back(point.coords[0]);
            y.push_back(point.coords[1]);
        },
        [](const auto &) {}
    );

    std::for_each(replay_data.begin(), replay_data.end(), [&](const auto & part) { part->call_func(collect_points); });

    std::stringstream ss;
    print(ss << "a = ", x) << ";\n";
    print(ss << "b = ", y) << ";\n";

    double len = max - min;
    min -= len * 0.15;
    max += len * 0.15;

    ss << "x = linspace(" << min << ", " << max << ", 100);\n";
    ss << "y = linspace(" << min << ", " << max << ", 100);\n";
    ss << "[X,Y] = meshgrid(x,y);\n";

    auto func_for_matlab = util::replace_all(as_string(func), " ^ ", ".^");
    func_for_matlab = util::replace_all(std::move(func_for_matlab), "x0", "X");
    func_for_matlab = util::replace_all(std::move(func_for_matlab), "x1", "Y");
    ss << "Z = " << func_for_matlab << ";\n";

    ss << "contour(X,Y,Z,20);hold on\n" << "plot(a,b,'-o');";

    return std::move(ss).str();
}

std::ostream & print(std::ostream & out, const util::ReplayData & replay_data)
{
    auto to_out = util::overload(
        [&out](const util::VdPoint & point) { print(out << point.version() << ": ", point.coords) << '\n'; },
        [&out](const util::VdValue & value) { out << value.version() << ": " << value.val << '\n'; },
        [&out](const util::VdComment & comm) { out << comm.version() << ": " << comm.comment << '\n'; }
    );

    for (const auto & part : replay_data) {
        part->call_func(to_out);
    }
    return out;
}

void count_and_print_newton(NewtonMethods & newtone, const Function & func, const std::vector<double> & init = {})
{
    auto print_replay = [&](const auto & res) {
        print(std::cout, newtone.replay_data()) << '\n';
        print(std::cout, res) << "\n";

        std::cout << "For matlab:\n\n" << format_for_matlab(newtone.replay_data(), func) << "\n\n";
    };

    std::cout << "Classic:\n";
    print_replay(newtone.classic(func, init));

    std::cout << "With single-dimensional search:\n";
    print_replay(newtone.with_sd_search(func, init));

    std::cout << "With descent direction:\n";
    print_replay(newtone.with_desc_dir(func, init));
}

void count_and_print_quasi(const Function & func, const std::vector<double> & init = {})
{
    QuasiNewton qn(0.000001);

    auto print_replay = [&](const auto & res) {
        print(std::cout, qn.replay_data()) << '\n';
        print(std::cout, res) << "\n";
        std::cout << "f(x) = " << qn.last_func()(res) << "\n\n";

        std::cout << "For matlab:\n\n" << format_for_matlab(qn.replay_data(), func) << "\n\n";
    };

    std::cout << "Broyder-Fletcher-Sheno:\n";
    print_replay(qn.search_bfs(func, init));

    std::cout << "Powell:\n";
    print_replay(qn.search_powell(func, init));
}

int main() {
    // NewtonMethods newtone(0.000001);

    // Function f1(2, {
    //     {{{0, 2}}, 1.},
    //     {{{1, 2}}, 1.},
    //     {{{0, 1}, {1, 1}}, 1.2}
    // });
    // std::vector init1{4., 1.};

    // std::cout << "func 1: " << f1 << "\n\n";
    // count_and_print_newton(newtone, f1, init1);

    // Function f2(2, {
    //     {{{1, 2}}, 100.},
    //     {{{0, 2}, {1, 1}}, -200.},
    //     {{{0, 4}}, 100.},

    //     {{}, 1.},
    //     {{{0, 1}}, -2},
    //     {{{0, 2}}, 1.}
    // });
    // std::vector init2{-1.2, 1.};

    // std::cout << "func 2: " << f2 << "\n\n";
    // count_and_print_newton(newtone, f2, init2);

    Function f3(2, {
        {{{1, 2}}, 100.},
        {{{0, 2}, {1, 1}}, -200.},
        {{{0, 4}}, 100.},
        {{}, 1.},
        {{{0, 1}}, -2.},
        {{{0, 2}}, 1.}
    });

    std::cout << "func 3: " << f3 << "\n\n";
    count_and_print_quasi(f3, {3., 4.});

    Function f4(2, {
        {{{0, 4}}, 1.},
        {{{0, 2}, {1, 1}}, 2.},
        {{{0, 2}}, -22.},
        {{{1, 1}}, -22.},
        {{{1, 2}}, 1.},
        {{}, 121.},

        {{{0, 2}}, 1.},
        {{{1, 4}}, 1.},
        {{}, 49.},
        {{{0, 1}, {1, 2}}, 2.},
        {{{0, 1}}, -14.},
        {{{1, 2}}, -14.}
    });

    std::cout << "func 4:" << f4 << "\n\n";
    count_and_print_quasi(f4, {0, 1});
}
