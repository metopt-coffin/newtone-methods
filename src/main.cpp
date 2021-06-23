#include "methods/Newton.h"
#include "methods/QuasiNewton.h"
#include "nd_methods/FastestDescent.h"
#include "sd_methods/Brent.h"
#include "util/Function.h"
#include "util/Misc.h"
#include "util/ReplayData.h"
#include "util/VectorOps.h"
#include "util/VersionedData.h"

#include <algorithm>
#include <iomanip>
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

void count_and_print_fast_desc(min_nd::FastestDescent & fd, const Function & func, util::VectorT init = {})
{
    auto print_replay = [&](const auto & res) {
        print(std::cout, fd.replay_data()) << '\n';
        print(std::cout, res) << "\n";
        std::cout << "f(x) = " << fd.last_func()(res) << "\n\n";

        std::cout << "For matlab:\n\n" << format_for_matlab(fd.replay_data(), func) << "\n\n";
    };

    std::cout << "Fastest descend:\n";
    print_replay(fd.find_min_traced(func, std::move(init)));
}

struct HardcodedHesse : public Func<2>
{

};

struct HardcodedFunction : public Function
{
    CallRes operator()(const util::VectorT & vec) const noexcept override
    {
        double a = (vec[0] - 1.) / 2.;
        double b = (vec[1] - 1.) / 3.;
        double c = (vec[0] - 2.) / 2.;
        double d = (vec[1] - 1.) / 3.;

        double a_b_1 = a * a + b * b + 1.;
        double c_d_1 = c * c + d * d + 1.;
        return 100. - (2. / a_b_1) - (1. / c_d_1);
    }
};

auto get_f4_p()
{
    auto a = ((var(0) - 1) * 0.5) ^ 2;
    auto b = ((var(1) - 1) * (1. / 3.)) ^ 2;
    auto c = ((var(0) - 2) * 0.5) ^ 2;

    auto a_b_1 = std::move(a) + b->clone() + 1.;
    auto b_c_1 = std::move(b) + std::move(c) + 1.;

    return 100. - (cns(2.) / std::move(a_b_1)) - (cns(1.) / std::move(b_c_1));
}

int main() {
    

    std::cout << std::setprecision(8);

    NewtonMethods newtone(0.000001);

    auto f1_p = (100. * ((var(1) - (var(0)) ^ 2) ^ 2)) + ((1 - var(0)) ^ 2);
    Function & f1 = *f1_p;

    auto f2_p = (((var(0) ^ 2) + var(1) - 11.) ^ 2) + ((var(0) + (var(1) ^ 2) - 7.) ^ 2);
    Function & f2 = *f2_p;

    auto f3_p = ((var(0) + (10. * var(1))) ^ 2) +
        (5. * ((var(2) - var(3)) ^ 2)) +
        ((var(1) - (2 * var(2))) ^ 4) +
        (10 * ((var(0) - var(1)) ^ 4));

    Function & f3 = *f3_p;

    auto f4_p = get_f4_p();
    Function & f4 = *f4_p;

    util::VectorT init{1., 2.};

    // std::cout << "Z = 100 * (Y - X.^2).^2 + (1 - X).^2\n";
    std::cout << f4 << "\n";
    // count_and_print_newton(newtone, f2, {0.8, 0.8});
    count_and_print_newton(newtone, f4, init);
    count_and_print_quasi(f4, init);

    // min_nd::FastestDescent fast_d(0.000001);
    // fast_d.find_min_traced(f1);
    // count_and_print_fast_desc(fast_d, f2, {-1.2, 1.});
}
