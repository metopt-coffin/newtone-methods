#pragma once
#include <algorithm>
#include <climits>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <sstream>
#include <vector>

#include "sole-solver/ProfileMatrix.h"
#include "sole-solver/QuadMatrix.h"

namespace fs = std::filesystem;
/*
 * Struct, that has only static methods-generators
 */
struct Generator
{
    /*--------------------------------QUADRATIC MATRIX--------------------------------*/

    /*
     * Generates a multidiagonal matrix with size of "dim"
     * and position of non-zero elements closer than "width" to the main diagonal
     */
    static QuadMatrix generate_quad_matrix(std::size_t dim, std::size_t width = UINT_MAX)
    {
        std::vector<std::vector<double>> m(dim, std::vector<double>(dim, 0.));
        std::uniform_real_distribution<double> dist(-limit, limit);
        for (int j = 0; j < dim; j++)
        {
            for (int i = j; i < std::min(j + width, dim); i++)
            {
                m[i][j] = dist(rand_gen);
                m[j][i] = dist(rand_gen);
            }
        }

        return QuadMatrix(std::move(m));
    }

    /*
     * Generates a quadratic matrix and writes it down to file
     */
    static void create_quad_matrix(const std::string& dir, const std::string& filename, std::size_t dim, std::size_t width = UINT_MAX)
    {
        fs::path p = ".";
        p /= dir;
        fs::create_directory(p);
        std::ofstream os(p / filename);
        os << generate_quad_matrix(dim, width);
    }

    /*
     * Reads a qudratic matrix from file
     */
    static QuadMatrix read_quad_matrix(const std::string& dir, const std::string& filename)
    {
        fs::path p = ".";
        p /= dir;
        std::ifstream is(p / filename);
        std::string line;
        std::vector<std::vector<double>> m;
        while (std::getline(is, line))
        {
            m.emplace_back(std::vector<double>());
            std::istringstream strs(line);
            double val;
            while (strs >> val)
            {
                m.back().push_back(val);
            }
        }

        return m;
    }

    /*
     * Create a test (matrix, answer vector, right part vector) and writes to given directory with given prefix
     */
    static void create_quad_test(const std::string& dir, const std::string& testname, const std::size_t& dim, const std::size_t& width = UINT_MAX)
    {
        create_test(dir, testname, std::move(generate_quad_matrix(dim, width)), std::move(generate_vector(dim)));
    }

    /*
     * Create a test (matrix, answer vector, right part vector) for the fourth task and writes to given directory with given k prefix
     */
    static void create_quad_cond_test(const std::string& dir, const int& k, const std::size_t& dim, const std::size_t& width = UINT_MAX)
    {
        create_test(dir, std::to_string(k),
                    std::move(generate_test_quad_matrix(dim, [k, width](std::vector<std::vector<double>>& matrix)
                                                           {
                                                               const std::size_t dim = matrix.size();
                                                               std::uniform_int_distribution<int> dist(-4, 0);
                                                               for (int i = 0; i < dim; i++)
                                                               {
                                                                   double sum = 0;
                                                                   for (int j = 0; j < dim; j++)
                                                                   {
                                                                       if (i != j)
                                                                       {
                                                                           sum += matrix[i][j] = dist(Generator::rand_gen);
                                                                       }
                                                                   }
                                                                   matrix[i][i] = -sum;
                                                                   if (!i)
                                                                       matrix[i][i] += std::pow(10., -k);
                                                               }
                                                           })),
                    std::move(generate_vector(dim)));
    }

    /*
     * Create a test (matrix, answer vector, right part vector) for the fourth task and writes to given directory with given k prefix
     * Note, that as we have zero-numeration we have to use changed formula: (i + j - 1) -> ((i + 1) + (j + 1) - 1) = (i + j + 1)
     */
    static void create_quad_gilbert_test(const std::string& dir, const std::size_t& dim, const std::size_t& width = UINT_MAX)
    {
        create_test(dir, std::to_string(dim),
                    std::move(generate_test_quad_matrix(dim, [width](std::vector<std::vector<double>>& matrix)
                                                           {
                                                               const std::size_t dim = matrix.size();
                                                               for (int i = 0; i < dim; i++)
                                                               {
                                                                   for (int j = 0; j < dim; j++)
                                                                   {
                                                                        matrix[i][j] = 1. / (i + j + 1.);
                                                                   }
                                                               }
                                                           })),
                    std::move(generate_vector(dim)));
    }

    /*------------------------------------PROFILE MATRIX------------------------------------*/

    /*
     * Generates a multidiagonal matrix in a profile format with size of "dim"
     * and position of non-zero elements closer than "width" to the main diagonal
     */
    static ProfileMatrix generate_profile_matrix(const std::size_t dim, const std::size_t width = UINT_MAX)
    {
        std::uniform_real_distribution<double> dist(-limit, limit);

        std::vector<double> diag(dim);
        for (std::size_t i = 0; i < dim; i++)
        { diag[i] = dist(rand_gen); }

        std::vector<std::size_t> profile(dim + 1);
        if (dim != 0)
        { profile[0] = profile[1] = 1; }
        for (std::size_t i = 2; i <= dim; i++)
        {
            profile[i] = profile[i - 1] + rand_gen() % (std::min(width - 1, i - 1) + 1);
        }

        std::vector<double> a_low(profile.back() - 1);
        std::vector<double> a_up(profile.back() - 1);
        for (std::size_t i = 0; i < a_low.size(); i++)
        {
            a_low[i] = dist(rand_gen);
            a_up[i] = dist(rand_gen);
        }

        return ProfileMatrix(std::move(diag), std::move(a_low), std::move(a_up), std::move(profile));
    }

    /*
     * Generates a profile matrix and writes it down to file
     */
    static void create_profile_matrix(const std::string& dir, const std::string& filename, std::size_t dim, std::size_t width = UINT_MAX)
    {
        fs::path p = ".";
        p /= dir;
        fs::create_directory(p);
        std::ofstream os(p / filename);
        os << std::setw(10) << generate_profile_matrix(dim, width);
    }

    /*
     * Reads a profile matrix from file
     */
    static ProfileMatrix read_profile_matrix(const std::string& dir, const std::string& filename)
    {
        fs::path p = ".";
        p /= dir;
        std::ifstream is(p / filename);

        std::vector<double> diag;
        std::vector<double> low;
        std::vector<double> up;
        std::vector<std::size_t> prof;
        read_vec(is, diag);
        read_vec(is, low);
        read_vec(is, up);
        read_vec(is, prof);

        return ProfileMatrix(std::move(diag), std::move(low), std::move(up), std::move(prof));
    }

    /*
     * Create a test (matrix, answer vector, right part vector) for the second task and writes to given directory with given k prefix
     */
    static void create_profile_cond_test(const std::string& dir, const int& k, const std::size_t& dim, const std::size_t& width = UINT_MAX)
    {
        create_test(dir, std::to_string(k),
                    std::move(generate_test_profile_matrix(dim, [k, width](std::vector<std::vector<double>>& matrix)
                                                           {
                                                               const std::size_t dim = matrix.size();
                                                               std::uniform_int_distribution<int> dist(-4, 0);
                                                               for (int i = 0; i < dim; i++)
                                                               {
                                                                   double sum = 0;
                                                                   for (int j = 0; j < dim; j++)
                                                                   {
                                                                       if (i != j && std::abs(i - j) <= width)
                                                                       {
                                                                           sum += matrix[i][j] = dist(Generator::rand_gen);
                                                                       }
                                                                   }
                                                                   matrix[i][i] = -sum;
                                                                   if (!i)
                                                                       matrix[i][i] += std::pow(10., -k);
                                                               }
                                                           })),
                    std::move(generate_vector(dim)));
    }

    /*
     * Create a test (matrix, answer vector, right part vector) for the third task and writes to given directory with given k prefix
     * Note, that as we have zero-numeration we have to use changed formula: (i + j - 1) -> ((i + 1) + (j + 1) - 1) = (i + j + 1)
     */
    static void create_profile_gilbert_test(const std::string& dir, const std::size_t& dim, const std::size_t& width = UINT_MAX)
    {
        create_test(dir, std::to_string(dim), 
                    std::move(generate_test_profile_matrix(dim, [width](std::vector<std::vector<double>>& matrix)
                                                           {
                                                               const std::size_t dim = matrix.size();
                                                               for (int i = 0; i < dim; i++)
                                                               {
                                                                   for (int j = 0; j < dim; j++)
                                                                   {
                                                                       if (std::abs(i - j) <= width)
                                                                       {
                                                                           matrix[i][j] = 1. / (i + j + 1.);
                                                                       }
                                                                   }
                                                               }
                                                           })),
                    std::move(generate_vector(dim)));
    }

    /*-----------------------------------MISCELLANEOUS-----------------------------------*/

    /*
     * Reads a vector from file
     */
    static std::vector<double> read_vector(const std::string& dir, const std::string& filename)
    {
        fs::path p = ".";
        p /= dir;
        std::ifstream is(p / filename);

        std::vector<double> read;
        read_vec(is, read);
        return read;
    }

    /*
     * Writes a vector from file
     */
    template<class T>
    static void print_vector(const std::string& dir, const std::string& filename, const std::vector<T>& v)
    {
        fs::path p = ".";
        p /= dir;
        fs::create_directory(p);
        std::ofstream os(p / filename);

        print_vec(os, v);
    }

    /*
     * Generates vector of values with given length
     */
    static std::vector<double> generate_vector(std::size_t dim)
    {
        std::uniform_real_distribution<double> dist(-limit, limit);
        std::vector<double> answer(dim);
        for (double& elem : answer)
        {
            elem = dist(rand_gen);
        }
        return answer;
    }

    static inline std::mt19937 rand_gen{std::random_device{}()};
    static constexpr double limit{ 1000. };

private:
    template <class T>
    static void read_vec(std::istream& is, std::vector<T>& vec)
    {
        std::string line;
        std::getline(is, line);
        std::istringstream strs(line);
        T val;
        while (strs >> val)
        {
            vec.push_back(val);
        }
    }

    template <class T>
    static void print_vec(std::ostream& os, const std::vector<T>& vec)
    {
        os << std::fixed << std::setprecision(20);
        for (const T& elem : vec)
        {
            os << elem << '\t';
        }
    }

    /*
     * Generate a quad matrix with the given filling rule
     */
    template <class Rule>
    static QuadMatrix generate_test_quad_matrix(const std::size_t& dim, Rule rule = Rule())
    {
        std::vector<std::vector<double>> m(dim, std::vector<double>(dim, 0.));
        rule(m);
        return QuadMatrix(m);
    }

    /*
     * Generate a profile matrix with the given filling rule
     */
    template <class Rule>
    static ProfileMatrix generate_test_profile_matrix(const std::size_t& dim, Rule rule = Rule())
    {
        std::vector<std::vector<double>> m(dim, std::vector<double>(dim, 0.));
        rule(m);
        return ProfileMatrix(QuadMatrix(m));
    }

    /*
     * Helper for creating and writing test files
     */
    template <class T>
    static void create_test(const std::string& dir, const std::string& testname, T&& matrix, std::vector<double> answer)
    {
        fs::path p = ".";
        p /= dir;
        fs::create_directory(p);

        std::string matrix_filename = testname + "_matrix.txt";
        std::string answer_filename = testname + "_answer.txt";
        std::string right_part_filename = testname + "_right.txt";
        std::vector<double> right_part = matrix * answer;

        std::ofstream os_matrix(p / matrix_filename);
        std::ofstream os_answer(p / answer_filename);
        std::ofstream os_right(p / right_part_filename);
        os_matrix << matrix;
        print_vec(os_answer, answer);
        print_vec(os_right, right_part);
    }
};
