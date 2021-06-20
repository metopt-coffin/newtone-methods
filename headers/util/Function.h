#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <ostream>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

/**
 * Term implements idea: coeff * (x_1 ^ pow_1) * (x_2 ^ pow_2) * .. (x_dim ^ pow_dim)
 */
struct Term
{
    Term(unsigned dims, std::vector<std::pair<unsigned, int>> idx_pow, double coeff)
        : m_pows(create_pows(dims, std::move(idx_pow)))
        , m_coeff(coeff)
        , m_hash(count_hash())
    {}

private:
    Term(double coeff, std::vector<int> pows)
        : m_pows(std::move(pows))
        , m_coeff(coeff)
        , m_hash(count_hash())
    {}

public:

    double operator()(const std::vector<double> x) const noexcept;

    Term partial_der(unsigned idx) const
    {
        if (m_pows[idx] == 0) {
            return Term(0., std::vector<int>(dims(), 0));
        } else {
            auto new_pows = m_pows;
            new_pows[idx]--;
            return Term(m_coeff * m_pows[idx], std::move(new_pows));
        }
    }

    unsigned dims() const noexcept { return m_pows.size(); }
    std::size_t hash() const noexcept { return m_hash; }

    double & coeff() noexcept { return m_coeff; }
    double coeff() const noexcept { return m_coeff; }

    bool operator==(const Term & rhs) const noexcept { return m_pows == rhs.m_pows; }
    bool operator!=(const Term & rhs) const noexcept { return m_pows != rhs.m_pows; }

    friend std::ostream & operator<<(std::ostream & out, const Term & term )
    {
        out << term.coeff();
        for (unsigned i = 0; i < term.m_pows.size(); ++i) {
            if (term.m_pows[i] != 0) {
                out << " * (x" << i << " ^ " << term.m_pows[i] << ")";
            }
        }
        return out;
    }

public:
    struct Hasher;

private:
    std::size_t count_hash() const noexcept;

    static std::vector<int> create_pows(unsigned dims, std::vector<std::pair<unsigned, int>> idx_pows);

private:
    std::vector<int> m_pows;
    double m_coeff;

    std::size_t m_hash;
};

struct Term::Hasher
{
    std::size_t operator()(const Term & term) const noexcept { return term.hash(); }
};

template <unsigned Depth>
struct Func
{
    using CallRes = std::vector<typename Func<Depth - 1>::CallRes>;
    using GradRes = Func<Depth + 1>;
    using Underlying = std::vector<Func<Depth - 1>>;

    explicit Func(Underlying funcs = {})
        : m_funcs(std::move(funcs))
    {}

    CallRes operator()(const std::vector<double> & x) const
    {
        CallRes res(m_funcs.size());
        std::transform(m_funcs.begin(), m_funcs.end(), res.begin(),
            [&x](const auto & func) { return func(x); });
        return res;
    }

    GradRes grad() const
    {
        std::vector<Func<Depth>> res;
        res.reserve(dims());
        for (const auto & func : m_funcs) {
            res.push_back(func.grad());
        }
        return GradRes(std::move(res));
    }

    unsigned dims() const noexcept { return m_funcs.size(); }

    friend std::ostream & operator<<(std::ostream & out, const Func<Depth> f) {
        out << "[\n";
        for (const auto & func : f.m_funcs) {
            out << func << '\n';
        }
        out << "]\n";
        return out;
    }

private:
    Underlying m_funcs;
};

template <>
struct Func<0>
{
    using CallRes = double;

    explicit Func(unsigned dims = 0, std::vector<std::pair<std::vector<std::pair<unsigned, int>>, double>> terms = {})
        : m_dims(dims)
    {
        for (auto & [term, coeff] : terms) {
            add_term(Term(dims, std::move(term), coeff));
        }
    }

    CallRes operator()(const std::vector<double> & x) const noexcept;

    Func<0> partial_der(unsigned idx) const
    {
        Func res(dims());
        for (const auto & term : m_terms) {
            res.add_term(term.partial_der(idx));
        }
        return res;
    }

    Func<1> grad() const
    {
        std::vector<Func<0>> res(dims());
        for (unsigned i = 0; i < dims(); ++i) {
            res[i] = partial_der(i);
        }
        return Func<1>(std::move(res));
    }

    unsigned dims() const noexcept { return m_dims; }

    friend std::ostream & operator<<(std::ostream & out, const Func<0> & func)
    {
        if (!func.m_terms.empty()) {
            out << '(' << *func.m_terms.begin() << ')';
            std::for_each(++func.m_terms.begin(), func.m_terms.end(), [&out](const Term & term) { out << " + (" << term << ')'; });
        }
        return out;
    }

private:
    using Terms = std::unordered_set<Term, Term::Hasher>;

    void add_term(Term term)
    {
        auto [it, emplaced] = m_terms.emplace(std::move(term));
        if (!emplaced) {
            term_mut(it).coeff() += term.coeff();
        }
    }

    static Term & term_mut(const Terms::iterator & it) noexcept { return const_cast<Term &>(*it); }
private:
    Terms m_terms;
    unsigned m_dims;
};

using Function = Func<0>;
