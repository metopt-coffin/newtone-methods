#pragma once

#include "util/VectorOps.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iterator>
#include <ostream>
#include <memory>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

template <unsigned Depth>
struct Func
{
    using CallRes = std::vector<typename Func<Depth - 1>::CallRes>;
    using GradRes = Func<Depth + 1>;
    using SubFn = std::conditional_t<Depth == 1, std::unique_ptr<Func<0>>, Func<Depth - 1>>;
    using Underlying = std::vector<SubFn>;

    explicit Func(Underlying funcs = {})
        : m_funcs(std::move(funcs))
    {}

    virtual CallRes operator()(const std::vector<double> & x) const
    {
        CallRes res(m_funcs.size());
        std::transform(m_funcs.begin(), m_funcs.end(), res.begin(),
            [&x](const auto & func) { return deref(func)(x); });
        return res;
    }

    GradRes grad() const
    {
        std::vector<Func<Depth>> res;
        res.reserve(dims());
        for (const auto & func : m_funcs) {
            res.emplace_back(deref(func).grad());
        }
        return GradRes(std::move(res));
    }

    unsigned dims() const noexcept { return m_funcs.size(); }

    friend std::ostream & operator<<(std::ostream & out, const Func<Depth> & f) {
        out << "[\n";
        for (const auto & func : f.m_funcs) {
            out << deref(func) << '\n';
        }
        out << "]\n";
        return out;
    }
private:
    static const Func<Depth - 1> & deref(const SubFn & fn)
    {
        if constexpr (Depth == 1) {
            return *fn;
        } else {
            return fn;
        }
    }

private:
    Underlying m_funcs;
};

template<>
struct Func<0>
{
    using CallRes = double;

protected:
    Func(unsigned dims)
        : m_dims(dims)
    {}

public:
    virtual ~Func() = default;

public:
    virtual double operator()(const util::VectorT & x) const noexcept = 0;
    virtual std::unique_ptr<Func> part_der(unsigned idx) const noexcept = 0;
    virtual std::unique_ptr<Func> clone() const noexcept = 0;
    virtual std::ostream & print(std::ostream & out) const = 0;

    friend std::ostream & operator<<(std::ostream & out, const Func<0> & fn) { return fn.print(out); }

    unsigned dims() const noexcept { return m_dims; }

    Func<1> grad() const noexcept
    {
        std::vector<std::unique_ptr<Func<0>>> derivs(dims());
        for (unsigned i = 0; i < dims(); ++i) {
            derivs[i] = part_der(i);
            derivs[i]->m_dims = dims();
        }
        return Func<1>(std::move(derivs));
    }

    double as_const() const noexcept
    {
        assert(dims() == 0);
        static util::VectorT empty;
        return (*this)(empty);
    }

protected:
    unsigned m_dims;
};

using Fn = Func<0>;

std::unique_ptr<Fn> cns(double value);
std::unique_ptr<Fn> var(unsigned idx);
std::unique_ptr<Fn> add(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r);
std::unique_ptr<Fn> sub(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r);
std::unique_ptr<Fn> mul(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r);
std::unique_ptr<Fn> pow(std::unique_ptr<Fn> base, int p);

inline std::unique_ptr<Fn> operator*(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r) { return mul(std::move(l), std::move(r)); }
inline std::unique_ptr<Fn> operator*(std::unique_ptr<Fn> l, double r) { return std::move(l) * cns(r); }
inline std::unique_ptr<Fn> operator*(double l, std::unique_ptr<Fn> r) { return cns(l) * std::move(r); }

inline std::unique_ptr<Fn> operator+(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r) { return add(std::move(l), std::move(r)); }
inline std::unique_ptr<Fn> operator+(std::unique_ptr<Fn> l, double r) { return std::move(l) + cns(r); }
inline std::unique_ptr<Fn> operator+(double l, std::unique_ptr<Fn> r) { return cns(l) + std::move(r); }

inline std::unique_ptr<Fn> operator-(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r) { return sub(std::move(l), std::move(r)); }
inline std::unique_ptr<Fn> operator-(std::unique_ptr<Fn> l, double r) { return std::move(l) - cns(r); }
inline std::unique_ptr<Fn> operator-(double l, std::unique_ptr<Fn> r) { return cns(l) - std::move(r); }

inline std::unique_ptr<Fn> operator^(std::unique_ptr<Fn> l, int r) { return pow(std::move(l), r); }

inline std::unique_ptr<Fn> operator/(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r) { return mul(std::move(l), pow(std::move(r), -1)); }

struct Variable : Fn
{
    Variable(unsigned idx)
        : Fn(idx + 1)
        , index(idx)
    {}

    double operator()(const util::VectorT & x) const noexcept override { return x[index]; }
    std::unique_ptr<Fn> part_der(unsigned idx) const noexcept override;
    std::unique_ptr<Fn> clone() const noexcept override { return std::make_unique<Variable>(index); }
    std::ostream & print(std::ostream & out) const override { return out << "x" << index; }

    unsigned index;
};

struct Const : Fn
{
    Const(double val)
        : Fn(0)
        , value(val)
    {}

    double operator()(const util::VectorT & x) const noexcept override { return value; }
    std::unique_ptr<Fn> part_der(unsigned idx) const noexcept override;
    std::unique_ptr<Fn> clone() const noexcept override { return std::make_unique<Const>(value); }
    std::ostream & print(std::ostream & out) const override { return out << value; }

    double value;
};

template <class Oper>
struct BinOp : Fn, Oper
{
    BinOp(std::unique_ptr<Fn> l, std::unique_ptr<Fn> r)
        : Fn(std::max(l->dims(), r->dims()))
        , m_l(std::move(l))
        , m_r(std::move(r))
    {}

    double operator()(const util::VectorT & x) const noexcept override
    {
        double l_res = (*m_l)(x);
        double r_res = (*m_r)(x);
        return as_op()(l_res, r_res);
    }

private:
    const Oper & as_op() const noexcept { return *static_cast<const Oper *>(this); }

protected:
    std::unique_ptr<Fn> m_l;
    std::unique_ptr<Fn> m_r;
};

struct Add : BinOp<std::plus<double>>
{
    using Super = BinOp<std::plus<double>>;
    using Super::Super;

    std::unique_ptr<Fn> clone() const noexcept override { return add(m_l->clone(), m_r->clone()); }

    std::unique_ptr<Fn> part_der(unsigned idx) const noexcept override;
    std::ostream & print(std::ostream & out) const override { return out << '(' << (*m_l) << " + " << (*m_r) << ')'; }
};

struct Sub : BinOp<std::minus<double>>
{
    using Super = BinOp<std::minus<double>>;
    using Super::Super;

    std::unique_ptr<Fn> clone() const noexcept override { return sub(m_l->clone(), m_r->clone()); }

    std::unique_ptr<Fn> part_der(unsigned idx) const noexcept override;
    std::ostream & print(std::ostream & out) const override { return out << '(' << (*m_l) << " - " << (*m_r) << ')'; }
};

struct Mul : BinOp<std::multiplies<double>>
{
    using Super = BinOp<std::multiplies<double>>;
    using Super::Super;

    std::unique_ptr<Fn> clone() const noexcept override { return mul(m_l->clone(), m_r->clone()); }

    std::unique_ptr<Fn> part_der(unsigned idx) const noexcept override;
    std::ostream & print(std::ostream & out) const override { return out << '(' << (*m_l) << " * " << (*m_r) << ')'; }
};

struct Pow : Fn
{
    Pow(std::unique_ptr<Fn> base, int pow)
        : Fn(base->dims())
        , m_base(std::move(base))
        , m_pow(pow)
    {}

    std::unique_ptr<Fn> clone() const noexcept override { return pow(m_base->clone(), m_pow); }
    double operator()(const util::VectorT & x) const noexcept override { return std::pow((*m_base)(x), m_pow); }
    std::unique_ptr<Fn> part_der(unsigned idx) const noexcept override;
    std::ostream & print(std::ostream & out) const override { return out << '(' << (*m_base) << " ^ " << m_pow << ')'; }

private:
    std::unique_ptr<Fn> m_base;
    int m_pow;
};

using Function = Fn;
