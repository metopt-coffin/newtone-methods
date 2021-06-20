#pragma once

#include <functional>
#include <memory>
#include <string_view>
#include <type_traits>

namespace util {

using CalculateFunc = std::function<double(double)>;

template <class Container>
struct is_container
{
    template <class Cont>
    static auto check(const Cont * c) -> std::integral_constant<bool, std::is_same_v<decltype(c->begin()), decltype(c->begin())>>;
    static auto check(...) -> std::false_type;

    using ContPtr = const Container *;
    static constexpr bool value = decltype(check(ContPtr{}))::value;
};
template <class Container>
inline constexpr bool is_container_v = is_container<Container>::value;

template <class Container, class T>
struct is_container_of : std::integral_constant<bool, std::is_same_v<typename Container::value_type, T>>
{};
template <class Container, class T>
inline constexpr bool is_container_of_v = is_container_of<Container, T>::value;

template <class T>
struct is_unique_ptr : std::false_type
{};
template <class T>
struct is_unique_ptr<std::unique_ptr<T>> : std::true_type
{};
template <class T>
static constexpr bool is_unique_ptr_v = is_unique_ptr<T>::value;

namespace detail {
template <class RefFrom, class RefTo> // RefFrom - ref, from which we get ref category, RefTo - to which we want to assign it
struct to_same_ref_as_t_impl
{
    using Refless = std::remove_reference_t<RefTo>; // here we save 'const', btw, so 'const T &&' is possible
    using type = std::conditional_t<std::is_rvalue_reference_v<RefFrom>, Refless &&, Refless &>;
};

} // namespace detail

template <class RefFrom, class RefTo>
using to_same_ref_as_t = typename detail::to_same_ref_as_t_impl<RefFrom, RefTo>::type;

/*
 * First is one, to which we want to assign new type, so it's RefTo
 */
template <class RefFrom, class RefTo>
decltype(auto) to_same_ref_as(RefTo && casted, RefFrom &&) { return static_cast<to_same_ref_as_t<RefFrom, RefTo>>(casted); }

template <class... Funcs>
struct Overloaded : public Funcs...
{
    using Funcs::operator() ...;

    template <class... Fs>
    constexpr Overloaded(Fs && ... fs)
        : Funcs(std::forward<Fs>(fs))...
    {}
};

template <class... Funcs>
constexpr auto overload(Funcs &&... funcs) { return Overloaded<std::decay_t<Funcs>...>(std::forward<Funcs>(funcs)...); }

template <class... Types>
struct TypeList;

namespace detail {

template <template <class... DestTypes> class Getter, class List>
struct FromTypeListImpl;

template <template <class... DestTypes> class Getter, class... SrcTypes>
struct FromTypeListImpl<Getter, TypeList<SrcTypes...>>
{
    using type = Getter<SrcTypes...>;
};
}

template <template <class... DestTypes> class Getter, class List>
using FromTypeList = typename detail::FromTypeListImpl<Getter, List>::type;


std::string replace_all(std::string inout, std::string_view what, std::string_view with);
} // namespace util
