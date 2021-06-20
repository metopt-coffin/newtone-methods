#pragma once

#include "util/Misc.h"
#include "util/VersionedData.h"

#include "sd_methods/Function.h"

#include <cstdint>
#include <iterator>
#include <memory>
#include <type_traits>

namespace util {

struct ReplayData
{
private:
    template <class Container>
    using ContValType = typename std::decay_t<Container>::value_type;


    template <class T>
    static constexpr bool IsVersionedData = std::is_base_of_v<VersionedData, T>;

public:
    using VdDataPtr = std::unique_ptr<const VersionedData>;
    using iterator = std::vector<VdDataPtr>::const_iterator;

    iterator begin() const { return m_data.begin(); }
    iterator end() const { return m_data.end(); }

    template <class Container>
    auto push_back(Container && cont) -> std::enable_if_t<
            is_container_v<std::decay_t<Container>> &&
            IsVersionedData<ContValType<Container>>>
    {
        if constexpr (std::is_rvalue_reference_v<Container &&>) {
            push_back_cont(Container{std::move(cont)}); // this way 'cont' will get moved from and thus cleared
        } else {
            push_back_cont(std::forward<Container>(cont));
        }
    }

    template <class Container>
    auto push_back(Container && cont) -> std::enable_if_t<
            is_container_v<std::decay_t<Container>> &&
            is_unique_ptr_v<ContValType<Container>> &&
            IsVersionedData<typename ContValType<Container>::element_type>>
    {
        static_assert(std::is_rvalue_reference_v<Container &&>, "Container of unique_ptr's should be moved");
        push_back_cont(Container{std::move(cont)}); // this way 'cont' will get moved from and thus cleared
    }

    template <class VdData>
    auto push_back(VdData && vd_data) -> std::enable_if_t<std::is_base_of_v<VersionedData, std::decay_t<VdData>>>
    {
        using ValueType = std::decay_t<VdData>;

        m_data.emplace_back(std::make_unique<const ValueType>(std::forward<VdData>(vd_data)));
        m_total_versions = std::max(m_total_versions, vd_data.version());
    }


    template <class VdDataType, class... Args>
    void emplace_back(Args &&... args) { m_data.emplace_back(std::make_unique<VdDataType>(std::forward<Args>(args)...)); }

    void clear()
    {
        m_data.clear();
        m_total_versions = 0;
    }

private:
    template <class Container>
    void push_back_cont(Container && cont)
    {
        m_data.reserve(m_data.size() + cont.size());
        for (auto & el : cont) {
            m_data.emplace_back(to_unique_ptr(to_same_ref_as(el, std::forward<Container>(cont))));
            m_total_versions = std::max(m_total_versions, el.version());
        }
    }

    template <class T>
    static auto to_unique_ptr(T && from)
    {
        if constexpr (is_unique_ptr_v<T>) {
            std::move(from);
        } else {
            return std::make_unique<std::decay_t<T>>(std::forward<T>(from));
        }
    }

private:
    std::vector<VdDataPtr> m_data;
    uint m_total_versions = 0;
};

} // namespace min1d
