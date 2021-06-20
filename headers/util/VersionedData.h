#pragma once

#include "util/Misc.h"

#include <cassert>
#include <string>

namespace util {

using uint = unsigned;

enum struct VdDataKind
{
    VdPointKind,
    VdValueKind,
    VdCommentKind,
};

struct VdPoint;
struct VdValue;
struct VdComment;

struct VersionedData
{
    explicit VersionedData(uint version)
        : m_version(version)
    {}
    virtual ~VersionedData() = default;

    virtual VdDataKind get_kind() const noexcept = 0;

    template <class Func>
    auto call_func(Func && func) const
    {
#define M(type) \
    case VdDataKind::type##Kind: return func(reinterpret_cast<const type &>(*this))

        switch (get_kind()) {
            M(VdPoint);
            M(VdValue);
            M(VdComment);
        default: assert(false && "There is no such VersionedData kind");
        }

#undef M
    }

    uint version() const noexcept { return m_version; }

protected:
    uint m_version;
};

struct VdPoint : VersionedData
{
    VdPoint(uint version, std::vector<double> coords)
        : VersionedData(version)
        , coords(std::move(coords))
    {}

    VdDataKind get_kind() const noexcept override { return VdDataKind::VdPointKind; }

    std::vector<double> coords;
};

struct VdValue : VersionedData
{
    VdValue(uint version, double val)
        : VersionedData(version)
        , val(val)
    {}

    VdDataKind get_kind() const noexcept override { return VdDataKind::VdValueKind; }

    double val;
};

struct VdComment : VersionedData
{
    VdComment(uint version, std::string comment)
        : VersionedData(version)
        , comment(std::move(comment))
    {}

    VdDataKind get_kind() const noexcept override { return VdDataKind::VdCommentKind; }

    std::string comment;
};

} // namespace util
