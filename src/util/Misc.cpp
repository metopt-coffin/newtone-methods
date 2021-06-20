#include "util/Misc.h"

#include <string>

namespace util {

std::string replace_all(std::string inout, std::string_view what, std::string_view with)
{
    for (std::string::size_type pos{};
         inout.npos != (pos = inout.find(what.data(), pos, what.length()));
         pos += with.length()) {
        inout.replace(pos, what.length(), with.data(), with.length());
    }
    return std::move(inout);
}

} // namespace util
