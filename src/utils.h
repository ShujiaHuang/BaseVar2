// Author: Shujia Huang
// Date: 2021-08-18
#ifndef __INCLUDE_NGSLIB_UTILS_H__
#define __INCLUDE_NGSLIB_UTILS_H__

#include <string>
#include <sstream>

namespace ngslib {

    // Template function can only be defined in C++ header file
    template<typename T>
    std::string tostring(T d) {
        std::stringstream ss;
        ss << d;
        return ss.str();
    }

    /** Check if a file is readable and exists.
     * @param name Name of a file to test
     * @return a bool type for file is readable and exists or not.
     */
    bool is_readable(const char *name);
    inline bool is_readable(const std::string &name) { return is_readable(name.c_str()); }

}  // namespace ngslib

#endif  // #ifndef __INCLUDE_NGSLIB_UTILS_H__
