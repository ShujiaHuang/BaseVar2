// Author: Shujia Huang
// Date: 2021-08-18
#ifndef __INCLUDE_NGSLIB_UTILS_H__
#define __INCLUDE_NGSLIB_UTILS_H__

#include <vector>
#include <tuple>
#include <string>
#include <sstream>


namespace ngslib {
    // define a genome region data style
    typedef std::tuple<std::string, uint32_t, uint32_t> GenomeRegionTuple;

    // Template function can only be defined in C++ header file
    template<typename T>
    std::string tostring(T d) {
        std::stringstream ss;
        ss << d;
        return ss.str();
    }

    /** Check if a file is readable and exists.
     * 
     * @param name Name of a file to test
     * @return a bool type for file is readable and exists or not.
     * 
     */
    bool is_readable(const char *name);
    inline bool is_readable(const std::string &name) { return is_readable(name.c_str()); }

    std::string join(std::vector<std::string> &input, const std::string delim="\t");
    void split(std::string in_str, std::vector<std::string> &out, const char *delim, bool is_append=false);
    void split(std::string in_str, std::vector<uint32_t> &out, const char *delim, bool is_append=false);

}  // namespace ngslib

#endif  // #ifndef __INCLUDE_NGSLIB_UTILS_H__
