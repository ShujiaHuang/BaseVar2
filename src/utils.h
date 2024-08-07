// Author: Shujia Huang
// Date: 2021-08-18
#ifndef __INCLUDE_NGSLIB_UTILS_H__
#define __INCLUDE_NGSLIB_UTILS_H__

#include <sstream>
#include <string>
#include <vector>
#include <tuple>


namespace ngslib {
    // define a genome region data style
    typedef std::tuple<std::string, uint32_t, uint32_t> GenomeRegionTuple;

    // Get file name from a path
    std::string basename(const std::string &path, const std::string delims = "/\\");
    
    // remove the filename extension
    std::string remove_filename_extension(const std::string &filename);

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

    /**
     * @brief Make a folder if it doesn't exist, handling concurrent race conditions.
     * @param path  The directry path
     * 
     */
    bool safe_mkdir(std::string folder_path);
    bool safe_remove(std::string file_path);

    /**
     * @brief Get the last modification file object
     * 
     * @param directory_path 
     * @return std::string 
     * 
     */
    std::string get_last_modification_file(std::string directory_path);

    template<typename T>
    std::string join(std::vector<T> &input, const std::string delim="\t") {
        if (input.empty()) 
            return "";

        std::string out_str = tostring(input[0]);
        for (size_t i(1); i < input.size(); ++i) {
            out_str += (delim + tostring(input[i]));
        }

        return out_str;
    }

    template<typename T>
    void split(std::string in_str, std::vector<T> &out, const char *delim, bool is_append=false) {
        if (!is_append) out.clear();

        std::stringstream ss;
        T d;
        while(1) {
            //erase delimiter
            int i = in_str.find_first_not_of(delim);
            if (i == std::string::npos) break;
            in_str.erase(0, i);

            i = in_str.find_first_of(delim);
            if (i == std::string::npos) {
                ss << in_str; ss >> d;  // data type conversion
                out.push_back(d);
                break;  // hit the end of input string, break the loop
            } else {
                std::string tok = in_str.substr(0, i);
                ss << tok; ss >> d;  // data type conversion
                out.push_back(d);
                in_str.erase(0, i);
            }
            ss.clear();  // must clear the stream before next loop!!! 
        }

        return;
    }

    // Template function to slice a vector from range X to Y
    template <typename T>
    std::vector<T> vector_slicing(const std::vector<T> &v, int x, int y) {
    
        if (x > y) 
            throw std::invalid_argument(
                "[ERROR] input value error in 'std::vector<T> vector_slicing"
                "(const std::vector<T> &v, int x, int y)', start (x) larger "
                "than end(y)!");

        // Begin and End iterator (C++17), [first, last)
        auto first = x < v.size() ? v.begin()+x : v.end();
        auto last  = y < v.size() ? v.begin()+y : v.end();
    
        // Copy the element
        std::vector<T> new_v(first, last);

        // Return the results
        return new_v;
    }
    

}  // namespace ngslib

#endif  // #ifndef __INCLUDE_NGSLIB_UTILS_H__
