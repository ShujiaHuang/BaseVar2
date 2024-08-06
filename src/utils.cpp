#include <unistd.h>
#include "utils.h"

namespace ngslib {

    // http://c.biancheng.net/cpp/html/303.html
    bool is_readable(const char *name) {
        return (access(name, R_OK) == 0);
    }

    std::string basename(const std::string &path, const std::string delims) {
        return path.substr(path.find_last_of(delims) + 1);
    }
    
    std::string remove_filename_extension(const std::string &filename) {
        size_t p(filename.find_last_of('.'));
        return p > 0 && p != std::string::npos ? filename.substr(0, p) : filename;
    }

    std::string join(std::vector<std::string> &input, const std::string delim) {

        if(input.empty()) return "";

        std::string out_str = input[0];
        for (size_t i(1); i < input.size(); ++i) {
            out_str += (delim + input[i]);
        }

        return out_str;
    }

    std::string join(std::vector<size_t> &input, const std::string delim) {

        if(input.empty()) return "";

        std::string out_str = tostring(input[0]);
        for (size_t i(1); i < input.size(); ++i) {
            out_str += (delim + tostring(input[i]));
        }

        return out_str;
    }

    void split(std::string in_str, std::vector<std::string> &out, const char *delim, bool is_append) {

        if (!is_append) out.clear();
        while(1) {
            //erase delimiter
            int i = in_str.find_first_not_of(delim);
            if(i == std::string::npos) break;

            in_str.erase(0, i);

            i = in_str.find_first_of(delim);
            if(i == std::string::npos) {
                out.push_back(in_str);
                break;
            } else {
                std::string tok = in_str.substr(0, i);
                out.push_back(tok);
                in_str.erase(0, i);
            }
        }
    }

    void split(std::string in_str, std::vector<uint32_t> &out, const char *delim, bool is_append) {

        if (!is_append) out.clear();
        while(1) {
            //erase delimiter
            int i = in_str.find_first_not_of(delim);
            if(i == std::string::npos) break;

            in_str.erase(0, i);

            i = in_str.find_first_of(delim);
            if(i == std::string::npos) {
                out.push_back(std::stoi(in_str));
                break;
            } else {
                std::string tok = in_str.substr(0, i);
                out.push_back(std::stoi(tok));
                in_str.erase(0, i);
            }
        }
    }

}  // namespae ngslib
