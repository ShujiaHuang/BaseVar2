#include <unistd.h>
#include "utils.h"

namespace ngslib {

    // http://c.biancheng.net/cpp/html/303.html
    bool is_readable(const char *name) {
        return (access(name, R_OK) == 0);
    }

    std::string join(std::vector<std::string> &goven, const std::string delim) {

        if(goven.empty()) return "";

        std::string line = goven[0];
        for (size_t i(1); i < goven.size(); ++i) {
            line += (delim + goven[i]);
        }

        return line;
    }

    void split(std::string line, std::vector<std::string> &token, const char *delim, bool is_append) {

        if (!is_append) {
            token.clear();
        }

        while(1) {
            //erase delimiter
            int i = line.find_first_not_of(delim);
            if(i == std::string::npos) break;

            line.erase(0, i);

            i = line.find_first_of(delim);
            if(i == std::string::npos) {
                token.push_back(line);
                break;
            } else {
                std::string tok = line.substr(0, i);
                token.push_back(tok);
                line.erase(0, i);
            }
        }
    }

}  // namespae ngslib
