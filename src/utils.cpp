#include <unistd.h>

#include "utils.h"

namespace ngslib {

    // http://c.biancheng.net/cpp/html/303.html
    bool is_readable(const char *name) {
        return (access(name, R_OK) == 0);
    }

}  // namespae ngslib
