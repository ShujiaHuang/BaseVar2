#include <unistd.h>
#include <filesystem>  // C++17 library
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

    bool safe_mkdir(std::string folder_path) {
        // set the folder_path tobe 'filesystem::path'
        std::filesystem::path dir_path(folder_path);
        return std::filesystem::create_directories(dir_path);
    }

    bool safe_remove(std::string filepath) {
        std::filesystem::path fname(filepath);
        return std::filesystem::remove(fname);
    }

    std::string get_last_modification_file(std::string directory_path) {
        // Find the last modification file in a directory and return it.
        std::filesystem::path lmf;        
        std::filesystem::file_time_type lmf_time, cf_time;  

        bool flag(true);
        for (const auto &fn : std::filesystem::directory_iterator(directory_path)) {
            cf_time = std::filesystem::last_write_time(fn);
            if (flag) {
                flag = false;
                lmf = fn;
                lmf_time = cf_time;
            } else if (cf_time > lmf_time) {
                // find last modification file.
                lmf = fn;
                lmf_time = cf_time;
            }
        }

        return lmf.string();
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
