#include <iostream>

#include <vector>
#include "utils.h"

int main(int argc, char *argv[]) {

    std::vector<std::string> a, b;
    a.push_back("Hello");
    a.push_back("World!");
    std::string s = "Hello world. There are two needles in this haystack with needles.";
    std::cout << "Find first 'are': " << s.find("are") << " - " << (s.find("are", 100) ==std::string::npos) << "\n";
    std::cout << s.size() << " - " << s.length() << "\n";
    
    ngslib::split(".I am hello\tworld,w", b, ".");
    std::cout << "'I am thello\\tworld,w' : " << b.size() << " : " << ngslib::join(b, "--") << std::endl;

    ngslib::split("hello,world::ww", b, ":");
    std::cout << "'hello,world::ww'    : " << b.size() << " " << ngslib::join(b, "--") << std::endl;

    ngslib::split(",hello,world,,ww,", b, ",");
    std::cout << "',hello,world,,ww,'  : " << b.size() << " " << ngslib::join(b, "--") << std::endl;

    ngslib::split("hello::world", b, ">>", false);
    ngslib::split("hello::world", b, "h", false);
    std::cout << b[0] << ": " << ngslib::join(b, "--") << std::endl;

    std::vector<int> c = {1,2,3,4};
    std::cout << ngslib::join(c) << std::endl;
    std::cout << ngslib::join(c, ":") << std::endl;

    std::vector<float> d = {0.1,2.5,3,4.0};
    std::cout << ngslib::join(d, " - ") << std::endl;

    ngslib::split("1,2,3,4", c, ",", false);
    ngslib::split(",1,2,3,4", d, ",", false);
    std::cout << ngslib::join(c, " ") << std::endl;
    std::cout << ngslib::join(d, " ") << std::endl;

    std::cout << "Get filename from Linux path: " << ngslib::basename("some/path/file.ext") << "\t" 
              << ngslib::remove_filename_extension(ngslib::basename("some/path/file.ext")) << "\n";
    std::cout << "Get filename from window path " << ngslib::basename("C:\\\\MyDirectory\\\\MyFile.bat") << "\t" 
              << ngslib::remove_filename_extension("C:\\\\MyDirectory\\\\MyFile.bat") << "\n";

    std::cout << "Get dirname from Linux path: " << ngslib::dirname("some/path/file.ext") << "\t" 
              << ngslib::dirname("./a") << "\n";
    std::cout << "Get dirname from window path " << ngslib::dirname("C:\\\\MyDirectory\\\\MyFile.bat") << "\n";

    std::cout << "is_readable: " << ngslib::is_readable("some") << " - test_threadpool - " << ngslib::is_readable("test_bamrecord") << "\n";
    std::cout << ngslib::safe_mkdir("a/b/c")  << "\n";
    std::cout << ngslib::safe_remove("a/b/c") << "\n";
    std::cout << ngslib::safe_remove("t_t_")  << "\n";

    std::string lmf = ngslib::get_last_modification_file("../data/");
    std::cout << "-- LMF: " << lmf << "\n";
    std::cout << "-- LMF a/b: " << ngslib::get_last_modification_file("a/b") << "\n";
    std::cout << "-- LMF ../data/bam100: " << ngslib::get_last_modification_file("../data/bam100") << "\n";

    std::vector<int> vc = ngslib::vector_slicing(c, 1, 3);  // [start, end)
    std::cout << "vector slicing `vc` size: " << vc.size() << "\t" << ngslib::join(vc, "-") << "\n";

    std::cout << "Get Linux absolute path: " << ngslib::abspath("some/path/file.ext") << "\n";
    std::cout << "Get linux absolute path: " << ngslib::abspath("/file.ext") << "\n";

    return 0;
}