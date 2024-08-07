#include <iostream>

#include <vector>
#include "utils.h"

int main(int argc, char *argv[]) {

    // // create thread pool with 4 worker threads
    // ThreadPool pool(4);

    // // enqueue and store future
    // auto result = pool.enqueue([](int answer) { return answer; }, 42);

    // // get result from future
    // std::cout << result.get() << std::endl;

    std::vector<std::string> a, b;
    std::vector<int> c = {1,2,3,4};
    std::vector<float> d = {0.1,2.5,3,4.0};
    a.push_back("Hello");
    a.push_back("World!");

    std::cout << ngslib::join(a, ",") << std::endl;
    ngslib::split("hello::world", b, ":");
    std::cout << ngslib::join(b, " -- ") << std::endl;
    std::cout << ngslib::join(c) << std::endl;
    std::cout << ngslib::join(c, ":") << std::endl;
    std::cout << ngslib::join(d, " - ") << std::endl;

    ngslib::split("hello::world", b, ">>", true);
    ngslib::split("hello::world", b, ">>h", true);
    std::cout << ngslib::join(b, " -- ") << std::endl;

    ngslib::split("1,2,3,4", c, ",", false);
    ngslib::split("1,2,3,4", d, ",", true);
    std::cout << ngslib::join(c, " ") << std::endl;
    std::cout << ngslib::join(d, " ") << std::endl;

    std::cout << "Get filename from Linux path: " << ngslib::basename("some/path/file.ext") << "\t" 
              << ngslib::remove_filename_extension(ngslib::basename("some/path/file.ext")) << "\n";
    std::cout << "Get filename from window path " << ngslib::basename("C:\\MyDirectory\\MyFile.bat") << "\t" 
              << ngslib::remove_filename_extension(ngslib::basename("C:\\MyDirectory\\MyFile.bat")) << "\n";

    std::cout << "is_readable: " << ngslib::is_readable("some") << "\n";
    std::cout << ngslib::safe_mkdir("a/b/c")  << "\n";
    std::cout << ngslib::safe_remove("a/b/c") << "\n";
    std::cout << ngslib::safe_remove("t_t_")  << "\n";

    std::string lmf = ngslib::get_last_modification_file("../data/");
    std::cout << "-- LMF: " << lmf << "\n";
    std::cout << "-- LMF a/b: " << ngslib::get_last_modification_file("a/b") << "\n";
    std::cout << "-- LMF ../data/bam100: " << ngslib::get_last_modification_file("../data/bam100") << "\n";

    std::vector<int> vc = ngslib::vector_slicing(c, 1, 3);  // [start, end)
    std::cout << "vector slicing `vc` size: " << vc.size() << "\t" << ngslib::join(vc, "-") << "\n";

    return 0;
}