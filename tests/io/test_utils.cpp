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
    a.push_back("Hello");
    a.push_back("World!");

    std::cout << ngslib::join(a, ",") << std::endl;
    ngslib::split("hello::world", b, ":");
    std::cout << ngslib::join(b, " -- ") << std::endl;

    ngslib::split("hello::world", b, ">>", true);
    ngslib::split("hello::world", b, ">>h", true);
    std::cout << ngslib::join(b, " -- ") << std::endl;

    std::cout << "Get filename from Linux path: " << ngslib::basename("some/path/file.ext") << "\t" 
              << ngslib::remove_filename_extension(ngslib::basename("some/path/file.ext")) << "\n";
    std::cout << "Get filename from window path " << ngslib::basename("C:\\MyDirectory\\MyFile.bat") << "\t" 
              << ngslib::remove_filename_extension(ngslib::basename("C:\\MyDirectory\\MyFile.bat")) << "\n";
    
    std::cout << "is_readable: " << ngslib::is_readable("some") << "\n";

    ngslib::safe_mkdir("a/b/");
    ngslib::safe_remove("a/b/c");
    ngslib::safe_remove("test_fasta");

    return 0;

}