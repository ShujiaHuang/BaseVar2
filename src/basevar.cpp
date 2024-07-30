/*
 *  The main program of BaseVar.
 *
 *  Created on: Jul 30, 2024
 *      Author: Shujia Huang
 */
#include <iostream>
#include <time.h>


int main(int argc, char* argv[]) {
    clock_t start_time = clock();
    std::cout << "BaseVar: A software for calling variants efficiently "
              << "from low-pass whole genome sequencing data.\n" << std::endl;
    
    

    std::cerr << "** Process done. " << (double)(clock() - start_time) / CLOCKS_PER_SEC
              << "seconds elapsed **\n";
    return 0;
}