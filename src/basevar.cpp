/*
 *  The main program of BaseVar.
 *
 *  Created on: Jul 30, 2018
 *      Author: Shujia Huang
 */
#include <iostream>
#include <getopt.h>
#include <ctime>

static const char* AUTHOR = "Shujia Huang (hshujia@qq.com)";
static const char* VERSION = "version 1.0.0";


static int usage() {
    // Usage discription
    std::cerr << "Author: " << AUTHOR << "\n";
    std::cerr << "Version (" << VERSION << ")\n";
    std::cerr << "Usage: basevar <command> [options] \n" << std::endl;


    return 1;
}


int basetype(int argc, char *argv[]) {
    // The main function for basetype
    static const struct option loptions[] = {
        {"reference", required_argument, NULL, 'R'},
        {"input", optional_argument, NULL, 'I'},
        {"align-file-list", optional_argument, NULL, 0},
        {"mapq", optional_argument, NULL, 'q'},
        {"min-af", optional_argument, NULL, 'm'},
        {"regions", optional_argument, NULL, 0},
        {"positions", optional_argument, NULL, 0},

        {"batch-count", optional_argument, NULL, 'B'},
        {"nCPU", optional_argument, NULL, 0},

        // special parameter for calculating specific population allele frequence
        {"pop-group", optional_argument, NULL, 0},
        {"--filename-has-samplename", optional_argument, NULL, 0},

        {"output-vcf", required_argument, NULL, 0},
        {"output-cvg", required_argument, NULL, 0},
        {"output-batch-file", required_argument, NULL, 0},
        {"smart-rerun", optional_argument, NULL, 0},

        {"help", no_argument, NULL, 'h'},
        {0,0,0,0}
    };

    char c;
    while((c = getopt_long(argc, argv, "I", loptions, NULL)) >= 0) {
        // 编写参数模式
        switch (c) {
            case 'I':
                /* code */
                break;
            
            default: usage(); break;
        }

    }

    return 0;
}


int main(int argc, char *argv[]) {
    clock_t start_time = clock();
    std::cerr << "BaseVar: A software for calling variants efficiently "
              << "from low-pass whole genome sequencing data.\n";    
    if (argc < 2)
        return usage();

    time_t now = time(0);
    std::cerr << "Program start on " << ctime(&now) << "\n";

    /* Coding here */

    now = time(0);
    std::cerr << "** Processing done on " << ctime(&now) << ", " 
              << (double)(clock() - start_time) / CLOCKS_PER_SEC
              << " seconds elapsed **\n" << std::endl;
    return 0;
}