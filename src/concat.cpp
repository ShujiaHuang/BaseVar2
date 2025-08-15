/**
 * @file concat.h
 * @brief Concat the BaseVar VCF files.
 * 
 * @author Shujia Huang
 * @date 2018-09-14
 * 
 */
#include "concat.h"

int _concat_basevar_outfile(const std::vector<std::string> &infiles, const std::string outfile) {
    if (infiles.empty()) return 1;
    // Get header information from the first file would be enough.
    ngslib::BGZFile f(infiles[0], "r");
    std::vector<std::string> h;
    std::string line;
    while (f.readline(line)) {
        if (line[0] != '#') break; // head char is '#'
        h.push_back(line);
    }
    f.close();
    merge_file_by_line(infiles, outfile, ngslib::join(h, "\n"), false);

    return 0;
}

int concat_runner(int argc, char *argv[]) {
    static const std::string CONCAT_USAGE = 
        "About: Concatenate or combine BaseVar's VCF files (You may use `bcftools concat` as an alternative).\n" 
        "Usage: basevar concat [options] <-O output VCFfile> [-L vcf.list] in1.vcf.gz [in2.vcf.gz ...]\n\n"
        "CAUTION: This function does not sort the postions, user should take care the concat order by themself.\n\n"
         
        "Required arguments:\n"
        "  -o, --output=FILE      Write output to a file.\n\n"

        "optional arguments:\n" 
        "  -L, --file-list=FILE   Input VCF files list, one file per row.\n"
        "  -h, --help             Show this help message and exit.";

    if (argc < 2) {
        std::cout << CONCAT_USAGE << "\n" << std::endl;
        exit(1);
    }

    std::vector<std::string> input_files;
    std::string in_filelist, output_file;

    // Parsing the commandline options. 
    static const struct option CONCAT_CMDLINE_LOPTS[] = {
        {"align-file-list", optional_argument, NULL, 'L'},
        {"output",          required_argument, NULL, 'o'},
        {"help",                  no_argument, NULL, 'h'},

        // must set this value
        {0, 0, 0, 0}
    };

    char c;
    while((c = getopt_long(argc, argv, "L:o:h", CONCAT_CMDLINE_LOPTS, NULL)) >= 0) {
        // 字符流解决命令行参数转浮点等类型的问题
        std::stringstream ss(optarg ? optarg: "");  
        switch (c) {
            case 'L': input_files = ngslib::get_firstcolumn_from_file(optarg); break;
            case 'o': output_file = optarg;                                    break;
            case 'h': 
                std::cout << CONCAT_USAGE << std::endl; 
                exit(EXIT_SUCCESS);
            default: 
                std::cerr << "Unknown argument: " << c << std::endl; 
                exit(EXIT_FAILURE);
        }
    }

    // Collect input VCF files
    while (optind < argc) {
        input_files.push_back(argv[optind++]);
    }

    /* Make sure we set valid arguments */
    if (input_files.empty() && in_filelist.empty())
        throw std::invalid_argument("[ERROR] Missing required VCF files.");

    if (output_file.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-O/--output'");
    
    std::cout << "[INFO] Finish loading arguments and we have " << input_files.size()
              << " files to concat.\n" << std::endl;
    
    return _concat_basevar_outfile(input_files, output_file);
}
