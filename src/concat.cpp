/**
 * @file concat.cpp
 * @brief Concatenate or combine BaseVar VCF files.
 *
 * Supports two modes:
 *   1. Default mode: line-by-line decompression/recompression (safe, with validation).
 *   2. --naive mode: BGZF block-level raw concatenation without recompression (extremely fast).
 *
 * The --naive mode is ideal when all input VCFs are produced by the same basevar run
 * (same header, non-overlapping regions, sorted order).
 *
 * @author Shujia Huang
 * @date 2018-09-14
 *
 */

#include "concat.h"
#include "io/iobgzf.h"
#include "io/bgzf_concat.h"  // bgzf_copy_blocks_to
#include "io/utils.h"      // ngslib::get_firstcolumn_from_file

#include <stdexcept>
#include <vector>
#include <string>

// ---------------------------------------------------------------------------
//  Naive concat: BGZF block-level concatenation (no recompression)
// ---------------------------------------------------------------------------
//
// The core block-level copy engine and skip_vcf_gz_header() are in
// io/bgzf_concat.h (shared with variant_caller.cpp).
//

/**
 * @brief Check that all input VCF files have compatible headers.
 *
 * For basevar's use case all files are produced by the same run, so we
 * verify:
 *   - Same number of samples
 *   - Same sample names in the same order
 */
static void naive_concat_check_headers(const std::vector<std::string> &infiles) {
    std::cerr << "[naive_concat] Checking the headers of " << infiles.size() << " files.\n";

    // Read header lines (all '#' lines) from the first file as reference
    ngslib::BGZFile ref(infiles[0], "r");
    std::vector<std::string> ref_header;
    std::vector<std::string> ref_samples;
    std::string line;
    while (ref.readline(line)) {
        if (line[0] != '#') break;
        ref_header.push_back(line);
        // Extract sample names from #CHROM line
        if (line.size() > 6 && line[1] == 'C' && line[2] == 'H' &&
            line[3] == 'R' && line[4] == 'O' && line[5] == 'M') {
            // Split by tab, samples start at column 9
            size_t col = 0, start = 0;
            for (size_t i = 0; i <= line.size(); ++i) {
                if (i == line.size() || line[i] == '\t') {
                    if (col >= 9) ref_samples.push_back(line.substr(start, i - start));
                    col++;
                    start = i + 1;
                }
            }
        }
    }
    ref.close();

    if (ref_header.empty())
        throw std::runtime_error("[naive_concat] Error: no header found in " + infiles[0]);

    // Check remaining files
    for (size_t i = 1; i < infiles.size(); ++i) {
        ngslib::BGZFile f(infiles[i], "r");
        std::vector<std::string> samples;
        std::string line;
        while (f.readline(line)) {
            if (line[0] != '#') break;
            if (line.size() > 6 && line[1] == 'C' && line[2] == 'H' &&
                line[3] == 'R' && line[4] == 'O' && line[5] == 'M') {
                size_t col = 0, start = 0;
                for (size_t j = 0; j <= line.size(); ++j) {
                    if (j == line.size() || line[j] == '\t') {
                        if (col >= 9) samples.push_back(line.substr(start, j - start));
                        col++;
                        start = j + 1;
                    }
                }
            }
        }
        f.close();

        if (samples.size() != ref_samples.size())
            throw std::runtime_error(
                "[naive_concat] Error: cannot concatenate, different number of samples: "
                + std::to_string(ref_samples.size()) + " vs " + std::to_string(samples.size())
                + " in " + infiles[0] + " vs " + infiles[i]);

        for (size_t j = 0; j < ref_samples.size(); ++j) {
            if (ref_samples[j] != samples[j])
                throw std::runtime_error(
                    "[naive_concat] Error: cannot concatenate, different sample names at index "
                    + std::to_string(j) + ": '" + ref_samples[j] + "' vs '" + samples[j]
                    + "' in " + infiles[0] + " vs " + infiles[i]);
        }
    }
    std::cerr << "[naive_concat] Done, the headers are compatible.\n";
}

/**
 * @brief Perform naive BGZF block-level concatenation.
 *
 * Instead of decompressing and recompressing each line, this function copies
 * raw BGZF blocks directly from input to output, stripping intermediate EOF
 * blocks.  This is orders of magnitude faster than line-by-line concat.
 *
 * @param infiles     Vector of input .vcf.gz file paths.
 * @param outfile     Output .vcf.gz file path.
 * @param trust_headers  If true, skip header compatibility check (--naive-force).
 */
static void naive_concat(const std::vector<std::string> &infiles,
                         const std::string &outfile,
                         bool trust_headers)
{
    if (infiles.empty()) return;

    // Pre-flight: verify all inputs are BGZF-compressed VCF (not BCF, not plain text)
    for (const auto &fn : infiles) {
        ngslib::BGZFile tmp(fn, "r");
        if (!tmp.is_open()) {
            std::cerr << "[naive_concat] Error: cannot open input file: " << fn << std::endl;
            exit(1);
        }
        if (!tmp.is_compressed()) {
            std::cerr << "[naive_concat] Error: input file is not BGZF-compressed: " << fn
                      << "\n  The --naive mode requires BGZF-compressed (.vcf.gz) input files."
                      << "\n  Use the default mode (without --naive) for plain-text VCF files."
                      << std::endl;
            exit(1);
        }
        // Check that the decompressed content is VCF text (starts with '#'), not BCF binary
        if (tmp.read_block() != 0 || !tmp.block_length()) {
            std::cerr << "[naive_concat] Error: failed to read first block of " << fn << std::endl;
            exit(1);
        }
        if (tmp.block_length() > 0 && tmp.uncompressed_data()[0] != '#') {
            std::cerr << "[naive_concat] Error: input file is not a VCF text file: " << fn
                      << "\n  The --naive mode only supports BGZF-compressed VCF (.vcf.gz) files."
                      << "\n  BCF binary format is not supported."
                      << "\n  Convert BCF to VCF.gz first: bcftools view -Oz input.bcf -o output.vcf.gz"
                      << std::endl;
            exit(1);
        }
    }

    if (!trust_headers)
        naive_concat_check_headers(infiles);

    // Verify output file suffix: naive mode always produces BGZF-compressed output
    if (ngslib::suffix_name(outfile) != ".gz") {
        std::cerr << "[naive_concat] Error: output file must have a .gz suffix for --naive/--naive-force mode: " << outfile
                  << "\n  The --naive/--naive-force mode produces BGZF-compressed output and requires a .gz suffix."
                  << "\n  Use a .vcf.gz (or .gz) output file, or switch to the default mode for plain-text VCF."
                  << std::endl;
        exit(1);
    }

    ngslib::BGZFile bgzf_out(outfile, "w");
    if (!bgzf_out.is_open())
        throw std::runtime_error("[naive_concat] Error: cannot open output file: " + outfile);

    // Core block-level merge: header-aware BGZF block copying (shared engine from bgzf_concat.h)
    ngslib::bgzf_naive_concat(infiles, bgzf_out, /*remove_inputs=*/false);
    // bgzf_out is closed here by ~BGZFile; bgzf_close writes the final EOF block.
}

// ---------------------------------------------------------------------------
//  Default (line-by-line) concat
// ---------------------------------------------------------------------------

/**
 * @brief Merge multiple VCF files line-by-line with a unified header.
 *
 * Supports both BGZF-compressed (.gz) and plain text output.
 * Skips '#' header lines from all input files and writes the provided header.
 * Optionally removes input files after merging.
 */
static void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile,
                               std::string header, bool is_remove_tempfile)
{
    if (infiles.empty()) return;

    bool is_compress = (ngslib::suffix_name(outfile) == ".gz") ? true : false;
    ngslib::BGZFile OUT(outfile, is_compress ? "wb" : "uw");
    OUT << header << "\n";

    for (const auto &fn : infiles) {
        ngslib::BGZFile f(fn, "r");
        std::string line;

        while (f.readline(line)) {
            if (line[0] == '#') continue;
            OUT << line << "\n";
        }
        OUT.flush();

        if (is_remove_tempfile) ngslib::safe_remove(fn);
    }

    OUT.close();
}

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

// ---------------------------------------------------------------------------
//  CLI runner
// ---------------------------------------------------------------------------
static int concat_runner_impl(int argc, char *argv[]) {
    static const std::string CONCAT_USAGE = 
        "About: Concatenate or combine BaseVar's VCF files.\n"
        "Usage: basevar concat [options] -o <output.vcf.gz> in1.vcf.gz [in2.vcf.gz ...]\n\n"
        "CAUTION: This function does not sort the positions, user should take care the concat order by themself.\n\n"
        "  With --naive, files are concatenated at the BGZF block level without recompression,\n"
        "  which is extremely fast.  All input files must have the same header (same samples).\n\n"
         
        "Required arguments:\n"
        "  -o, --output=FILE      Write output to a file. The output format is determined by the file\n"
        "                         suffix: '.vcf.gz' or '.gz' for BGZF-compressed VCF, '.vcf' for plain\n"
        "                         text VCF. Note: --naive mode always produces BGZF-compressed output.\n\n"

        "Optional arguments:\n" 
        "  -L, --file-list=FILE   Input VCF files list, one file per row.\n"
        "  -n, --naive            Concatenate without recompression (BGZF block-level, very fast).\n"
        "      --naive-force      Same as --naive but skip header compatibility check.\n"
        "  -h, --help             Show this help message and exit.";

    if (argc < 2) {
        std::cout << CONCAT_USAGE << "\n" << std::endl;
        exit(1);
    }

    std::vector<std::string> input_files;
    std::string output_file;
    bool naive_mode = false;
    bool naive_force = false;

    // Parsing the commandline options. 
    static const struct option CONCAT_CMDLINE_LOPTS[] = {
        {"file-list",   required_argument, NULL, 'L'},
        {"output",      required_argument, NULL, 'o'},
        {"naive",             no_argument, NULL, 'n'},
        {"naive-force",       no_argument, NULL,  1 },  // long-only option
        {"help",              no_argument, NULL, 'h'},
        {0, 0, 0, 0}
    };

    char c;
    while((c = getopt_long(argc, argv, "L:o:nh", CONCAT_CMDLINE_LOPTS, NULL)) >= 0) {
        switch (c) {
            case 'L': input_files = ngslib::get_firstcolumn_from_file(optarg); break;
            case 'o': output_file = optarg;                                    break;
            case 'n': naive_mode = true;                                       break;
            case  1 : naive_mode = true; naive_force = true;                   break;
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
    if (input_files.empty())
        throw std::invalid_argument("[ERROR] Missing required VCF files.");

    if (output_file.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-o/--output'");
    
    std::cout << "[INFO] Finish loading arguments and we have " << input_files.size()
              << " files to concat." << std::endl;

    if (naive_mode) {
        std::cout << "[INFO] Using naive mode (BGZF block-level concat, no recompression)."
                  << (naive_force ? " Header check skipped (--naive-force)." : "")
                  << std::endl;
        naive_concat(input_files, output_file, naive_force);
    } else {
        std::cout << "[INFO] Using default mode (line-by-line concat)." << std::endl;
        return _concat_basevar_outfile(input_files, output_file);
    }

    return 0;
}

int concat_runner(int argc, char *argv[]) {
    try {
        return concat_runner_impl(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
