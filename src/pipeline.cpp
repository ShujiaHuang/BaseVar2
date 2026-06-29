/**
 * @file pipeline.cpp
 *
 * @brief Implementation of the `basevar pipeline` subcommand.
 *
 *  Functionally equivalent to scripts/create_pipeline.py.
 *
 *  @author Shujia Huang
 *  @date   2026-05-19
 */
#include "pipeline.h"

#include <getopt.h>
#include <unistd.h>
#include <climits>

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <filesystem>

namespace basevar {
namespace pipeline {

// ---------------------------------------------------------------------------
// .fai loading
// ---------------------------------------------------------------------------
std::vector<ngslib::GenomeRegion> load_reference_fai(
    const std::string& fai_path,
    const std::vector<std::string>& chroms_filter)
{
    std::ifstream fh(fai_path);
    if (!fh) {
        throw std::runtime_error("[ERROR] Cannot open .fai file: " + fai_path);
    }

    // Build a set for O(log n) membership tests when a whitelist is given.
    std::vector<std::string> filter_sorted(chroms_filter);
    std::sort(filter_sorted.begin(), filter_sorted.end());
    const bool has_filter = !filter_sorted.empty();

    std::vector<ngslib::GenomeRegion> out;
    std::string line;
    size_t lineno = 0;
    while (std::getline(fh, line)) {
        ++lineno;
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string chrom;
        uint64_t length = 0;
        if (!(iss >> chrom >> length)) {
            throw std::runtime_error(
                "[ERROR] Malformed .fai line " + std::to_string(lineno) +
                " in " + fai_path + ": '" + line + "'");
        }

        if (has_filter &&
            !std::binary_search(filter_sorted.begin(), filter_sorted.end(), chrom)) {
            continue;
        }
        out.emplace_back(chrom, 1, static_cast<uint32_t>(length));
    }

    return out;
}

// ---------------------------------------------------------------------------
// region string parsing
// ---------------------------------------------------------------------------
ngslib::GenomeRegion parse_region_string(
    const std::string& region_str_in,
    const std::map<std::string, uint32_t>& ref_fai_map)
{
    // trim leading/trailing whitespace
    std::string region_str = region_str_in;
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    region_str.erase(region_str.begin(),
                     std::find_if(region_str.begin(), region_str.end(), not_space));
    region_str.erase(std::find_if(region_str.rbegin(), region_str.rend(), not_space).base(),
                     region_str.end());

    if (region_str.empty()) {
        throw std::invalid_argument("[ERROR] Empty region string.");
    }

    // Split chrom and coords on the FIRST ':' only.
    std::string chrom, coord_str;
    size_t colon = region_str.find(':');
    if (colon == std::string::npos) {
        chrom = region_str;
    } else {
        chrom     = region_str.substr(0, colon);
        coord_str = region_str.substr(colon + 1);
    }

    auto it = ref_fai_map.find(chrom);
    if (it == ref_fai_map.end()) {
        throw std::invalid_argument("[ERROR] Chromosome '" + chrom +
                                    "' not found in reference fai.");
    }
    uint32_t chrom_len = it->second;

    if (coord_str.empty()) {
        // 'chr'
        return ngslib::GenomeRegion(chrom, 1, chrom_len);
    }

    // Parse 'start' or 'start-end'.
    size_t dash = coord_str.find('-');
    std::string start_s = (dash == std::string::npos) ? coord_str
                                                      : coord_str.substr(0, dash);
    std::string end_s   = (dash == std::string::npos) ? ""
                                                      : coord_str.substr(dash + 1);

    long long start_ll = 0, end_ll = 0;
    try {
        start_ll = std::stoll(start_s);
    } catch (const std::exception&) {
        throw std::invalid_argument(
            "[ERROR] Invalid start position in region '" + region_str + "'.");
    }

    if (end_s.empty()) {
        end_ll = chrom_len;
    } else {
        try {
            end_ll = std::stoll(end_s);
        } catch (const std::exception&) {
            throw std::invalid_argument(
                "[ERROR] Invalid end position in region '" + region_str + "'.");
        }
    }

    if (start_ll < 1) {
        throw std::invalid_argument(
            "[ERROR] Start position must be >= 1 in region '" + region_str + "'.");
    }
    if (start_ll > end_ll) {
        throw std::invalid_argument(
            "[ERROR] Start position > end position in region '" + region_str + "'.");
    }
    if (static_cast<uint64_t>(end_ll) > chrom_len) {
        throw std::invalid_argument(
            "[ERROR] End position " + std::to_string(end_ll) +
            " exceeds chromosome length " + std::to_string(chrom_len) +
            " for '" + chrom + "' in region '" + region_str + "'.");
    }

    return ngslib::GenomeRegion(chrom,
                                static_cast<uint32_t>(start_ll),
                                static_cast<uint32_t>(end_ll));
}

std::vector<ngslib::GenomeRegion> parse_regions(
    const std::string& regions_str,
    const std::map<std::string, uint32_t>& ref_fai_map,
    const std::vector<std::string>& chroms_filter)
{
    std::vector<std::string> tokens;
    ngslib::split(regions_str, tokens, ",");

    // Optional whitelist
    std::vector<std::string> filter_sorted(chroms_filter);
    std::sort(filter_sorted.begin(), filter_sorted.end());
    const bool has_filter = !filter_sorted.empty();

    std::vector<ngslib::GenomeRegion> out;
    out.reserve(tokens.size());
    for (auto& tok : tokens) {
        // Skip empty tokens (e.g. trailing comma).
        auto not_space = [](unsigned char c) { return !std::isspace(c); };
        std::string t = tok;
        t.erase(t.begin(), std::find_if(t.begin(), t.end(), not_space));
        t.erase(std::find_if(t.rbegin(), t.rend(), not_space).base(), t.end());
        if (t.empty()) continue;

        ngslib::GenomeRegion r = parse_region_string(t, ref_fai_map);
        if (has_filter &&
            !std::binary_search(filter_sorted.begin(), filter_sorted.end(), r.chrom)) {
            continue;
        }
        out.push_back(r);
    }
    return out;
}

// ---------------------------------------------------------------------------
// per sub-region command string
// ---------------------------------------------------------------------------
std::string build_command(const std::string& exe_prog,
                          const std::string& passthrough_str,
                          const std::string& chrom,
                          uint32_t start,
                          uint32_t end,
                          const std::string& outdir)
{
    std::string prefix = chrom + "_" + std::to_string(start) + "_" + std::to_string(end);
    std::filesystem::path out_dir(outdir);
    std::string out_vcf = (out_dir / (prefix + ".vcf.gz")).string();
    std::string out_log = (out_dir / (prefix + ".log")).string();

    std::ostringstream oss;
    oss << "time " << exe_prog;
    if (!passthrough_str.empty()) oss << " " << passthrough_str;
    oss << " -r " << chrom << ":" << start << "-" << end
        << " -o " << out_vcf
        << " > " << out_log
        << " && echo \"** " << prefix << " done **\"";
    return oss.str();
}

// ---------------------------------------------------------------------------
// help text
// ---------------------------------------------------------------------------
static const std::string PIPELINE_USAGE =
    "About: Generate per-region `basevar caller` commands for whole-genome\n"
    "       variant calling.  All non-pipeline options are passed through to\n"
    "       `basevar caller` verbatim, so any new caller option works\n"
    "       automatically without changing this subcommand.\n\n"
    "Usage: basevar pipeline [pipeline options] [caller pass-through options]\n\n"
    "Pipeline options:\n"
    "  -o, --outdir=DIR        Output directory for per-region VCF files\n"
    "                          and logs (required).\n"
    "      --ref_fai=FILE      Reference FASTA index file (.fai) used to\n"
    "                          determine chromosome lengths (required).\n"
    "  -d, --delta=INT         Size of each sub-region in bp [2000000].\n"
    "  -c, --chrom=STR         Only process these comma-delimited\n"
    "                          chromosomes (e.g. chr1,chr2).\n"
    "  -h, --help              Show this help message and exit.\n\n"
    "Common caller pass-through options (forwarded as-is):\n"
    "  -f, --reference FILE         Reference FASTA file.\n"
    "  -L, --align-file-list FILE   BAM/CRAM list, one file per line.\n"
    "  -r, --regions REG[,...]      Restrict to these regions; the pipeline\n"
    "                               will further split them by --delta.\n"
    "  -Q, --min-BQ INT             Minimum base quality.\n"
    "  -q, --mapq INT               Minimum mapping quality.\n"
    "  -B, --batch-count INT        Samples per batch file.\n"
    "  -t, --thread INT             Number of threads.\n"
    "  -G, --pop-group FILE         Per-population allele frequency file.\n"
    "  --filename-has-samplename    BAM filename starts with sample ID.\n"
    "  --gt-mode=STRING             Genotype calling mode.\n"
    "  --ref-bias=FLOAT             Reference bias coefficient (β) for genotype likelihood calculation.\n"
    "  --max-alleles=INT            Maximum number of active alleles allowed at a site.\n"
    "  --smart-rerun                Resume interrupted runs.\n\n"
    "Example:\n"
    "  basevar pipeline \\\n"
    "      -o /path/to/outdir --ref_fai reference.fasta.fai \\\n"
    "      -d 2000000 -c chr20 \\\n"
    "      -f reference.fasta -L bam.list -Q 20 -q 30 -B 500 -t 4 \\\n"
    "      --filename-has-samplename > basevar_wgs.sh";

// ---------------------------------------------------------------------------
// argv helper: pop element at index `i` from args (shifts the tail down).
// ---------------------------------------------------------------------------
static void erase_at(std::vector<std::string>& args, size_t i) {
    args.erase(args.begin() + i);
}

// ---------------------------------------------------------------------------
// subcommand entry point
// ---------------------------------------------------------------------------
int pipeline_runner(int argc, char* argv[], const std::string& basevar_executable) {

    if (argc < 2) {
        std::cout << PIPELINE_USAGE << "\n" << std::endl;
        return 1;
    }

    // ----- 1) Move all argv into a std::vector for easier manipulation -----
    // We KEEP argv[0] in place so getopt_long reports the right program name.
    std::vector<std::string> args;
    args.reserve(argc);
    for (int i = 0; i < argc; ++i) args.emplace_back(argv[i]);

    // ----- 2) Extract pipeline-specific options FIRST, before getopt_long  -----
    //         We walk the vector manually so unknown options are preserved
    //         intact for the pass-through string.  This mirrors the behaviour
    //         of argparse.parse_known_args() in the original Python script.
    PipelineArgs pargs;
    bool show_help = false;

    auto consume_value = [&](size_t& i, const std::string& name,
                             const std::string& cur) -> std::string {
        // Handles "--flag=value", "--flag value" and short "-X value".
        size_t eq = cur.find('=');
        if (eq != std::string::npos) {
            std::string v = cur.substr(eq + 1);
            erase_at(args, i);
            return v;
        }
        if (i + 1 >= args.size()) {
            throw std::invalid_argument(
                "[ERROR] Option '" + name + "' requires a value.");
        }
        std::string v = args[i + 1];
        erase_at(args, i);  // remove flag
        erase_at(args, i);  // remove value (same index after the first erase)
        return v;
    };

    size_t i = 1;  // skip argv[0]
    while (i < args.size()) {
        const std::string& a = args[i];

        if (a == "-h" || a == "--help") {
            show_help = true;
            erase_at(args, i);
            continue;
        }
        if (a == "-o" || a == "--outdir" || a.rfind("--outdir=", 0) == 0) {
            pargs.outdir = consume_value(i, "--outdir", a);
            continue;
        }
        if (a == "--ref_fai" || a.rfind("--ref_fai=", 0) == 0) {
            pargs.ref_fai = consume_value(i, "--ref_fai", a);
            continue;
        }
        if (a == "-d" || a == "--delta" || a.rfind("--delta=", 0) == 0) {
            std::string v = consume_value(i, "--delta", a);
            try {
                long long d = std::stoll(v);
                if (d <= 0) {
                    throw std::invalid_argument(
                        "[ERROR] --delta must be a positive integer, got '" + v + "'.");
                }
                pargs.delta = static_cast<uint32_t>(d);
            } catch (const std::invalid_argument&) {
                throw;
            } catch (const std::exception&) {
                throw std::invalid_argument(
                    "[ERROR] --delta must be an integer, got '" + v + "'.");
            }
            continue;
        }
        if (a == "-c" || a == "--chrom" || a.rfind("--chrom=", 0) == 0) {
            pargs.chrom = consume_value(i, "--chrom", a);
            continue;
        }

        // Unknown option/positional: leave for the pass-through path.
        ++i;
    }

    if (show_help) {
        std::cout << PIPELINE_USAGE << "\n" << std::endl;
        return 0;
    }

    // ----- 3) Validate pipeline-specific args -----
    if (pargs.outdir.empty()) {
        std::cerr << "[ERROR] Missing required option '-o/--outdir'." << std::endl;
        return 1;
    }
    if (pargs.ref_fai.empty()) {
        std::cerr << "[ERROR] Missing required option '--ref_fai'." << std::endl;
        return 1;
    }
    if (!ngslib::is_readable(pargs.ref_fai)) {
        std::cerr << "[ERROR] --ref_fai file is not readable: "
                  << pargs.ref_fai << std::endl;
        return 1;
    }

    // Normalise outdir to an absolute path (matches the Python script).
    pargs.outdir = std::filesystem::absolute(pargs.outdir).string();

    // ----- 4) Build the chromosome-length map from .fai -----
    std::vector<std::string> chroms_filter;
    if (!pargs.chrom.empty()) {
        std::vector<std::string> tmp;
        ngslib::split(pargs.chrom, tmp, ",");
        for (auto& c : tmp) {
            auto not_space = [](unsigned char ch) { return !std::isspace(ch); };
            c.erase(c.begin(), std::find_if(c.begin(), c.end(), not_space));
            c.erase(std::find_if(c.rbegin(), c.rend(), not_space).base(), c.end());
            if (!c.empty()) chroms_filter.push_back(c);
        }
    }

    std::vector<ngslib::GenomeRegion> ref_fai_entries;
    try {
        ref_fai_entries = load_reference_fai(pargs.ref_fai, chroms_filter);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    if (ref_fai_entries.empty()) {
        std::cerr << "[ERROR] No chromosomes loaded from --ref_fai "
                  << "(check --chrom filter or fai contents)." << std::endl;
        return 1;
    }

    std::map<std::string, uint32_t> ref_fai_map;
    for (const auto& e : ref_fai_entries) ref_fai_map[e.chrom] = e.end;

    // ----- 5) Intercept -r/--regions from the pass-through args -----
    std::string regions_str;
    size_t j = 1;  // skip argv[0]
    while (j < args.size()) {
        const std::string& a = args[j];
        if ((a == "-r" || a == "--regions") && j + 1 < args.size()) {
            regions_str = args[j + 1];
            erase_at(args, j);
            erase_at(args, j);
            continue;
        }
        if (a.rfind("--regions=", 0) == 0) {
            regions_str = a.substr(std::string("--regions=").size());
            erase_at(args, j);
            continue;
        }
        ++j;
    }

    // ----- 6) Determine the working region list -----
    std::vector<ngslib::GenomeRegion> region_list;
    if (!regions_str.empty()) {
        try {
            region_list = parse_regions(regions_str, ref_fai_map, chroms_filter);
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] Failed to parse -r/--regions: " << e.what() << std::endl;
            return 1;
        }
        if (region_list.empty()) {
            std::cerr << "[ERROR] -r/--regions produced no valid regions "
                      << "(check --chrom filter or region strings)." << std::endl;
            return 1;
        }
    } else {
        region_list = ref_fai_entries;  // already filtered by --chrom
    }

    // ----- 7) Build the pass-through string (everything except argv[0]) -----
    std::ostringstream pass;
    for (size_t k = 1; k < args.size(); ++k) {
        if (k > 1) pass << " ";
        pass << args[k];
    }
    std::string passthrough_str = pass.str();

    // ----- 8) Resolve the executable program for the generated commands -----
    //         Priority:
    //           1. $basevar environment variable (matches Python script)
    //           2. caller-supplied path (argv[0] from main)
    //           3. literal "basevar"
    std::string exe_prog;
    if (const char* env_p = std::getenv("basevar")) {
        if (env_p[0] != '\0') exe_prog = std::string(env_p) + " caller";
    }
    if (exe_prog.empty()) {
        if (!basevar_executable.empty()) {
            exe_prog = basevar_executable + " caller";
        } else {
            exe_prog = "basevar caller";
        }
    }

    // ----- 9) Emit one command per sub-region -----
    for (const auto& r : region_list) {
        for (uint64_t s = r.start; s <= r.end; s += pargs.delta) {
            uint64_t e = std::min<uint64_t>(s + pargs.delta - 1, r.end);
            std::cout << build_command(exe_prog, passthrough_str, r.chrom,
                                       static_cast<uint32_t>(s),
                                       static_cast<uint32_t>(e),
                                       pargs.outdir)
                      << "\n";
        }
    }

    return 0;
}

}  // namespace pipeline
}  // namespace basevar
