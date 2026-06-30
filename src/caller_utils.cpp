#include "caller_utils.h"
#include <cmath>
#include <cstdio>  // snprintf

#include "io/fasta.h"
#include "io/iobgzf.h"
#include "algorithm.h"

std::string format_double(double value, int precision) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
}

StrandBiasInfo strand_bias(const std::string &major_base, 
                           const std::string &alt_base,
                           const std::vector<std::string> &bases,
                           const std::vector<char> &strands)  // strands 和 bases 是配对的，一一对应 
{
    int maj_fwd = 0, maj_rev = 0;
    int alt_fwd = 0, alt_rev = 0;
    for (size_t i(0); i < bases.size(); ++i) {
        if (bases[i][0] == 'N' || bases[i][0] == 'n') continue;

        if (strands[i] == '+') {
            if (bases[i] == major_base) {
                ++maj_fwd;
            } else if (alt_base == bases[i]) {
                ++alt_fwd;
            }

        } else if (strands[i] == '-') {
            if (bases[i] == major_base) {
                ++maj_rev;
            } else if (alt_base == bases[i]) {
                ++alt_rev;
            }

        } else {
            throw std::runtime_error("[ERROR] in 'strand_bias' Get strange strand symbol: " + bases[i] + " " + std::to_string(strands[i]));
        }
    }

    if (alt_base == major_base) { 
        // 如果 alt_base 刚好是 major_base 那么就按照 50% 的比例构造 alt_fwd 和 alt_rev 的理论值
        alt_fwd = alt_rev = int(std::round(0.5 * (maj_fwd + maj_rev)));
    }

    double fs = -10 * std::log10(fisher_exact_test(maj_fwd, maj_rev, alt_fwd, alt_rev));
    // 应对'全 major allele' 或者是 '全 alt allele'
    if (std::isinf(fs)) {
        fs = 10000;
    } else if (fs == 0) {
        fs = 0.0;
    }

    // Symmetric Odds Ratio (GATK-compatible SOR)
    // Reference: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
    //
    // SOR = ln(symRatio) * (refRatio - altRatio)
    //   symRatio = (Rf*Ar)/(Rr*Af) + (Rr*Af)/(Rf*Ar)
    //   refRatio = Rr/Rf,  altRatio = Ar/Af
    //
    // where Rf=maj_fwd, Rr=maj_rev, Af=alt_fwd, Ar=alt_rev
    // When any cell count is 0, add pseudocount of 1 to all cells (GATK convention)
    double Rf = static_cast<double>(maj_fwd);
    double Rr = static_cast<double>(maj_rev);
    double Af = static_cast<double>(alt_fwd);
    double Ar = static_cast<double>(alt_rev);

    // When any cell count is 0, add pseudocount of 1 to all cells (GATK convention)
    if (maj_fwd == 0 || maj_rev == 0 || alt_fwd == 0 || alt_rev == 0) {
        Rf += 1.0; Rr += 1.0; Af += 1.0; Ar += 1.0;
    }

    double sym_ratio = (Rf * Ar) / (Rr * Af) + (Rr * Af) / (Rf * Ar);
    double ref_ratio = Rr / Rf;
    double alt_ratio = Ar / Af;
    double sor = std::log(sym_ratio) * (ref_ratio - alt_ratio);

    StrandBiasInfo sbi;
    sbi.fwd = (alt_base == major_base) ? maj_fwd : alt_fwd;
    sbi.rev = (alt_base == major_base) ? maj_rev : alt_rev;
    sbi.fs  = fs; 
    sbi.sor = sor;

    return sbi;
}

AlleleInfo collect_and_normalized_allele_info(VariantInfo &variant, std::vector<BaseType::BatchInfo> &smps_bi_v) {
    // Find longest REF
    std::string shared_ref;
    for (const auto& ref : variant.ref_bases) {  // ref_bases are all upper case already
        if (ref.length() > shared_ref.length()) {
            shared_ref = ref;
        }
    }

    AlleleInfo ai; 
    ai.ref = shared_ref;  // set raw ref and it's upper case already

    // Normalize ALTs and update variants
    std::set<std::string> unique_alts;
    std::map<std::string, std::string> raw_bases_map;  // raw bases to normalized bases
    for (size_t j = 0; j < variant.ale_bases.size(); j++) {
        std::string raw_base = variant.ale_bases[j];
        std::string alt(raw_base); // it will be modified based on shared_ref
        std::string ref(variant.ref_bases[j]);

        // rebase if Indels
        if (alt[0] == '-') {
            alt = ref[0];               // replace by the first ref bases for DEL seq
            variant.ale_bases[j] = alt; // Update the ALT sequence
        } else if (alt[0] == '+') {
            alt = ref + alt.substr(1);  // replace the first base('+') by ref bases
            variant.ale_bases[j] = alt; // rewrite deletion seq
        }
        
        if (ref != shared_ref && shared_ref.length() > ref.length()) {
            alt += shared_ref.substr(ref.length()); // Adding suffix seq to fit with shared_ref
            variant.ale_bases[j] = alt;             // Update the ALT sequence
        }

        unique_alts.insert(alt);
        raw_bases_map[raw_base] = alt;  // Store raw bases mapping

        ai.total_dp += variant.depths[j];  // accumulate total depth
        if (alt != shared_ref) {
            // Set allele count and allele frequency
            ai.allele_counts[alt] = 0; // Initialize allele count
            ai.allele_depths[alt] = variant.depths[j];
            ai.allele_freqs[alt]  = variant.freqs[j];
        }
        ai.strand_bias_info[alt] = variant.strand_bias[j];  // Store strand bias information
    }
    unique_alts.erase(shared_ref);  // Remove REF from ALTs

    // Set ALT field: Unique and sorted ALT sequences by length and then by ASCII
    ai.alts = ngslib::get_unique_strings(std::vector<std::string>(unique_alts.begin(), unique_alts.end()));
   
    // Update align bases in smps_bi_v and calculate allele counts
    // Replace align bases with normalized alt bases
    for (auto& smp_bi : smps_bi_v) {
        for (size_t i = 0; i < smp_bi.align_bases.size(); ++i) {
            std::string ab = smp_bi.align_bases[i];
            if (raw_bases_map.find(ab) != raw_bases_map.end()) {
                smp_bi.align_bases[i] = raw_bases_map[ab];  // Update align base to normalized alt
            }
        }
    }

    return ai;
}

// Helper function: Convert PL index to genotype pair
// VCF 4.2 standard PL ordering: j=0..n-1 (outer), k=0..j (inner)
// Genotypes: (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), ...
std::pair<size_t, size_t> pl_index_to_genotype(size_t pl_idx, size_t n_alleles) {
    size_t idx = 0;
    for (size_t j = 0; j < n_alleles; ++j) {
        for (size_t k = 0; k <= j; ++k) {
            if (idx == pl_idx) return {k, j};
            ++idx;
        }
    }
    return {0, 0}; // Default to homozygous reference genotype
}

VCFSampleAnnotation process_sample_variant(const std::string& ref_base,          // It's upper case already
                                           const std::vector<std::string>& alts, // It's upper case already
                                           const BaseType::BatchInfo& smp_bi,    // smp_bi.align_bases should be upper case already
                                           double af,                            // population AF prior (-1 = disable Bayesian)
                                           double ref_bias)                      // reference bias coefficient β (0.5 = no bias)
{
    VCFSampleAnnotation sa;
    
    // Calculate PL values for all possible genotypes in one call
    sa.PL = calculatePL(ref_base, alts, smp_bi.align_bases, smp_bi.align_base_quals, ref_bias);

    size_t n_alleles = alts.size() + 1;  // Total number of alleles including REF

    if (af >= 0) {
        // Bayesian mode: compute genotype posterior using population AF prior
        GenotypePosterior gp = compute_genotype_posterior(sa.PL, af);

        size_t min_idx = gp.best_gt_idx;
        sa.GQ = gp.gq;
        sa.posterior = gp.posteriors;
        sa.dosage = gp.dosage;
        sa.per_allele_dosage = gp.per_allele_dosage;

        // Convert PL index to genotype codes
        auto [a1, a2] = pl_index_to_genotype(min_idx, n_alleles);
        sa.gtcode = {a1, a2};
    } else {
        // Original mode: pure likelihood (no prior)
        size_t min_idx = argmin(sa.PL.begin(), sa.PL.end());

        // Find the second minimum PL value and set GQ
        int second_min_pl = std::numeric_limits<int>::max();
        for (size_t i(0); i < sa.PL.size(); ++i) {
            if (i != min_idx && sa.PL[i] < second_min_pl) {
                second_min_pl = sa.PL[i];
            }
        }
        sa.GQ = static_cast<double>(second_min_pl);

        // Convert PL index to genotype codes
        auto [a1, a2] = pl_index_to_genotype(min_idx, n_alleles);
        sa.gtcode = {a1, a2};
    }

    // Store allele depths in order (REF first, followed by ALTs)
    std::vector<std::string> ref_alts_order = {ref_base};
    ref_alts_order.insert(ref_alts_order.end(), alts.begin(), alts.end());

    // Initialize allele depth counting
    std::map<std::string, int> allele_depths;
    for (const auto& alt : ref_alts_order) {
        allele_depths[alt] = 0;
    }
    
    // Count depth for each allele
    for (const auto& base : smp_bi.align_bases) {
        if (base[0] == 'N' || base[0] == 'n') continue;
        if (allele_depths.find(base) != allele_depths.end()) {
            allele_depths[base]++;
        }
    }
    
    for (const auto& b : ref_alts_order) {
        sa.allele_depths.push_back(allele_depths[b]);
        if (allele_depths[b] > 0) {
            sa.sample_alts.push_back(b);
        }
    }

    return sa;
}

std::string format_sample_string(const VCFSampleAnnotation& sa) {
    if (sa.sample_alts.empty()) {
        return "./.";  // No variants found
    }

    int dp = std::accumulate(sa.allele_depths.begin(), sa.allele_depths.end(), 0);

    // Format GQ: always output as integer (VCF spec: Type=Integer)
    // Legacy mode: GQ = PL-gap (already integer-valued)
    // Bayesian mode: GQ = -10*log10(1-P_best), rounded to nearest integer
    std::string gq_str = std::to_string(static_cast<int>(std::round(sa.GQ)));
    std::string sample_info = ngslib::join(sa.gtcode, "/")        + ":" +  // GT, genotype
                              gq_str                              + ":" +  // GQ, Genotype Quality
                              ngslib::join(sa.PL, ",")            + ":" +  // PL, Phred-scaled likelihoods
                              ngslib::join(sa.allele_depths, ",") + ":" +  // AD, active allele depth, so sum(AD) <= PD
                              std::to_string(dp);                          // DP, high quality mapping total depth

    return sample_info;
}

std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &addition_info, 
                              const std::vector<std::string> &samples, const std::string other_comment)
{
    std::vector<std::string> header = {
        "##fileformat=VCFv4.2",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        
        // FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality: Phred-scaled confidence (posterior-based in Bayesian mode, PL-gap in legacy mode)\">",
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Normalized, Phred-scaled genotype likelihoods rounded to the closest integer\">",
        "##FORMAT=<ID=AD,Number=R,Type=String,Description=\"Allelic depths for the ref and alt alleles in the order listed\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth for specific sample (reads with bad mapped quality or with bad mates are filtered)\">",
        
        // INFO fields
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency: dosage-based (AC/AN) in posterior mode (Recommend by BaseVar); LRT EM frequency in legacy mode\">",
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count: expected allele count derived from genotype posterior dosage in posterior mode (dosage-based); reads-based in legacy mode\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles: 2*N_samples in posterior mode (dosage-based count); reads-based in legacy mode\">",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate total read depth (high-quality); some reads may have been filtered\">",
        "##INFO=<ID=DP4,Number=A,Type=Integer,Description=\"A list of number of high-quality ref-forward, ref-reverse, alt1-forward, alt1-reverse, alt2-forward, alt2-reverse, ... ,bases\">",
        "##INFO=<ID=FS,Number=A,Type=Float,Description=\"Phred-scaled P-value using Fisher's exact test to detect strand bias\">",
        "##INFO=<ID=SOR,Number=A,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">",
        "##INFO=<ID=AC_obs,Number=A,Type=Integer,Description=\"Observed allele count from hard genotype calls (GT-based), posterior mode only\">",
        "##INFO=<ID=AN_obs,Number=1,Type=Integer,Description=\"Total number of observed alleles from called (non-missing) genotype calls, posterior mode only\">",
        "##INFO=<ID=AF_obs,Number=A,Type=Float,Description=\"Observed allele frequency = AC_obs / AN_obs, posterior mode only\">",
        "##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Hardy-Weinberg equilibrium chi-square p-value (dosage-based, suitable for low-coverage data)\">",
        "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">",
        "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">",
        "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">"
    };  // initial by common information of header
    if (!addition_info.empty()) header.insert(header.end(), addition_info.begin(), addition_info.end());

    ngslib::Fasta fa = ref_file_path;
    std::vector<std::string> contigs;
    for (size_t i(0); i < fa.nseq(); ++i) {
        std::string seqname = fa.iseq_name(i);
        uint32_t seqlen = fa.seq_length(seqname);
        contigs.push_back("##contig=<ID=" + seqname + ",length=" + std::to_string(seqlen) + 
                          ",assembly=" + ref_file_path + ">");
    }

    header.insert(header.end(), contigs.begin(), contigs.end());
    header.push_back("##reference=file://" + ngslib::abspath(ref_file_path));
    if (!other_comment.empty()) header.push_back(other_comment);
    header.push_back("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + ngslib::join(samples, "\t"));

    return ngslib::join(header, "\n");
}

std::string cvg_header_define(const std::vector<std::string> &group_info, const std::vector<std::string> &BASES) {

    std::string h = "#CHROM\tPOS\tREF\tDepth\t" + ngslib::join(BASES, "\t") + "\t" +
                    "Indels\tFS\tSOR\tStrand_Coverage(REF_FWD,REF_REV,ALT_FWD,ALT_REV)";
    // group cvg 的计算有 bug， 暂时不输出
    // if (!group_info.empty())
    //     h += "\t" + ngslib::join(group_info, "\t");

    std::vector<std::string> header = {
        "##fileformat=CVGv1.0",
        "##Group information is the depth of A:C:G:T:Indel", 
        h
    };

    return ngslib::join(header, "\n");
}

void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile, 
                        std::string header, bool is_remove_tempfile)
{
    if (infiles.empty()) return;

    bool is_compress = (ngslib::suffix_name(outfile) == ".gz") ? true : false;
    ngslib::BGZFile OUT(outfile, is_compress ? "wb" : "uw"); 
    OUT << header << "\n";

    /* Merge all files here */
    for (auto fn: infiles) {
        ngslib::BGZFile f(fn, "r");
        std::string line;

        while (f.readline(line)) {
            if (line[0] == '#') continue;
            OUT << line << "\n";
        }
        OUT.flush(); // 确保数据被写入

        if (is_remove_tempfile) ngslib::safe_remove(fn);
    }

    OUT.close();
    return;
}
