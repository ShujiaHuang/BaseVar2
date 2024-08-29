/**
 * @file basetype.cpp
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 */

#include "basetype.h"

//////////////////////////////////////////////////////////////
//// The codes for the member function of BaseType class /////
//////////////////////////////////////////////////////////////
BaseType::BaseType(const BatchInfo *smp_bi, double af) : _min_af(af), _only_call_snp(true) {
    ind_allele_likelihood.reserve(smp_bi->n);

    for (size_t i(0); i < smp_bi->n; ++i) {

        char fb = smp_bi->align_bases[i][0];  // a string, get first base
        if (fb != 'N' && (!_only_call_snp || (fb != '+' && fb != '-'))) { 
            // ignore all the 'N' bases and indels if only call snp

            if (smp_bi->align_bases[i].size() != 1)
                throw std::runtime_error("[ERROR] Why dose the size of aligned base is not 1? Check: " + 
                                            smp_bi->align_bases[i]); 

            std::vector<double> allele_lh = {0, 0, 0, 0};  // match allele likelihood for [A, C, G, T]
            for(size_t j(0); j < BASES.size(); ++j) {
                // convert the quality phred scale to be the base confident probabilty value
                if (fb == BASES[j]) {
                    allele_lh[j] = 1.0 - exp((smp_bi->align_base_quals[i] - 33) * MLN10TO10);
                } else {
                    allele_lh[j] = exp((smp_bi->align_base_quals[i] - 33) * MLN10TO10) / 3;
                }
std::cout << "Qual: " << j << " - " << allele_lh[j] << "\n";
            }
            ind_allele_likelihood.push_back(allele_lh);
        }
    }
}