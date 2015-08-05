#ifndef __HET_VARIATION_HPP
#define __HET_VARIATION_HPP

#include <iostream>
#include <string>
#include <vector>
#include <htslib/vcf.h>

#include "RC_Sequence.hpp"

class Phased_Set;

class Het_Variation
{
public:
    typedef rc_sequence::Sequence< std::string > DNA_Sequence;

    Het_Variation() = default;
    Het_Variation(int rf_start) : _rf_start(rf_start) {}
    Het_Variation(const bcf_hdr_t * hdr_p, bcf1_t * rec_p, int *& dat, int & dat_size, bool assign_random_phase)
    {
        _chr_name = bcf_hdr_id2name(hdr_p, rec_p->rid);
        _rf_start = rec_p->pos;
        _rf_len = rec_p->rlen;
        for (int i = 0; i < rec_p->n_allele; ++i)
        {
            _allele_seq_v.emplace_back(rec_p->d.allele[i]);
        }
        _n_gt = bcf_get_genotypes(hdr_p, rec_p, &dat, &dat_size);
        _is_het = (_n_gt == 2
                   and not bcf_gt_is_missing(dat[0])
                   and not bcf_gt_is_missing(dat[1])
                   and bcf_gt_allele(dat[0]) != bcf_gt_allele(dat[1]));
        if (_is_het)
        {
            for (int i = 0; i < 2; ++i)
            {
                _gt[i] = bcf_gt_allele(dat[i]);
            }
            _is_phased = bcf_gt_is_phased(dat[1]);
            _is_snp = (_allele_seq_v[0].size() == 1
                       and _allele_seq_v[_gt[0]].size() == 1
                       and _allele_seq_v[_gt[1]].size() == 1);
            if (not _is_phased and assign_random_phase)
            {
                // assign a random phase
                int flip_phase = lrand48() % 2;
                if (flip_phase) std::swap(_gt[0], _gt[1]);
            }
        }
        else
        {
            // keep raw genotypes to be able to output them
            _raw_gt.assign(dat, dat + dat_size);
        }
        frag_total = 0;
        frag_supp_allele[0] = 0;
        frag_supp_allele[1] = 0;
        frag_conflicting = 0;
        ps_phase = 0;
        ps_start_1 = rf_start() + 1;
    }

    const std::string & chr_name() const { return _chr_name; }
    const DNA_Sequence & flank_seq(int i) const { return _flank_seq[i]; }
    const DNA_Sequence & allele_seq(int i) const { return _allele_seq_v.at(i); }
    int n_gt() const { return _n_gt; }
    const int & rf_start() const { return _rf_start; }
    int & rf_start() { return _rf_start; }
    int rf_len() const { return _rf_len; }
    int rf_end() const { return rf_start() + rf_len(); }
    int gt(int i) const { return _gt[i]; }
    const std::vector< int > & raw_gt() const { return _raw_gt; }
    size_t n_alleles() const { return _allele_seq_v.size(); }
    bool is_het() const { return _is_het; }
    bool is_phased() const { return _is_phased; }
    bool is_snp() const { return _is_snp; }
    void reset_phase() { _is_phased = false; if (_gt[0] > _gt[1]) std::swap(_gt[0], _gt[1]); }

    void load_flanks(const faidx_t * faidx_p, int len)
    {
        char * seq;
        int seq_size;
        seq = faidx_fetch_seq(faidx_p, chr_name().c_str(),
                              rf_start() - len, rf_start() - 1,
                              &seq_size);
        _flank_seq[0] = seq;
        free(seq);
        seq = faidx_fetch_seq(faidx_p, chr_name().c_str(),
                              rf_end(), rf_end() + len - 1,
                              &seq_size);
        _flank_seq[1] = seq;
        free(seq);
    }

    friend bool operator < (const Het_Variation& lhs, const Het_Variation& rhs)
    {
        return lhs.rf_start() < rhs.rf_start();
    }

    friend std::ostream& operator << (std::ostream& os, const Het_Variation& v)
    {
        os << v.chr_name() << "\t" << v.rf_start() + 1 << "\t";
        for (size_t i = 0; i < v.n_alleles(); ++i)
        {
            os << v.allele_seq(i) << ((i > 0 and i < v.n_alleles() - 1)? "," : "\t");
        }
        os << v.gt(0) << (v.is_phased()? "|" : "/") << v.gt(1) << "\t"
           << v.flank_seq(0) << "\t" << v.flank_seq(1) << "\t"
           << (v.is_snp()? "snp" : "indel") << "\t"
           << v.frag_total << "\t"
           << v.frag_supp_allele[0] << "\t" << v.frag_supp_allele[1] << "\t"
           << (v.frag_total - (v.frag_supp_allele[0] + v.frag_supp_allele[1])) << "\t"
           << v.frag_conflicting;
        return os;
    }

    mutable size_t frag_total;
    mutable size_t frag_supp_allele[2];
    mutable size_t frag_conflicting;
    mutable Phased_Set * phased_set_ptr;
    mutable int ps_phase;
    mutable int ps_start_1;

private:
    std::string _chr_name;
    DNA_Sequence _flank_seq[2];
    std::vector< DNA_Sequence > _allele_seq_v;
    std::vector< int > _raw_gt;
    int _n_gt;
    int _rf_start;
    int _rf_len;
    int _gt[2];
    bool _is_het;
    bool _is_phased;
    bool _is_snp; // true iff all of (potentially 3) alleles: ref, gt[0], gt[1] are length 1

}; // class Het_Variation

struct Het_Variation_Ptr_Comp
{
    bool operator () (const Het_Variation * lhs, const Het_Variation * rhs) const
    {
        if (!lhs) return true;
        else if (!rhs) return false;
        else return *lhs < *rhs;
    }
}; // struct Het_Variation_Ptr_Comp

#endif

