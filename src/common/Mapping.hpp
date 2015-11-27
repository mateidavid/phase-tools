#ifndef __MAPPING_HPP
#define __MAPPING_HPP

#include <string>
#include <sstream>
#include <htslib/sam.h>
#include "RC_Sequence.hpp"
#include "Cigar.hpp"

class Mapping
{
public:
    Mapping() : _rec_p(nullptr) {}
    Mapping(const bam_hdr_t * hdr_p, bam1_t * rec_p)
        : _rec_p(rec_p)
    {
        _query_name = bam_get_qname(rec_p);
        _flag = rec_p->core.flag;
        for (int i = 0; i < rec_p->core.l_qseq; ++i)
        {
            _seq += seq_nt16_str[bam_seqi(bam_get_seq(rec_p), i)];
        }
        if (is_mapped())
        {
            _chr_name = hdr_p->target_name[rec_p->core.tid];
            _rf_start = rec_p->core.pos;
            _cigar = Cigar(bam_get_cigar(rec_p), rec_p->core.n_cigar);
            _rf_len = _cigar.rf_len();
        }
        if (is_paired() and mp_is_mapped())
        {
            _mp_chr_name = hdr_p->target_name[rec_p->core.mtid];
            _mp_rf_start = rec_p->core.mpos;
        }
    }

    bam1_t * rec_p() const { return _rec_p; }
    bool is_paired() const { return _flag & BAM_FPAIRED; }
    bool is_primary() const { return not (_flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)); }
    bool is_mapped() const { return not (_flag & BAM_FUNMAP); }
    bool mp_is_mapped() const { return not (_flag & BAM_FMUNMAP); }
    bool treat_as_paired() const
    {
        return is_paired()
            and (not is_mapped() or not mp_is_mapped() or chr_name() == mp_chr_name());
    }

    const std::string & query_name () const { return _query_name; }
    const std::string & chr_name () const { return _chr_name; }
    const std::string & mp_chr_name () const { return _mp_chr_name; }
    const rc_sequence::Sequence< std::string > & seq() const { return _seq; }
    const Cigar & cigar() const { return _cigar; }
    int rf_start() const { return _rf_start; }
    int rf_len() const { return _rf_len; }
    int rf_end() const { return rf_start() + rf_len(); }
    int mp_rf_start() const { return _mp_rf_start; }

    friend std::ostream & operator << (std::ostream & os, const Mapping & m)
    {
        os
            << m._query_name << "\t"
            << (m.is_paired()? "paired" : "unpaired") << "\t"
            << (m.is_primary()? "primary" : "non-primary") << "\t"
            << (m.is_mapped()? m._chr_name : std::string()) << "\t"
            << (m.is_mapped()? m._rf_start : -1) << "\t"
            << (m.is_mapped()? m._rf_len : -1) << "\t"
            << (m.is_paired() and m.mp_is_mapped()? m._mp_chr_name : std::string()) << "\t"
            << (m.is_paired() and m.mp_is_mapped()? m._mp_rf_start : -1) << "\t"
            << (m.is_mapped()? m._cigar.to_string() : std::string()) << "\t"
            << m._seq;
        return os;
    }

private:
    bam1_t * _rec_p;
    std::string _query_name;
    std::string _chr_name;
    std::string _mp_chr_name;
    rc_sequence::Sequence< std::string > _seq;
    Cigar _cigar;
    int _flag;
    int _rf_start;
    int _rf_len;
    int _mp_rf_start;
}; // class Mapping

#endif
