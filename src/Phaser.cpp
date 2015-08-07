#include "Phaser.hpp"
#include "overlapper.h"

using namespace std;

int Phaser::phase_snp(const Mapping & m, const Het_Variation & v)
{
    assert(v.is_snp());
    assert(m.rf_start() <= v.rf_start() and v.rf_end() <= m.rf_start() + m.rf_len());
    // compute the mapped range, on query, of the single reference position
    auto rg = m.cigar().mapped_range(make_pair(v.rf_start() - m.rf_start(), v.rf_start() - m.rf_start() + 1), true, true);
    if (rg.second != rg.first + 1)
    {
        // the changed reference base is mapped to something other than a match; we give up
        return -1;
    }
    else
    {
        // the changed reference base is mapped to seq[rg.first]
        assert(0 <= rg.first and rg.second <= static_cast< int >(m.seq().size()));
        if (m.seq()[rg.first] == v.allele_seq(v.gt(0))[0])
            return 0;
        else if (m.seq()[rg.first] == v.allele_seq(v.gt(1))[0])
            return 1;
        else
            return -1;
    }
}

int Phaser::phase_indel(const Mapping & m, const Het_Variation & v)
{
    // We first estimate the mapping of the query to the flanking regions of the
    // indel.
    // delta = estimate of the max indel size in flank mapping
    // We use delta to compute a band limit for the DP.
    pair< int, int > rf_rg[2];
    pair< int, int > qr_rg[2];
    int delta[2];
    int flank_len = max(v.flank_seq(0).size(), v.flank_seq(1).size());
    assert(v.rf_start() <= m.rf_end());
    assert(m.rf_start() <= v.rf_end());
    rf_rg[0] = make_pair(max(v.rf_start() - flank_len - m.rf_start(), 0),
                         max(v.rf_start() - m.rf_start(), 0));
    rf_rg[1] = make_pair(min(v.rf_end() - m.rf_start(), m.rf_len()),
                         min(v.rf_end() + flank_len - m.rf_start(), m.rf_len()));
    qr_rg[0] = m.cigar().mapped_range(rf_rg[0], true, false);
    qr_rg[1] = m.cigar().mapped_range(rf_rg[1], true, false);
    delta[0] = abs((qr_rg[0].second - qr_rg[0].first) - (rf_rg[0].second - rf_rg[0].first));
    delta[1] = abs((qr_rg[1].second - qr_rg[1].first) - (rf_rg[1].second - rf_rg[1].first));

    DNA_Sequence rf_allele[2];
    rf_allele[0] = v.flank_seq(0) + v.allele_seq(v.gt(0)) + v.flank_seq(1);
    rf_allele[1] = v.flank_seq(0) + v.allele_seq(v.gt(1)) + v.flank_seq(1);
    int max_allele_seq_size = max(v.allele_seq(v.gt(0)).size(), v.allele_seq(v.gt(1)).size());
    DNA_Sequence qr;
    bool anchor_head = true;
    int band_width;
    if (rf_rg[0].second - rf_rg[0].first == flank_len
        and (rf_rg[1].second - rf_rg[1].first < flank_len
             or delta[0] < delta[1]))
    {
        // use 5p flank
        int qr_end = qr_rg[0].second;
        qr_end = min(qr_end + max_allele_seq_size + flank_len, static_cast< int >(m.seq().size()));
        qr = m.seq().substr(qr_rg[0].first, qr_end - qr_rg[0].first);
        anchor_head = true;
        band_width = delta[0];
    }
    else if (rf_rg[1].second - rf_rg[1].first == flank_len)
    {
        // use 3p flank
        int qr_start = qr_rg[1].first;
        qr_start = max(qr_start - max_allele_seq_size - flank_len, 0);
        qr = m.seq().substr(qr_start, qr_rg[1].second - qr_start);
        anchor_head = false;
        band_width = delta[1];
    }
    else
    {
        // neither flank is fully captured by the mapping, we give up
        //++global::indel_phasing_neither_flank;
        return -1;
    }
    SequenceOverlap res;
    int score[2];
    res = Overlapper::extendMatch(rf_allele[0], qr,
                                  (anchor_head? 0 : rf_allele[0].size() - 1),
                                  (anchor_head? 0 : qr.size() - 1),
                                  band_width);
    score[0] = res.score;
    res = Overlapper::extendMatch(rf_allele[1], qr,
                                  (anchor_head? 0 : rf_allele[1].size() - 1),
                                  (anchor_head? 0 : qr.size() - 1),
                                  band_width);
    score[1] = res.score;
    return (score[0] == score[1]? -1 : (score[0] > score[1]? 0 : 1));
}

int Phaser::phase(const Mapping & m, const Het_Variation & v)
{
    // for now, only deal with SNPs
    int res = v.is_snp()
        ? phase_snp(m, v)
        : phase_indel(m, v);
    // update counts
    ++v.frag_total;
    if (res >= 0) ++v.frag_supp_allele[res];
    return res;
}
