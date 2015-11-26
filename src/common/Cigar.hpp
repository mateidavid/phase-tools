#ifndef __CIGAR_HPP
#define __CIGAR_HPP

#include <tuple>
#include <vector>
#include <string>
#include <sstream>
#include <htslib/sam.h>

class Cigar
{
public:
    struct Pos
    {
        int rf_offset;
        int qr_offset;
        size_t op_idx;
        int op_offset;
    }; // struct Pos

    Cigar() = default;
    Cigar(uint32_t * cigar_ptr, int cigar_len)
        : _v(cigar_ptr, cigar_ptr + cigar_len),
          _rf_len(bam_cigar2rlen(cigar_len, cigar_ptr)),
          _qr_len(bam_cigar2qlen(cigar_len, cigar_ptr))
    {}

    int rf_len() const { return _rf_len; }
    int qr_len() const { return _qr_len; }

    Pos start_pos() const { return Pos{0, 0, 0, 0}; }
    Pos end_pos() const { return Pos{rf_len(), qr_len(), _v.size(), 0}; }

    int op(size_t i) const { return bam_cigar_op(_v.at(i)); }
    int op_len(size_t i) const { return bam_cigar_oplen(_v.at(i)); }
    char op_chr(size_t i) const { return bam_cigar_opchr(_v.at(i)); }
    bool op_act_on_rf(size_t i) const { return bam_cigar_type(_v.at(i)) & 0x2; }
    bool op_act_on_qr(size_t i) const { return bam_cigar_type(_v.at(i)) & 0x1; }
    bool op_act_on(size_t i, bool on_ref) const { return on_ref? op_act_on_rf(i) : op_act_on_qr(i); }

    std::string to_string() const
    {
        std::ostringstream oss;
        for (size_t i = 0; i < _v.size(); ++i)
        {
            oss << op_len(i) << op_chr(i);
        }
        return oss.str();
    }

    /**
     * Starting at given cigar position @pos, advance @len on reference if @on_ref, else on query.
     * NOTE: With len>0: position returned is leftmost possible (i.e., not including length-0 ops);
     * with len==0: position is rightmost possible (i.e., skipping past length-0 ops).
     */
    Pos advance_pos(Pos pos, int len, bool on_ref) const
    {
        assert(len >= 0);
        assert(pos.op_idx <= _v.size());
        assert(pos.op_idx == _v.size() or pos.op_offset < op_len(pos.op_idx));
        assert(pos.op_idx < _v.size() or pos.op_offset == 0);
        bool skip_zero_len_ops = len == 0;
        while (pos.op_idx < _v.size())
        {
            int op_leftover = op_len(pos.op_idx) - pos.op_offset;
            assert(op_leftover > 0);
            //assert(op_act_on(pos.op_idx, true) or op_act_on(pos.op_idx, false));
            if (not op_act_on(pos.op_idx, on_ref) or op_leftover <= len)
            {
                // fully advance past current op
                if (op_act_on(pos.op_idx, on_ref)) len -= op_leftover;
                pos.rf_offset += op_act_on_rf(pos.op_idx)? op_leftover : 0;
                pos.qr_offset += op_act_on_qr(pos.op_idx)? op_leftover : 0;
                pos.op_offset = 0;
                ++pos.op_idx;
                if (len == 0 and not skip_zero_len_ops) break;
            }
            else
            {
                // partly advance in current op
                pos.rf_offset += op_act_on_rf(pos.op_idx)? len : 0;
                pos.qr_offset += op_act_on_qr(pos.op_idx)? len : 0;
                pos.op_offset += len;
                break;
            }
        }
        return pos;
    }

    /**
     * Given a range @rg of offsets on reference if @on_ref, else on query,
     * compute the mapped range of offsets under the cigar mapping.
     */
    std::pair< int, int > mapped_range(const std::pair< int, int >& rg, bool on_ref,
                                       bool include_end_indels) const
    {
        assert(rg.first <= rg.second);
        auto rg_start_pos = advance_pos(start_pos(), rg.first, on_ref);
        if (not include_end_indels) rg_start_pos = advance_pos(rg_start_pos, 0, on_ref);
        auto rg_end_pos = advance_pos(rg_start_pos, rg.second - rg.first, on_ref);
        if (include_end_indels) rg_end_pos = advance_pos(rg_end_pos, 0, on_ref);
        return on_ref
            ? std::make_pair(rg_start_pos.qr_offset, rg_end_pos.qr_offset)
            : std::make_pair(rg_start_pos.rf_offset, rg_end_pos.rf_offset);
    }

private:
    std::vector< uint32_t > _v;
    int _rf_len;
    int _qr_len;
}; // class Cigar

#endif
