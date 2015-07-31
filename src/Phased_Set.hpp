#ifndef __PHASED_SET_HPP
#define __PHASED_SET_HPP

#include <iostream>
#include <string>
#include <vector>
#include <htslib/vcf.h>

#include "Het_Variation.hpp"

namespace detail
{

struct Phased_Set_Comparator
{
    bool operator () (const std::pair< const Het_Variation *, bool > & lhs,
                      const std::pair< const Het_Variation *, bool > & rhs)
    {
        if (Het_Variation_Ptr_Comp()(lhs.first, rhs.first)) return true;
        else if (Het_Variation_Ptr_Comp()(rhs.first, lhs.first)) return false;
        else return lhs.second < rhs.second;
    }
}; // struct Phased_Set_Comparator

} // namespace detail


class Phased_Set
{
public:
    typedef std::set< std::pair< const Het_Variation *, bool >, detail::Phased_Set_Comparator > Het_Set;

    Phased_Set() = default;
    Phased_Set(const Het_Variation & het)
    {
        _het_set.insert(std::make_pair(&het, false));
    }

    const Het_Set & het_set() const { return _het_set; }
    Het_Set & het_set() { return _het_set; }

    friend std::ostream& operator << (std::ostream& os, const Phased_Set& s)
    {
        os << "(";
        for (const auto & p : s.het_set())
        {
            const auto & v = *p.first;
            bool orientation = p.second;
            //os << v.chr_name() << '\t' << v.rf_start() << '\t' << v.rf_end() << '\t'
            //   << v.gt(orientation) << '|' << v.gt(1 - orientation) << std::endl;
            os << (&p != &*s.het_set().begin()? "," : "")
               << v.chr_name() << ":" << v.rf_start() + 1 << ":" << orientation;
        }
        os << ")";
        return os;
    }

private:
    Het_Set _het_set;

}; // class Phased_Set

#endif

