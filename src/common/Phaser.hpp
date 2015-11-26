#ifndef __PHASER_HPP
#define __PHASER_HPP

#include "Mapping.hpp"
#include "Het_Variation.hpp"

class Phaser
{
public:
    typedef Het_Variation::DNA_Sequence DNA_Sequence;
    static int phase(const Mapping & m, const Het_Variation & v);
private:
    static int phase_snp(const Mapping & m, const Het_Variation & v);
    static int phase_indel(const Mapping & m, const Het_Variation & v);
}; // class Phaser

#endif
