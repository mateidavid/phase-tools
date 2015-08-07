#ifndef __PHASED_SET_CONNECTION_HPP
#define __PHASED_SET_CONNECTION_HPP

#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <boost/intrusive/set.hpp>

#include "Phased_Set.hpp"

class Phased_Set_Connection
    : public boost::intrusive::set_base_hook<>
{
public:
    Phased_Set_Connection()
        : _ps_ptr{nullptr, nullptr}, _count{0, 0, 0, 0}, _discordance(1.0)
    {}
    Phased_Set_Connection(const Phased_Set * ps0_ptr, const Phased_Set * ps1_ptr)
        : _ps_ptr{ps0_ptr, ps1_ptr}, _count{0, 0, 0, 0}, _discordance(1.0)
    {}

    const Phased_Set * ps_ptr(int i) const { return _ps_ptr[i]; }
    void set_counts(size_t (&cnt)[2][2])
    {
        _count[0][0] = cnt[0][0];
        _count[0][1] = cnt[0][1];
        _count[1][0] = cnt[1][0];
        _count[1][1] = cnt[1][1];
        recompute_discordance();
    }
    const size_t & count(int i, int j) const { return _count[i][j]; }
    void increment_count(int i, int j) { ++_count[i][j]; recompute_discordance(); }
    const double & discordance() const { return _discordance; }
    bool rel_phase() const { return _count[0][1] + _count[1][0] > _count[0][0] + _count[1][1]; }

    void swap_ps()
    {
        std::swap(_ps_ptr[0], _ps_ptr[1]);
        std::swap(_count[0][1], _count[1][0]);
    }
    void swap_phase(int i)
    {
        if (i == 0)
        {
            std::swap(_count[0][0], _count[1][0]);
            std::swap(_count[0][1], _count[1][1]);
        }
        else
        {
            std::swap(_count[0][0], _count[0][1]);
            std::swap(_count[1][0], _count[1][1]);
        }
    }

    void replace_ps_ptr(const Phased_Set * old_ps_ptr, const Phased_Set * new_ps_ptr, bool switch_phase)
    {
        assert(old_ps_ptr == _ps_ptr[0] or old_ps_ptr == _ps_ptr[1]);
        int i = (old_ps_ptr == _ps_ptr[1]);
        _ps_ptr[i] = new_ps_ptr;
        if (switch_phase)
        {
            swap_phase(i);
        }
    }

    void merge(Phased_Set_Connection & other, Phased_Set * common_ps_ptr, bool switch_phase)
    {
        if (common_ps_ptr == _ps_ptr[1]) swap_ps();
        if (common_ps_ptr == other._ps_ptr[1]) other.swap_ps();
        assert(common_ps_ptr == _ps_ptr[0]);
        assert(common_ps_ptr == other._ps_ptr[0]);
        if (switch_phase) other.swap_phase(1);
        _count[0][0] += other._count[0][0];
        _count[0][1] += other._count[0][1];
        _count[1][0] += other._count[1][0];
        _count[1][1] += other._count[1][1];
        recompute_discordance();
        other._ps_ptr[0] = nullptr;
        other._ps_ptr[1] = nullptr;
        other._count[0][0] = 0;
        other._count[0][1] = 0;
        other._count[1][0] = 0;
        other._count[1][1] = 0;
    }

    friend bool operator < (const Phased_Set_Connection& lhs, const Phased_Set_Connection& rhs)
    {
        return lhs.discordance() < rhs.discordance();
    }

    friend std::ostream& operator << (std::ostream& os, const Phased_Set_Connection& c)
    {
        /*
        os << static_cast< const void * >(c.ps_ptr(0)) << '\t'
           << static_cast< const void * >(c.ps_ptr(1)) << '\t'
           << c.count(0, 0) << '\t' << c.count(0, 1) << '\t' << c.count(1, 0) << '\t' << c.count(1, 1) << '\t'
           << c.discordance() << std::endl;
        */
        os << *c.ps_ptr(0) << ',' << *c.ps_ptr(1) << ':'
           << c.count(0, 0) << ',' << c.count(0, 1) << ',' << c.count(1, 0) << ',' << c.count(1, 1) << ':'
           << c.discordance();
        return os;
    }

    /*
     * Set maximum discordance to be reported for connections with zero contrary evidence.
     * This is needed in order to allow such connections to pass the discordance test,
     * but with the maximum allowable discordance, ensuring they are processed last.
     */
    static double & max_discordance()
    {
        static double _max_discordance = .5;
        return _max_discordance;
    }

private:
    const Phased_Set * _ps_ptr[2];
    size_t _count[2][2];
    double _discordance;

    void recompute_discordance()
    {
        size_t a = _count[0][0] + _count[1][1];
        size_t b = _count[0][1] + _count[1][0];

        if (a == 0 and b == 0)
        {
            _discordance = 1.0;
        }
        else if (a == 0 or b == 0)
        {
            _discordance = std::min(max_discordance(), 1.0 / (1.0 + a + b));
        }
        else
        {
            _discordance = (1.0 + std::min(a, b)) / (1.0 + a + b);
        }
    }
}; // class Phased_Set_Connection

typedef std::deque< Phased_Set_Connection > Connection_Store;
typedef boost::intrusive::multiset< Phased_Set_Connection > Connection_Set;
typedef std::map< const Phased_Set *, std::map< const Phased_Set *, Phased_Set_Connection * > > Connection_Map;

#endif
