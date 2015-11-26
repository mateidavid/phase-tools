//-----------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __RC_SEQUENCE_HPP
#define __RC_SEQUENCE_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <mutex>
#include <stdexcept>
#include <cctype>
#include <boost/iterator/transform_iterator.hpp>
#include "shortcuts.hpp"


namespace rc_sequence
{

/** DNA complementer class.
 *
 * The class provides static methods base() and sequence() which return
 * the complement of a base and of a sequence, respectively.
 *
 * The class initializes the complement table internally, in a thread-safe
 * manner using std::mutex from C++11.
 */
class DNA_Complementer;

/** Sequence proxy class.
 *
 * Its primary role is to accumulate rev(), comp(), and substr() operations
 * applied on an existing Sequence object. The Sequence object must not be
 * destroyed while one of its Proxy objects is used!
 *
 * One can use a proxy object to:
 * - accumulate rev(), comp(), and substr() operations
 * - retrieve a given position using at() and operator[]()
 * - compare with other sequence proxy objects
 * - print object (ostream& operator <<)
 * - construct a new Sequence object
 */
template < typename String_Type, typename Complementer >
class Sequence_Proxy;

/** Sequence class.
 *
 * Implemented as a non-orthodox derivation from a regular string class.
 * The derived class adds no new data members, but provides new functions
 * rev(), comp(), revcomp(), and overloads substr(). All of these produce
 * Sequence_Proxy objects that accumulate the corresponding restrictions.
 */
template < typename String_Type, typename Complementer >
class Sequence;

class DNA_Complementer
{
public:
    /** Function that complements one base. */
    static char base(char c)
    {
        int res = table().at(static_cast< unsigned >(c));
        assert(res >= 0);
        return static_cast< char >(res);
    }
    /** Function that complements a sequence. */
    template < typename Sequence_Type >
    static Sequence_Type sequence(const Sequence_Type& s)
    {
        return Sequence_Type(iterator(s.cbegin()), iterator(s.cend()));
    }
private:
    /** Iterator modifier that complements the output of the iterator. */
    template < typename Iterator >
    static boost::transform_iterator< decltype(&base), Iterator > iterator(const Iterator& it)
    {
        return boost::make_transform_iterator(it, &base);
    }
    /** Complement table holder (thread-safe initialization). */
    static const std::vector< int >& table()
    {
        static std::vector< int > _complement_table(128, -1);
        static const std::vector< std::vector< char > > _complement_table_pairs =
            { { 'A', 'T' },
              { 'C', 'G' },
              { 'N', 'N' },
              { 'X', 'X' },
              { '.', '.' } };
        static bool _inited = false;
        if (not _inited)
        {
            static std::mutex _mutex;
            {
                std::lock_guard< std::mutex > _lock(_mutex);
                if (not _inited)
                {
                    for (const auto& p : _complement_table_pairs)
                    {
                        _complement_table[static_cast< int >(tolower(p[0]))] = tolower(p[1]);
                        _complement_table[static_cast< int >(toupper(p[0]))] = toupper(p[1]);
                        _complement_table[static_cast< int >(tolower(p[1]))] = tolower(p[0]);
                        _complement_table[static_cast< int >(toupper(p[1]))] = toupper(p[0]);
                    }
                }
            }
            _inited = true;
        }
        return _complement_table;
    }
}; // class DNA_Complementer

namespace detail
{

/** Iterator for proxy objects */
template < typename Complementer >
class Sequence_Proxy_Iterator
    : public boost::iterator_adaptor< Sequence_Proxy_Iterator< Complementer >,
                                      const char*,
                                      char,
                                      boost::random_access_traversal_tag,
                                      char >
{
private:
    typedef boost::iterator_adaptor< Sequence_Proxy_Iterator< Complementer >,
                                     const char*,
                                     char,
                                     boost::random_access_traversal_tag,
                                     char > Adaptor_Base;
public:
    Sequence_Proxy_Iterator()
        : Sequence_Proxy_Iterator::iterator_adaptor_(0), _reversed(false), _complemented(false) {}
    explicit Sequence_Proxy_Iterator(const char* other, bool reversed = false, bool complemented = false)
        : Sequence_Proxy_Iterator::iterator_adaptor_(other), _reversed(reversed), _complemented(complemented) {}

    using typename Adaptor_Base::value_type;
    using typename Adaptor_Base::reference;
    using typename Adaptor_Base::difference_type;

private:
    friend class boost::iterator_core_access;

    char dereference() const
    {
        char base = not _reversed? *(this->base()) : *(this->base() - 1);
        return not _complemented? base : Complementer::base(base);
    }
    void increment() { increment_dir(true); }
    void decrement() { increment_dir(false); }
    void increment_dir(bool dir)
    {
        if (dir != _reversed)
        {
            ++this->base_reference();
        }
        else
        {
            --this->base_reference();
        }
    }
    void advance(difference_type n)
    {
        this->base_reference() += not _reversed? n : -n;
    }
    difference_type distance_to(const Sequence_Proxy_Iterator& other) const
    {
        return not _reversed? other.base() - this->base() : this->base() - other.base();
    }

    bool _reversed;
    bool _complemented;
}; // class Sequence_Proxy_Iterator

} // namespace detail

template < typename String_Type, typename Complementer = DNA_Complementer >
class Sequence_Proxy
{
private:
    typedef Sequence< String_Type, Complementer > sequence_type;
    typedef detail::Sequence_Proxy_Iterator< Complementer > iterator;

    friend class Sequence< String_Type, Complementer >;

    DEFAULT_COPY_CTOR(Sequence_Proxy);
    DELETE_COPY_ASOP(Sequence_Proxy);
    Sequence_Proxy(const String_Type* seq_p, size_t start, size_t len, bool reversed, bool complemented)
        : _seq_p(seq_p), _start(start), _len(len), _reversed(reversed), _complemented(complemented) {}

public:
    Sequence_Proxy() : _seq_p(), _start(0), _len(0), _reversed(false), _complemented(false) {}
    DEFAULT_MOVE_CTOR(Sequence_Proxy);
    DEFAULT_MOVE_ASOP(Sequence_Proxy);

    Sequence_Proxy(const String_Type& s) : _seq_p(&s), _start(0), _len(s.size()), _reversed(false), _complemented(false) {}

    size_t size() const { return _len; }
    bool empty() const { return size() > 0; }

    /** Further rev/comp/substr applied on this proxy */
    Sequence_Proxy rev(bool reversed = true) const
    { Sequence_Proxy res = *this; res._reversed = (_reversed != reversed); return res; }
    Sequence_Proxy comp(bool complemented = true) const
    { Sequence_Proxy res = *this; res._complemented = (_complemented != complemented); return res; }
    Sequence_Proxy revcomp(bool revcomped = true) const { return (*this).rev(revcomped).comp(revcomped); }
    Sequence_Proxy substr(size_t pos = 0, size_t len = String_Type::npos) const
    {
        Sequence_Proxy res = *this;
        if (pos > size())
        {
            throw std::out_of_range("pos > size");
        }
        if (len == String_Type::npos)
        {
            len = size() - pos;
        }
        else
        {
            if (pos + len > size())
            {
                throw std::out_of_range("pos + len > size");
            }
        }
        if (not _reversed)
        {
            res._start += pos;
        }
        else
        {
            res._start = (_start + _len - (pos + len));
        }
        res._len = len;
        return res;
    }

    /** Element access */
    char operator [] (size_t pos) const
    {
        return *(begin() + pos);
    }
    char at(size_t pos) const
    {
        if (pos >= _len)
        {
            throw std::out_of_range("pos >= size");
        }
        return (*this)[pos];
    }

    /** Iterators */
    iterator begin() const { return iterator(&(*_seq_p)[not _reversed? _start : _start + _len], _reversed, _complemented); }
    iterator end() const { return iterator(&(*_seq_p)[not _reversed? _start + _len : _start], _reversed, _complemented); }

    /** Convert back to string */
    operator String_Type () const
    {
        if (not _reversed and not _complemented)
        {
            return String_Type(&(*_seq_p)[_start], &(*_seq_p)[_start + _len]);
        }
        else
        {
            return String_Type(begin(), end());
        }
    }

    /** Formatted output operator */
    friend std::ostream& operator << (std::ostream& os, const Sequence_Proxy& rhs)
    {
        if (not rhs._reversed and not rhs._complemented)
        {
            std::copy(&(*rhs._seq_p)[rhs._start], &(*rhs._seq_p)[rhs._start + rhs._len],
                      std::ostream_iterator< char >(os));
        }
        else
        {
            std::copy(rhs.begin(), rhs.end(),
                      std::ostream_iterator< char >(os));
        }
        return os;
    }

    /** Lexicographical comparison */
    int compare(const Sequence_Proxy& rhs) const
    {
        if (not _reversed and not _complemented and not rhs._reversed and not rhs._complemented)
        {
            return _seq_p->compare(_start, _len, *rhs._seq_p, rhs._start, rhs._len);
        }
        else
        {
            auto it1 = begin();
            auto it1_end = end();
            auto it2 = rhs.begin();
            auto it2_end = rhs.end();
            while (it1 != it1_end and it2 != it2_end)
            {
                if (*it1 != *it2)
                {
                    return static_cast< int >(*it1) - static_cast< int >(*it2);
                }
                ++it1;
                ++it2;
            }
            if (it1 == it1_end and it2 == it2_end)
            {
                return 0;
            }
            else if (it1 == it1_end)
            {
                return -1;
            }
            else // it2 == it2_end
            {
                return 1;
            }
        }
    }
    friend bool operator == (const Sequence_Proxy& lhs, const Sequence_Proxy& rhs)
    {
        return lhs.compare(rhs) == 0;
    }
    friend bool operator != (const Sequence_Proxy& lhs, const Sequence_Proxy& rhs) { return !(lhs == rhs); }
    friend bool operator < (const Sequence_Proxy& lhs, const Sequence_Proxy& rhs)
    {
        return lhs.compare(rhs) < 0;
    }
    friend bool operator <= (const Sequence_Proxy& lhs, const Sequence_Proxy& rhs) { return lhs < rhs or lhs == rhs; }
    friend bool operator > (const Sequence_Proxy& lhs, const Sequence_Proxy& rhs)
    {
        return lhs.compare(rhs) > 0;
    }
    friend bool operator >= (const Sequence_Proxy& lhs, const Sequence_Proxy& rhs) { return lhs > rhs or lhs == rhs; }

private:
    const String_Type* _seq_p;
    size_t _start;
    size_t _len;
    bool _reversed;
    bool _complemented;
}; // class Sequence_Proxy

template < typename String_Type, typename Complementer = DNA_Complementer >
class Sequence
    : public String_Type
{
public:
    typedef Sequence_Proxy< String_Type, Complementer > proxy_type;

    DEFAULT_DEF_CTOR(Sequence);
    DEFAULT_COPY_CTOR(Sequence);
    DEFAULT_MOVE_CTOR(Sequence);
    DEFAULT_COPY_ASOP(Sequence);
    DEFAULT_MOVE_ASOP(Sequence);

    /** Construction from base */
    Sequence(const String_Type& seq) : String_Type(seq) {}
    Sequence(String_Type&& seq) : String_Type(std::move(seq)) {}
    Sequence(const proxy_type& seq) : String_Type(std::move(String_Type(seq))) {}
    /** Assignment from base */
    Sequence& operator = (const String_Type& seq)
    {
        if (&seq != static_cast< String_Type* >(this))
        {
            *static_cast< String_Type* >(this) = seq;
        }
        return *this;
    }
    Sequence& operator = (String_Type&& seq)
    {
        if (&seq != static_cast< String_Type* >(this))
        {
            *static_cast< String_Type* >(this) = std::move(seq);
        }
        return *this;
    }
    Sequence& operator = (const proxy_type& seq)
    {
        *this = std::move(String_Type(seq));
        return *this;
    }

    /** Functions producing proxies */
    proxy_type rev(bool reversed = true) const
    {
        return proxy_type(this, 0, this->size(), reversed, false);
    }
    proxy_type comp(bool complemented = true) const
    {
        return proxy_type(this, 0, this->size(), false, complemented);
    }
    proxy_type revcomp(bool revcomped = true) const
    {
        return proxy_type(this, 0, this->size(), revcomped, revcomped);
    }
    proxy_type substr(size_t pos = 0, size_t len = String_Type::npos) const
    {
        if (pos > this->size())
        {
            throw std::out_of_range("pos > size");
        }
        if (len == String_Type::npos)
        {
            len = this->size() - pos;
        }
        else
        {
            if (pos + len > this->size())
            {
                throw std::out_of_range("pos + len > size");
            }
        }
        return proxy_type(this, pos, len, false, false);
    }

    /** Functions using proxies */
    using String_Type::append;
    Sequence& append(const proxy_type& p)
    {
        if (not p._reversed and not p._complemented)
        {
            this->append(*p._seq_p, p._start, p._len);
        }
        else
        {
            this->append(p.begin(), p.end());
        }
        return *this;
    }
    using String_Type::operator +=;
    Sequence& operator += (const proxy_type& p) { return this->append(p); }

}; // class Sequence

} // namespace rc_sequence


#endif
