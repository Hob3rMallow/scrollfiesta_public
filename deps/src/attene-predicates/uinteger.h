#ifndef __UINTEGER_H_
#define __UINTEGER_H_

// Multi-precision unsigned integer backing BSNumber / BSRational, ported from
// GTE's UIntegerALU32 / UIntegerAP32. It is NOT a general-purpose bigint: words
// are 32-bit, little-endian, and the high word is always nonzero (num_bits
// tracks the leading 1-bit).
//
//   UIntegerAP -- arbitrary precision; words live in a core/buf.h stretchy
//                buffer (grows as needed).
//
// The arithmetic-logic-unit is a set of free functions (uint_add, uint_sub,
// ...) over UIntegerAP's small interface (num_bits/set_num_bits/size/words/
// get_back/set_back/set_all_zero, plus ctors from uint32_t/uint64_t).

#include <stdint.h>
#include <assert.h>
#include "bit_scan.h"

// =============================================================================
// UIntegerAP -- arbitrary precision (buf.h-backed; methods in uinteger.cpp)
// =============================================================================

struct UIntegerAP
{
    int32_t   m_num_bits;
    int32_t   m_size;
    uint32_t* m_words; // core/buf.h stretchy buffer (may be null when empty)

    UIntegerAP() : m_num_bits(0), m_size(0), m_words(nullptr) {}
    UIntegerAP(uint32_t number);
    UIntegerAP(uint64_t number);
    UIntegerAP(const UIntegerAP& number);
    UIntegerAP(UIntegerAP&& number);
    UIntegerAP& operator=(const UIntegerAP& number);
    UIntegerAP& operator=(UIntegerAP&& number);
    ~UIntegerAP();

    inline int32_t         num_bits() const { return m_num_bits; }
    inline int32_t         size()     const { return m_size; }
    inline uint32_t*       words()          { return m_words; }
    inline const uint32_t* words()    const { return m_words; }
    inline uint32_t        get_back() const { return m_words[m_size - 1]; }
    inline void            set_back(uint32_t v) { m_words[m_size - 1] = v; }

    void set_num_bits(int32_t numBits);
    void set_all_zero();
};

// =============================================================================
// Arithmetic-logic-unit (free functions over UIntegerAP). Comparisons assume
// the leading 1-bits are aligned (the BSNumber callers guarantee this). 'self'
// must not alias 'n0' or 'n1' (the callers guarantee this too).
// =============================================================================

inline bool uint_equal(const UIntegerAP& a, const UIntegerAP& b)
{
    int32_t numBits = a.num_bits();
    if (numBits != b.num_bits())
        return false;
    if (numBits > 0)
    {
        const uint32_t* ab = a.words();
        const uint32_t* bb = b.words();
        for (int32_t i = a.size() - 1; i >= 0; --i)
        {
            if (ab[i] != bb[i])
                return false;
        }
    }
    return true;
}

inline bool uint_less(const UIntegerAP& self, const UIntegerAP& number)
{
    int32_t nNumBits = number.num_bits();
    const uint32_t* nBits = number.words();

    int32_t numBits = self.num_bits();
    if (numBits > 0 && nNumBits > 0)
    {
        const uint32_t* bits = self.words();
        int32_t bitIndex0 = numBits - 1;
        int32_t bitIndex1 = nNumBits - 1;
        int32_t block0 = bitIndex0 / 32;
        int32_t block1 = bitIndex1 / 32;
        int32_t numBlockBits0 = 1 + (bitIndex0 % 32);
        int32_t numBlockBits1 = 1 + (bitIndex1 % 32);
        uint64_t n0shift = bits[block0];
        uint64_t n1shift = nBits[block1];
        while (block0 >= 0 && block1 >= 0)
        {
            uint32_t value0 = (uint32_t)((n0shift << (32 - numBlockBits0)) & 0xFFFFFFFFull);
            uint32_t value1 = (uint32_t)((n1shift << (32 - numBlockBits1)) & 0xFFFFFFFFull);

            if (--block0 >= 0)
            {
                n0shift = bits[block0];
                value0 |= (uint32_t)((n0shift >> numBlockBits0) & 0xFFFFFFFFull);
            }
            if (--block1 >= 0)
            {
                n1shift = nBits[block1];
                value1 |= (uint32_t)((n1shift >> numBlockBits1) & 0xFFFFFFFFull);
            }
            if (value0 < value1) return true;
            if (value0 > value1) return false;
        }
        return block0 < block1;
    }
    else
    {
        return nNumBits > 0;
    }
}

inline void uint_add(UIntegerAP& self, const UIntegerAP& n0, const UIntegerAP& n1)
{
    int32_t n0NumBits = n0.num_bits();
    int32_t n1NumBits = n1.num_bits();

    int32_t numBits = (n0NumBits >= n1NumBits ? n0NumBits : n1NumBits) + 1;
    self.set_num_bits(numBits);
    self.set_back(0);

    int32_t numElements0 = n0.size();
    int32_t numElements1 = n1.size();

    const uint32_t* u0 = (numElements0 >= numElements1 ? n0.words() : n1.words());
    const uint32_t* u1 = (numElements0 >= numElements1 ? n1.words() : n0.words());
    int32_t firstN  = (numElements0 <= numElements1 ? numElements0 : numElements1);
    int32_t secondN = (numElements0 >= numElements1 ? numElements0 : numElements1);

    uint32_t* bits = self.words();
    uint64_t carry = 0, sum;
    int32_t i;
    for (i = 0; i < firstN; ++i)
    {
        sum = (uint64_t)u0[i] + ((uint64_t)u1[i] + carry);
        bits[i] = (uint32_t)(sum & 0xFFFFFFFFull);
        carry = (sum >> 32);
    }
    if (carry > 0)
    {
        for (/**/; i < secondN; ++i)
        {
            sum = (uint64_t)u0[i] + carry;
            bits[i] = (uint32_t)(sum & 0xFFFFFFFFull);
            carry = (sum >> 32);
        }
        if (carry > 0)
            bits[i] = (uint32_t)(carry & 0xFFFFFFFFull);
    }
    else
    {
        for (/**/; i < secondN; ++i)
            bits[i] = u0[i];
    }

    uint32_t firstBitIndex = (numBits - 1) % 32;
    uint32_t mask = (1u << firstBitIndex);
    if ((mask & self.get_back()) == 0)
        self.set_num_bits(--numBits);
}

// Requires n0 > n1. Uses two's-complement of n1.
inline void uint_sub(UIntegerAP& self, const UIntegerAP& n0, const UIntegerAP& n1)
{
    int32_t n0NumBits = n0.num_bits();
    const uint32_t* n0Bits = n0.words();
    const uint32_t* n1Bits = n1.words();

    int32_t numElements0 = n0.size();
    int32_t numElements1 = n1.size();

    UIntegerAP n2;
    n2.set_num_bits(n0NumBits);
    uint32_t* n2Bits = n2.words();
    int32_t i;
    for (i = 0; i < numElements1; ++i)
        n2Bits[i] = ~n1Bits[i];
    for (/**/; i < numElements0; ++i)
        n2Bits[i] = ~0u;

    uint64_t carry = 1, sum;
    for (i = 0; i < numElements0; ++i)
    {
        sum = (uint64_t)n2Bits[i] + carry;
        n2Bits[i] = (uint32_t)(sum & 0xFFFFFFFFull);
        carry = (sum >> 32);
    }

    self.set_num_bits(n0NumBits + 1);
    self.set_back(0);

    uint32_t* bits = self.words();
    for (i = 0, carry = 0; i < numElements0; ++i)
    {
        sum = (uint64_t)n2Bits[i] + ((uint64_t)n0Bits[i] + carry);
        bits[i] = (uint32_t)(sum & 0xFFFFFFFFull);
        carry = (sum >> 32);
    }

    int32_t block;
    for (block = numElements0 - 1; block >= 0; --block)
    {
        if (bits[block] > 0)
            break;
    }

    if (block >= 0)
        self.set_num_bits(32 * block + bit_leading_u32(bits[block]) + 1);
    else
        self.set_num_bits(0);
}

inline void uint_mul(UIntegerAP& self, const UIntegerAP& n0, const UIntegerAP& n1)
{
    int32_t n0NumBits = n0.num_bits();
    int32_t n1NumBits = n1.num_bits();
    const uint32_t* n0Bits = n0.words();
    const uint32_t* n1Bits = n1.words();

    int32_t numBits = n0NumBits + n1NumBits;
    self.set_num_bits(numBits);
    self.set_all_zero();
    uint32_t* bits = self.words();

    UIntegerAP product;
    product.set_num_bits(numBits);
    product.set_all_zero();
    uint32_t* pBits = product.words();

    int32_t const numElements0 = n0.size();
    int32_t const numElements1 = n1.size();
    int32_t const numElements = self.size();

    int32_t i0, i1, i2;
    uint64_t term, sum;

    uint64_t block0 = n0Bits[0];
    uint64_t carry = 0;
    for (i1 = 0; i1 < numElements1; ++i1)
    {
        term = block0 * n1Bits[i1] + carry;
        bits[i1] = (uint32_t)(term & 0xFFFFFFFFull);
        carry = (term >> 32);
    }
    if (i1 < numElements)
        bits[i1] = (uint32_t)(carry & 0xFFFFFFFFull);

    for (i0 = 1; i0 < numElements0; ++i0)
    {
        block0 = n0Bits[i0];
        carry = 0;
        for (i1 = 0, i2 = i0; i1 < numElements1; ++i1, ++i2)
        {
            term = block0 * n1Bits[i1] + carry;
            pBits[i2] = (uint32_t)(term & 0xFFFFFFFFull);
            carry = (term >> 32);
        }
        if (i2 < numElements)
            pBits[i2] = (uint32_t)(carry & 0xFFFFFFFFull);

        carry = 0;
        for (i1 = 0, i2 = i0; i1 < numElements1; ++i1, ++i2)
        {
            sum = (uint64_t)pBits[i2] + ((uint64_t)bits[i2] + carry);
            bits[i2] = (uint32_t)(sum & 0xFFFFFFFFull);
            carry = (sum >> 32);
        }
        if (i2 < numElements)
        {
            sum = (uint64_t)pBits[i2] + carry;
            bits[i2] = (uint32_t)(sum & 0xFFFFFFFFull);
        }
    }

    uint32_t firstBitIndex = (numBits - 1) % 32;
    uint32_t mask = (1u << firstBitIndex);
    if ((mask & self.get_back()) == 0)
        self.set_num_bits(--numBits);
}

inline void uint_shift_left(UIntegerAP& self, const UIntegerAP& number, int32_t shift)
{
    int32_t nNumBits = number.num_bits();
    const uint32_t* nBits = number.words();

    self.set_num_bits(nNumBits + shift);

    uint32_t* bits = self.words();
    int32_t const shiftBlock = shift / 32;
    for (int32_t i = 0; i < shiftBlock; ++i)
        bits[i] = 0;

    int32_t const numInElements = number.size();
    int32_t const lshift = shift % 32;
    int32_t i, j;
    if (lshift > 0)
    {
        int32_t const rshift = 32 - lshift;
        uint32_t prev = 0, curr;
        for (i = shiftBlock, j = 0; j < numInElements; ++i, ++j)
        {
            curr = nBits[j];
            bits[i] = (curr << lshift) | (prev >> rshift);
            prev = curr;
        }
        if (i < self.size())
            bits[i] = (prev >> rshift);
    }
    else
    {
        for (i = shiftBlock, j = 0; j < numInElements; ++i, ++j)
            bits[i] = nBits[j];
    }
}

// 'number' is even and positive; shift it right to odd. Returns the shift count.
inline int32_t uint_shift_right_to_odd(UIntegerAP& self, const UIntegerAP& number)
{
    const uint32_t* nBits = number.words();

    int32_t const numElements = number.size();
    int32_t const numM1 = numElements - 1;
    int32_t firstBitIndex = 32 * numM1 + bit_leading_u32(nBits[numM1]);

    int32_t lastBitIndex = -1;
    for (int32_t block = 0; block < numElements; ++block)
    {
        uint32_t value = nBits[block];
        if (value > 0)
        {
            lastBitIndex = 32 * block + bit_trailing_u32(value);
            break;
        }
    }

    self.set_num_bits(firstBitIndex - lastBitIndex + 1);
    uint32_t* bits = self.words();
    int32_t const numBlocks = self.size();

    int32_t const shiftBlock = lastBitIndex / 32;
    int32_t rshift = lastBitIndex % 32;
    if (rshift > 0)
    {
        int32_t const lshift = 32 - rshift;
        int32_t i, j = shiftBlock;
        uint32_t curr = nBits[j++];
        for (i = 0; j < numElements; ++i, ++j)
        {
            uint32_t next = nBits[j];
            bits[i] = (curr >> rshift) | (next << lshift);
            curr = next;
        }
        if (i < numBlocks)
            bits[i] = (curr >> rshift);
    }
    else
    {
        for (int32_t i = 0, j = shiftBlock; i < numBlocks; ++i, ++j)
            bits[i] = nBits[j];
    }

    return rshift + 32 * shiftBlock;
}

// self = (self + 1) shifted right to odd; returns the shift count.
inline int32_t uint_round_up(UIntegerAP& self)
{
    UIntegerAP one(1u), rounded;
    uint_add(rounded, self, one);
    return uint_shift_right_to_odd(self, rounded);
}

// Leading numRequested (< 64) bits, left-aligned in the high bits of the result.
inline uint64_t uint_get_prefix(const UIntegerAP& self, int32_t numRequested)
{
    const uint32_t* bits = self.words();

    int32_t bitIndex = self.num_bits() - 1;
    int32_t blockIndex = bitIndex / 32;
    uint64_t prefix = bits[blockIndex];

    int32_t firstBitIndex = bitIndex % 32;
    int32_t numBlockBits = firstBitIndex + 1;

    int32_t targetIndex = 63;
    prefix <<= targetIndex - firstBitIndex;
    numRequested -= numBlockBits;
    targetIndex -= numBlockBits;

    if (numRequested > 0 && --blockIndex >= 0)
    {
        uint64_t nextBlock = bits[blockIndex];
        nextBlock <<= targetIndex - 31;
        prefix |= nextBlock;
        numRequested -= 32;
        targetIndex -= 32;

        if (numRequested > 0 && --blockIndex >= 0)
        {
            nextBlock = bits[blockIndex];
            nextBlock >>= 31 - targetIndex;
            prefix |= nextBlock;
        }
    }

    return prefix;
}

#endif
