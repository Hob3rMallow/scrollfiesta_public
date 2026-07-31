#ifndef __BSNUMBER_H_
#define __BSNUMBER_H_

// Exact arithmetic, ported from GTE's BSNumber.h / BSRational.h.
//
//   BSNumber   "binary scientific number" = sign * mantissa * 2^exponent, where
//              mantissa is an odd unsigned integer (trailing zero bits are
//              folded into the exponent). +, -, * are EXACT (no rounding).
//   BSRational = numerator / denominator, both BSNumber. +, -, *, / are EXACT.
//              Not reduced to lowest terms (only the power-of-two normalization
//              that BSNumber already does) -- matches GTE.
//
// The unsigned-integer backend is UIntegerAP (arbitrary precision). Float/double
// conversions use round-to-nearest-ties-to-even. This code is exact and finicky;
// see toaster/tests/test_bsnumber.cpp for the coverage that pins it down.

#include <stdint.h>
#include <string.h>
#include "uinteger.h"
#include "ieee_binary.h"

// Rounding modes for the precision/float conversions (subset of <cfenv>).
enum
{
    BS_ROUND_NEAREST    = 0, // round to nearest, ties to even
    BS_ROUND_UPWARD     = 1, // toward +infinity
    BS_ROUND_DOWNWARD   = 2, // toward -infinity
    BS_ROUND_TOWARDZERO = 3  // truncate
};

struct BSRational;

// =============================================================================
// BSNumber
// =============================================================================

struct BSNumber
{
    // 0 is { sign 0, biased_exponent 0, uinteger 0 }. Nonzero: sign != 0,
    // uinteger > 0 and odd.
    int32_t m_sign;
    int32_t m_biased_exponent;
    UIntegerAP m_uinteger;

    BSNumber() : m_sign(0), m_biased_exponent(0) {}

    BSNumber(int32_t number) { from_signed((int64_t)number); }
    BSNumber(int64_t number) { from_signed(number); }
    BSNumber(uint32_t number) { from_unsigned((uint64_t)number); }
    BSNumber(uint64_t number) { from_unsigned(number); }
    BSNumber(float number)  : m_sign(0), m_biased_exponent(0) { convert_from_f32(number); }
    BSNumber(double number) : m_sign(0), m_biased_exponent(0) { convert_from_f64(number); }
    BSNumber(const char* number) { from_string(number); }

    operator float()  const { return convert_to_f32(); }
    operator double() const { return convert_to_f64(); }

    // Member access.
    int32_t  sign()             const { return m_sign; }
    void     set_sign(int32_t s)      { m_sign = s; }
    int32_t  biased_exponent()  const { return m_biased_exponent; }
    void     set_biased_exponent(int32_t e) { m_biased_exponent = e; }
    int32_t  exponent()         const { return m_biased_exponent + m_uinteger.num_bits() - 1; }
    void     set_exponent(int32_t e)  { m_biased_exponent = e - m_uinteger.num_bits() + 1; }
    UIntegerAP const& uinteger()         const { return m_uinteger; }
    UIntegerAP&       uinteger()               { return m_uinteger; }

    // Comparisons.
    bool operator==(const BSNumber& n) const
    {
        return (m_sign == n.m_sign ? equal_ignore_sign(*this, n) : false);
    }
    bool operator!=(const BSNumber& n) const { return !operator==(n); }

    bool operator<(const BSNumber& n) const
    {
        if (m_sign > 0)
        {
            if (n.m_sign <= 0) return false;
            return less_than_ignore_sign(*this, n);
        }
        else if (m_sign < 0)
        {
            if (n.m_sign >= 0) return true;
            return less_than_ignore_sign(n, *this);
        }
        else
        {
            return n.m_sign > 0;
        }
    }
    bool operator<=(const BSNumber& n) const { return !n.operator<(*this); }
    bool operator>(const BSNumber& n)  const { return n.operator<(*this); }
    bool operator>=(const BSNumber& n) const { return !operator<(n); }

    // Unary.
    BSNumber operator+() const { return *this; }
    BSNumber operator-() const { BSNumber r = *this; r.m_sign = -r.m_sign; return r; }

    // Arithmetic.
    BSNumber operator+(const BSNumber& n1) const
    {
        const BSNumber& n0 = *this;
        if (n0.m_sign == 0) return n1;
        if (n1.m_sign == 0) return n0;

        if (n0.m_sign > 0)
        {
            if (n1.m_sign > 0)
                return add_ignore_sign(n0, n1, +1);
            if (!equal_ignore_sign(n0, n1))
                return less_than_ignore_sign(n1, n0) ? sub_ignore_sign(n0, n1, +1)
                                                     : sub_ignore_sign(n1, n0, -1);
        }
        else
        {
            if (n1.m_sign < 0)
                return add_ignore_sign(n0, n1, -1);
            if (!equal_ignore_sign(n0, n1))
                return less_than_ignore_sign(n1, n0) ? sub_ignore_sign(n0, n1, -1)
                                                     : sub_ignore_sign(n1, n0, +1);
        }
        return BSNumber(); // exact cancellation -> 0
    }

    BSNumber operator-(const BSNumber& n1) const { return operator+(-n1); }

    BSNumber operator*(const BSNumber& number) const
    {
        BSNumber result;
        int32_t sign = m_sign * number.m_sign;
        if (sign != 0)
        {
            result.m_sign = sign;
            result.m_biased_exponent = m_biased_exponent + number.m_biased_exponent;
            uint_mul(result.m_uinteger, m_uinteger, number.m_uinteger);
        }
        return result;
    }

    BSNumber& operator+=(const BSNumber& n) { *this = operator+(n); return *this; }
    BSNumber& operator-=(const BSNumber& n) { *this = operator-(n); return *this; }
    BSNumber& operator*=(const BSNumber& n) { *this = operator*(n); return *this; }

    // Serialization to/from a caller-owned byte blob. byte_size() gives the
    // exact number of bytes; write_bytes/read_bytes move them.
    int32_t byte_size() const
    {
        return 3 * (int32_t)sizeof(int32_t) + m_uinteger.size() * (int32_t)sizeof(uint32_t);
    }
    void write_bytes(uint8_t* dst) const
    {
        int32_t nb = m_uinteger.num_bits();
        memcpy(dst, &m_sign, sizeof(int32_t)); dst += sizeof(int32_t);
        memcpy(dst, &m_biased_exponent, sizeof(int32_t)); dst += sizeof(int32_t);
        memcpy(dst, &nb, sizeof(int32_t)); dst += sizeof(int32_t);
        memcpy(dst, m_uinteger.words(), m_uinteger.size() * sizeof(uint32_t));
    }
    void read_bytes(const uint8_t* src)
    {
        int32_t nb;
        memcpy(&m_sign, src, sizeof(int32_t)); src += sizeof(int32_t);
        memcpy(&m_biased_exponent, src, sizeof(int32_t)); src += sizeof(int32_t);
        memcpy(&nb, src, sizeof(int32_t)); src += sizeof(int32_t);
        m_uinteger.set_num_bits(nb);
        if (nb > 0)
            memcpy(m_uinteger.words(), src, m_uinteger.size() * sizeof(uint32_t));
    }

    // ---- internals -------------------------------------------------------

    void from_signed(int64_t number)
    {
        if (number == 0) { m_sign = 0; m_biased_exponent = 0; }
        else
        {
            uint64_t mag;
            if (number < 0) { m_sign = -1; mag = (0ull - (uint64_t)number); }
            else            { m_sign =  1; mag = (uint64_t)number; }
            m_biased_exponent = bit_trailing_u64(mag);
            m_uinteger = UIntegerAP(mag);
        }
    }

    void from_unsigned(uint64_t number)
    {
        if (number == 0) { m_sign = 0; m_biased_exponent = 0; }
        else
        {
            m_sign = 1;
            m_biased_exponent = bit_trailing_u64(number);
            m_uinteger = UIntegerAP(number);
        }
    }

    static bool equal_ignore_sign(const BSNumber& n0, const BSNumber& n1)
    {
        return n0.m_biased_exponent == n1.m_biased_exponent && uint_equal(n0.m_uinteger, n1.m_uinteger);
    }

    static bool less_than_ignore_sign(const BSNumber& n0, const BSNumber& n1)
    {
        int32_t e0 = n0.exponent(), e1 = n1.exponent();
        if (e0 < e1) return true;
        if (e0 > e1) return false;
        return uint_less(n0.m_uinteger, n1.m_uinteger);
    }

    static BSNumber add_ignore_sign(const BSNumber& n0, const BSNumber& n1, int32_t resultSign)
    {
        BSNumber result, temp;
        int32_t diff = n0.m_biased_exponent - n1.m_biased_exponent;
        if (diff > 0)
        {
            uint_shift_left(temp.m_uinteger, n0.m_uinteger, diff);
            uint_add(result.m_uinteger, temp.m_uinteger, n1.m_uinteger);
            result.m_biased_exponent = n1.m_biased_exponent;
        }
        else if (diff < 0)
        {
            uint_shift_left(temp.m_uinteger, n1.m_uinteger, -diff);
            uint_add(result.m_uinteger, n0.m_uinteger, temp.m_uinteger);
            result.m_biased_exponent = n0.m_biased_exponent;
        }
        else
        {
            uint_add(temp.m_uinteger, n0.m_uinteger, n1.m_uinteger);
            int32_t shift = uint_shift_right_to_odd(result.m_uinteger, temp.m_uinteger);
            result.m_biased_exponent = n0.m_biased_exponent + shift;
        }
        result.m_sign = resultSign;
        return result;
    }

    // n0 > n1 required.
    static BSNumber sub_ignore_sign(const BSNumber& n0, const BSNumber& n1, int32_t resultSign)
    {
        BSNumber result, temp;
        int32_t diff = n0.m_biased_exponent - n1.m_biased_exponent;
        if (diff > 0)
        {
            uint_shift_left(temp.m_uinteger, n0.m_uinteger, diff);
            uint_sub(result.m_uinteger, temp.m_uinteger, n1.m_uinteger);
            result.m_biased_exponent = n1.m_biased_exponent;
        }
        else if (diff < 0)
        {
            uint_shift_left(temp.m_uinteger, n1.m_uinteger, -diff);
            uint_sub(result.m_uinteger, n0.m_uinteger, temp.m_uinteger);
            result.m_biased_exponent = n0.m_biased_exponent;
        }
        else
        {
            uint_sub(temp.m_uinteger, n0.m_uinteger, n1.m_uinteger);
            int32_t shift = uint_shift_right_to_odd(result.m_uinteger, temp.m_uinteger);
            result.m_biased_exponent = n0.m_biased_exponent + shift;
        }
        result.m_sign = resultSign;
        return result;
    }

    static BSNumber convert_to_integer(const char* s)
    {
        size_t len = strlen(s);
        BSNumber x((int32_t)(s[len - 1] - '0'));
        if (len > 1)
        {
            BSNumber ten(10), pow10(10);
            for (size_t i = 1; i < len; ++i)
            {
                int32_t digit = (int32_t)(s[len - 1 - i] - '0');
                if (digit > 0)
                    x += BSNumber(digit) * pow10;
                pow10 = pow10 * ten;
            }
        }
        return x;
    }

    void from_string(const char* number)
    {
        int32_t sign = 1;
        const char* p = number;
        if (*p == '-') { sign = -1; ++p; }
        else if (*p == '+') { ++p; }
        *this = convert_to_integer(p);
        if (m_sign != 0)
            m_sign = sign;
    }

    // float -> BSNumber (IEEE binary32 decomposition).
    void convert_from_f32(float number)
    {
        IEEEBinary32 x(number);
        uint32_t s = x.GetSign();
        uint32_t e = x.GetBiased();
        uint32_t t = x.GetTrailing();

        if (e == 0)
        {
            if (t == 0) { m_sign = 0; m_biased_exponent = 0; }
            else
            {
                int32_t last = bit_trailing_u32(t);
                int32_t diff = IEEEBinary32::NUM_TRAILING_BITS - last;
                m_sign = (s > 0 ? -1 : 1);
                m_biased_exponent = IEEEBinary32::MIN_SUB_EXPONENT - diff;
                m_uinteger = UIntegerAP((uint32_t)(t >> last));
            }
        }
        else if ((int32_t)e < IEEEBinary32::MAX_BIASED_EXPONENT)
        {
            if (t > 0)
            {
                int32_t last = bit_trailing_u32(t);
                int32_t diff = IEEEBinary32::NUM_TRAILING_BITS - last;
                m_sign = (s > 0 ? -1 : 1);
                m_biased_exponent = (int32_t)e - IEEEBinary32::EXPONENT_BIAS - diff;
                m_uinteger = UIntegerAP((uint32_t)((t | IEEEBinary32::SUP_TRAILING) >> last));
            }
            else
            {
                m_sign = (s > 0 ? -1 : 1);
                m_biased_exponent = (int32_t)e - IEEEBinary32::EXPONENT_BIAS;
                m_uinteger = UIntegerAP((uint32_t)1);
            }
        }
        else
        {
            // Inf/NaN: BSNumber has no representation. Graceful fallback like GTE:
            // infinities -> +-2^(1+bias); NaN -> 0.
            if (t == 0)
            {
                m_sign = (s > 0 ? -1 : 1);
                m_biased_exponent = 1 + IEEEBinary32::EXPONENT_BIAS;
                m_uinteger = UIntegerAP((uint32_t)1);
            }
            else
            {
                m_sign = 0;
                m_biased_exponent = 0;
            }
        }
    }

    // double -> BSNumber (IEEE binary64 decomposition).
    void convert_from_f64(double number)
    {
        IEEEBinary64 x(number);
        uint64_t s = x.GetSign();
        uint64_t e = x.GetBiased();
        uint64_t t = x.GetTrailing();

        if (e == 0)
        {
            if (t == 0) { m_sign = 0; m_biased_exponent = 0; }
            else
            {
                int32_t last = bit_trailing_u64(t);
                int32_t diff = IEEEBinary64::NUM_TRAILING_BITS - last;
                m_sign = (s > 0 ? -1 : 1);
                m_biased_exponent = IEEEBinary64::MIN_SUB_EXPONENT - diff;
                m_uinteger = UIntegerAP((uint64_t)(t >> last));
            }
        }
        else if ((int32_t)e < IEEEBinary64::MAX_BIASED_EXPONENT)
        {
            if (t > 0)
            {
                int32_t last = bit_trailing_u64(t);
                int32_t diff = IEEEBinary64::NUM_TRAILING_BITS - last;
                m_sign = (s > 0 ? -1 : 1);
                m_biased_exponent = (int32_t)e - IEEEBinary64::EXPONENT_BIAS - diff;
                m_uinteger = UIntegerAP((uint64_t)((t | IEEEBinary64::SUP_TRAILING) >> last));
            }
            else
            {
                m_sign = (s > 0 ? -1 : 1);
                m_biased_exponent = (int32_t)e - IEEEBinary64::EXPONENT_BIAS;
                m_uinteger = UIntegerAP((uint64_t)1);
            }
        }
        else
        {
            // Inf/NaN: BSNumber has no representation. Graceful fallback like GTE:
            // infinities -> +-2^(1+bias); NaN -> 0.
            if (t == 0)
            {
                m_sign = (s > 0 ? -1 : 1);
                m_biased_exponent = 1 + IEEEBinary64::EXPONENT_BIAS;
                m_uinteger = UIntegerAP((uint64_t)1);
            }
            else
            {
                m_sign = 0;
                m_biased_exponent = 0;
            }
        }
    }

    // Rounded leading significand bits, computed in 64-bit and narrowed by the
    // caller. Shared by convert_to_f32/f64 (numSignificandBits = 24 or 53).
    uint64_t get_trailing(int32_t numSignificandBits, int32_t normal, int32_t sigma) const
    {
        int32_t const numRequested = numSignificandBits + normal;
        uint64_t prefix = uint_get_prefix(m_uinteger, numRequested);

        int32_t diff = numRequested - sigma;
        int32_t roundBitIndex = 64 - diff;
        uint64_t mask = (1ull << roundBitIndex);
        uint64_t round;
        if (prefix & mask)
        {
            if (m_uinteger.num_bits() == diff)
                round = (prefix & (mask << 1)) ? 1 : 0; // ties to even
            else
                round = 1;
        }
        else
        {
            round = 0;
        }

        uint64_t trailing = prefix >> (roundBitIndex + 1);
        trailing += round;
        return trailing;
    }

    // BSNumber -> float (IEEE binary32 composition).
    float convert_to_f32() const
    {
        uint32_t s = (m_sign < 0 ? 1 : 0);
        uint32_t e, t;

        if (m_sign != 0)
        {
            int32_t exp = exponent();
            if (exp < IEEEBinary32::MIN_EXPONENT)
            {
                if (exp < IEEEBinary32::MIN_EXPONENT - 1 || m_uinteger.num_bits() == 1) { e = 0; t = 0; }
                else { e = 0; t = 1; }
            }
            else if (exp < IEEEBinary32::MIN_SUB_EXPONENT)
            {
                t = (uint32_t)get_trailing(IEEEBinary32::NUM_SIGNIFICAND_BITS, 0, IEEEBinary32::MIN_SUB_EXPONENT - exp - 1);
                if (t & IEEEBinary32::SUP_TRAILING) { e = 1; t = 0; }
                else { e = 0; }
            }
            else if (exp <= IEEEBinary32::EXPONENT_BIAS)
            {
                e = (uint32_t)(exp + IEEEBinary32::EXPONENT_BIAS);
                t = (uint32_t)get_trailing(IEEEBinary32::NUM_SIGNIFICAND_BITS, 1, 0);
                if (t & (IEEEBinary32::SUP_TRAILING << 1)) { ++e; t >>= 1; }
                t &= ~IEEEBinary32::SUP_TRAILING;
            }
            else
            {
                e = (uint32_t)IEEEBinary32::MAX_BIASED_EXPONENT; // infinity
                t = 0;
            }
        }
        else { e = 0; t = 0; }

        IEEEBinary32 x(s, e, t);
        return x.number;
    }

    // BSNumber -> double (IEEE binary64 composition).
    double convert_to_f64() const
    {
        uint64_t s = (m_sign < 0 ? 1 : 0);
        uint64_t e, t;

        if (m_sign != 0)
        {
            int32_t exp = exponent();
            if (exp < IEEEBinary64::MIN_EXPONENT)
            {
                if (exp < IEEEBinary64::MIN_EXPONENT - 1 || m_uinteger.num_bits() == 1) { e = 0; t = 0; }
                else { e = 0; t = 1; }
            }
            else if (exp < IEEEBinary64::MIN_SUB_EXPONENT)
            {
                t = get_trailing(IEEEBinary64::NUM_SIGNIFICAND_BITS, 0, IEEEBinary64::MIN_SUB_EXPONENT - exp - 1);
                if (t & IEEEBinary64::SUP_TRAILING) { e = 1; t = 0; }
                else { e = 0; }
            }
            else if (exp <= IEEEBinary64::EXPONENT_BIAS)
            {
                e = (uint64_t)(exp + IEEEBinary64::EXPONENT_BIAS);
                t = get_trailing(IEEEBinary64::NUM_SIGNIFICAND_BITS, 1, 0);
                if (t & (IEEEBinary64::SUP_TRAILING << 1)) { ++e; t >>= 1; }
                t &= ~IEEEBinary64::SUP_TRAILING;
            }
            else
            {
                e = (uint64_t)IEEEBinary64::MAX_BIASED_EXPONENT; // infinity
                t = 0;
            }
        }
        else { e = 0; t = 0; }

        IEEEBinary64 x(s, e, t);
        return x.number;
    }
};

// Long division n/d -> a 'precision'-bit BSNumber, rounded per roundingMode.
// (GTE's Convert(BSRational, precision, mode, BSNumber&), taking the numerator
// and denominator directly to avoid a circular dependency on BSRational.)
inline void bsnumber_divide_round(const BSNumber& numerator, const BSNumber& denominator,
    int32_t precision, int32_t roundingMode, BSNumber& output)
{
    if (numerator.sign() == 0) { output = BSNumber(0); return; }

    BSNumber n = numerator;
    BSNumber d = denominator;

    int32_t sign = n.sign() * d.sign();
    n.set_sign(1);
    d.set_sign(1);
    int32_t pmq = n.exponent() - d.exponent();
    n.set_exponent(0);
    d.set_exponent(0);
    if (n < d) { n.set_exponent(n.exponent() + 1); --pmq; }

    BSNumber two(2);
    UIntegerAP& w = output.uinteger();
    w.set_num_bits(precision);
    w.set_all_zero();
    int32_t const size = w.size();
    int32_t const precisionM1 = precision - 1;
    int32_t leading = precisionM1 % 32;
    uint32_t mask = (1u << leading);
    uint32_t* bits = w.words();
    int32_t current = size - 1;
    int32_t lastBit = -1;
    for (int32_t i = precisionM1; i >= 0; --i)
    {
        if (n < d) { n = two * n; lastBit = 0; }
        else       { n = two * (n - d); bits[current] |= mask; lastBit = 1; }

        if (mask == 0x00000001u) { --current; mask = 0x80000000u; }
        else                     { mask >>= 1; }
    }

    if (roundingMode == BS_ROUND_NEAREST)
    {
        n = n - d;
        if (n.sign() > 0 || (n.sign() == 0 && lastBit == 1))
            pmq += uint_round_up(w);
    }
    else if (roundingMode == BS_ROUND_UPWARD)
    {
        if (n.sign() > 0 && sign > 0)
            pmq += uint_round_up(w);
    }
    else if (roundingMode == BS_ROUND_DOWNWARD)
    {
        if (n.sign() > 0 && sign < 0)
            pmq += uint_round_up(w);
    }
    // else BS_ROUND_TOWARDZERO: truncate (nothing to do)

    output.set_sign(sign);
    output.set_biased_exponent(pmq - precisionM1);
}

// =============================================================================
// BSRational
// =============================================================================

struct BSRational
{
    BSNumber m_numerator;
    BSNumber m_denominator;

    BSRational() : m_numerator(0), m_denominator(1) {}
    BSRational(int32_t n)  : m_numerator(n), m_denominator(1) {}
    BSRational(int64_t n)  : m_numerator(n), m_denominator(1) {}
    BSRational(uint32_t n) : m_numerator(n), m_denominator(1) {}
    BSRational(uint64_t n) : m_numerator(n), m_denominator(1) {}
    BSRational(float n)    : m_numerator(n), m_denominator(1) {}
    BSRational(double n)   : m_numerator(n), m_denominator(1) {}
    BSRational(const BSNumber& n) : m_numerator(n), m_denominator(1) {}

    BSRational(float n,   float d)   { init_fraction(BSNumber(n), BSNumber(d)); }
    BSRational(double n,  double d)  { init_fraction(BSNumber(n), BSNumber(d)); }
    BSRational(int32_t n, int32_t d) { init_fraction(BSNumber(n), BSNumber(d)); }
    BSRational(int64_t n, int64_t d) { init_fraction(BSNumber(n), BSNumber(d)); }
    BSRational(const BSNumber& n, const BSNumber& d) { init_fraction_normalized(n, d); }
    BSRational(const char* number) { from_string(number); }

    operator float()  const { BSNumber r; bsnumber_divide_round(m_numerator, m_denominator, 24, BS_ROUND_NEAREST, r); return (float)r; }
    operator double() const { BSNumber r; bsnumber_divide_round(m_numerator, m_denominator, 53, BS_ROUND_NEAREST, r); return (double)r; }

    int32_t sign() const { return m_numerator.sign() * m_denominator.sign(); }
    void    set_sign(int32_t s) { m_numerator.set_sign(s); m_denominator.set_sign(1); }
    BSNumber const& numerator()   const { return m_numerator; }
    BSNumber&       numerator()         { return m_numerator; }
    BSNumber const& denominator() const { return m_denominator; }
    BSNumber&       denominator()       { return m_denominator; }

    bool operator==(const BSRational& r) const
    {
        if (m_numerator.m_sign != r.m_numerator.m_sign) return false;
        if (m_numerator.m_sign == 0) return true;
        return m_numerator * r.m_denominator == m_denominator * r.m_numerator;
    }
    bool operator!=(const BSRational& r) const { return !operator==(r); }

    bool operator<(const BSRational& r) const
    {
        if (m_numerator.m_sign > 0)
        {
            if (r.m_numerator.m_sign <= 0) return false;
        }
        else if (m_numerator.m_sign == 0)
        {
            return r.m_numerator.m_sign > 0;
        }
        else // m_numerator.m_sign < 0
        {
            if (r.m_numerator.m_sign >= 0) return true;
        }
        return m_numerator * r.m_denominator < m_denominator * r.m_numerator;
    }
    bool operator<=(const BSRational& r) const { return !r.operator<(*this); }
    bool operator>(const BSRational& r)  const { return r.operator<(*this); }
    bool operator>=(const BSRational& r) const { return !operator<(r); }

    BSRational operator+() const { return *this; }
    BSRational operator-() const { return BSRational(-m_numerator, m_denominator); }

    BSRational operator+(const BSRational& r) const
    {
        BSNumber numerator = m_numerator * r.m_denominator + m_denominator * r.m_numerator;
        if (numerator.m_sign != 0)
            return BSRational(numerator, m_denominator * r.m_denominator);
        return BSRational(0);
    }

    BSRational operator-(const BSRational& r) const
    {
        BSNumber numerator = m_numerator * r.m_denominator - m_denominator * r.m_numerator;
        if (numerator.m_sign != 0)
            return BSRational(numerator, m_denominator * r.m_denominator);
        return BSRational(0);
    }

    BSRational operator*(const BSRational& r) const
    {
        BSNumber numerator = m_numerator * r.m_numerator;
        if (numerator.m_sign != 0)
            return BSRational(numerator, m_denominator * r.m_denominator);
        return BSRational(0);
    }

    BSRational operator/(const BSRational& r) const
    {
        // r != 0 required.
        BSNumber numerator = m_numerator * r.m_denominator;
        if (numerator.m_sign != 0)
        {
            BSNumber denominator = m_denominator * r.m_numerator;
            if (denominator.m_sign < 0) { numerator.m_sign = -numerator.m_sign; denominator.m_sign = 1; }
            return BSRational(numerator, denominator);
        }
        return BSRational(0);
    }

    BSRational& operator+=(const BSRational& r) { *this = operator+(r); return *this; }
    BSRational& operator-=(const BSRational& r) { *this = operator-(r); return *this; }
    BSRational& operator*=(const BSRational& r) { *this = operator*(r); return *this; }
    BSRational& operator/=(const BSRational& r) { *this = operator/(r); return *this; }

    // ---- internals -------------------------------------------------------

    void init_fraction(BSNumber n, BSNumber d)
    {
        if (d.m_sign < 0) { n.m_sign = -n.m_sign; d.m_sign = 1; }
        m_numerator = n;
        m_denominator = d;
    }

    void init_fraction_normalized(const BSNumber& nIn, const BSNumber& dIn)
    {
        BSNumber n = nIn, d = dIn;
        if (d.m_sign < 0) { n.m_sign = -n.m_sign; d.m_sign = 1; }
        // Drive the denominator exponent to a canonical position, compensating
        // on the numerator, so the two exponents do not grow together.
        n.m_biased_exponent -= d.exponent();
        d.m_biased_exponent = -(d.uinteger().num_bits() - 1);
        m_numerator = n;
        m_denominator = d;
    }

    static BSRational convert_to_fraction(const char* s)
    {
        BSRational y(0), ten(10), pow10(10);
        size_t len = strlen(s);
        for (size_t i = 0; i < len; ++i)
        {
            int32_t digit = (int32_t)(s[i] - '0');
            if (digit > 0)
                y += BSRational(digit) / pow10;
            pow10 = pow10 * ten;
        }
        return y;
    }

    void from_string(const char* number)
    {
        int32_t sign = 1;
        const char* p = number;
        if (*p == '-') { sign = -1; ++p; }
        else if (*p == '+') { ++p; }

        const char* dot = strchr(p, '.');
        if (dot)
        {
            size_t intLen = (size_t)(dot - p);
            const char* frac = dot + 1;
            if (intLen > 0)
            {
                // "x.y" or "x."
                char intbuf[512];
                memcpy(intbuf, p, intLen);
                intbuf[intLen] = '\0';
                BSRational intPart(BSNumber::convert_to_integer(intbuf));
                if (*frac != '\0')
                    *this = intPart + convert_to_fraction(frac);
                else
                    *this = intPart;
            }
            else
            {
                // ".y"
                *this = convert_to_fraction(frac);
            }
        }
        else
        {
            *this = BSRational(BSNumber::convert_to_integer(p));
        }
        // The parsed magnitude is nonnegative; apply the leading sign.
        if (m_numerator.m_sign != 0)
            m_numerator.m_sign = sign;
    }
};

#endif
