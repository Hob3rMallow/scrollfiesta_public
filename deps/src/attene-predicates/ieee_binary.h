#ifndef __IEEE_BINARY_H_
#define __IEEE_BINARY_H_

// IEEE-754 bit-level helpers. Ported from GTE's IEEEBinary.h / IEEEBinary16.h,
// reduced to the concrete float/double/half cases (no 4-parameter template).
//
//   IEEEBinary32 / IEEEBinary64 decompose a float/double into sign, biased
//   exponent and trailing significand; BSNumber uses them for exact
//   float<->BSNumber conversion.
//
//   IEEEBinary16 is a standalone binary16 (half-precision) value with
//   conversions to/from float and double.
//
// Each value type stores its bits and its float view in a union, so the same
// object can be inspected as raw bits or used arithmetically.

#include <stdint.h>

struct IEEEBinary32
{
    typedef float    FloatType;
    typedef uint32_t UIntType;

    static const int32_t  NUM_ENCODING_BITS    = 32;
    static const int32_t  NUM_SIGNIFICAND_BITS = 24;
    static const int32_t  NUM_TRAILING_BITS    = 23;
    static const int32_t  NUM_EXPONENT_BITS    = 8;
    static const int32_t  EXPONENT_BIAS        = 127;
    static const int32_t  MAX_BIASED_EXPONENT  = 255;
    static const int32_t  MIN_SUB_EXPONENT     = 1 - EXPONENT_BIAS;                    // -126
    static const int32_t  MIN_EXPONENT         = MIN_SUB_EXPONENT - NUM_TRAILING_BITS; // -149
    static const int32_t  SIGN_SHIFT           = 31;
    static const uint32_t SIGN_MASK     = 0x80000000u;
    static const uint32_t NOT_SIGN_MASK = 0x7FFFFFFFu;
    static const uint32_t TRAILING_MASK = 0x007FFFFFu;
    static const uint32_t EXPONENT_MASK = 0x7F800000u;
    static const uint32_t SUP_TRAILING  = 0x00800000u; // 1 << NUM_TRAILING_BITS
    static const uint32_t MAX_TRAILING  = 0x007FFFFFu;

    union { uint32_t encoding; float number; };

    IEEEBinary32() { encoding = 0; }
    IEEEBinary32(uint32_t inEncoding) { encoding = inEncoding; }
    IEEEBinary32(float inNumber) { number = inNumber; }
    IEEEBinary32(uint32_t sign, uint32_t biased, uint32_t trailing)
    {
        encoding = (sign << SIGN_SHIFT) | (biased << NUM_TRAILING_BITS) | trailing;
    }

    uint32_t GetSign()     const { return (encoding & SIGN_MASK) >> SIGN_SHIFT; }
    uint32_t GetBiased()   const { return (encoding & EXPONENT_MASK) >> NUM_TRAILING_BITS; }
    uint32_t GetTrailing() const { return encoding & TRAILING_MASK; }
};

struct IEEEBinary64
{
    typedef double   FloatType;
    typedef uint64_t UIntType;

    static const int32_t  NUM_ENCODING_BITS    = 64;
    static const int32_t  NUM_SIGNIFICAND_BITS = 53;
    static const int32_t  NUM_TRAILING_BITS    = 52;
    static const int32_t  NUM_EXPONENT_BITS    = 11;
    static const int32_t  EXPONENT_BIAS        = 1023;
    static const int32_t  MAX_BIASED_EXPONENT  = 2047;
    static const int32_t  MIN_SUB_EXPONENT     = 1 - EXPONENT_BIAS;                    // -1022
    static const int32_t  MIN_EXPONENT         = MIN_SUB_EXPONENT - NUM_TRAILING_BITS; // -1074
    static const int32_t  SIGN_SHIFT           = 63;
    static const uint64_t SIGN_MASK     = 0x8000000000000000ull;
    static const uint64_t NOT_SIGN_MASK = 0x7FFFFFFFFFFFFFFFull;
    static const uint64_t TRAILING_MASK = 0x000FFFFFFFFFFFFFull;
    static const uint64_t EXPONENT_MASK = 0x7FF0000000000000ull;
    static const uint64_t SUP_TRAILING  = 0x0010000000000000ull; // 1 << NUM_TRAILING_BITS
    static const uint64_t MAX_TRAILING  = 0x000FFFFFFFFFFFFFull;

    union { uint64_t encoding; double number; };

    IEEEBinary64() { encoding = 0; }
    IEEEBinary64(uint64_t inEncoding) { encoding = inEncoding; }
    IEEEBinary64(double inNumber) { number = inNumber; }
    IEEEBinary64(uint64_t sign, uint64_t biased, uint64_t trailing)
    {
        encoding = (sign << SIGN_SHIFT) | (biased << NUM_TRAILING_BITS) | trailing;
    }

    uint64_t GetSign()     const { return (encoding & SIGN_MASK) >> SIGN_SHIFT; }
    uint64_t GetBiased()   const { return (encoding & EXPONENT_MASK) >> NUM_TRAILING_BITS; }
    uint64_t GetTrailing() const { return encoding & TRAILING_MASK; }
};

// binary16 half-precision float. Conversions round-to-nearest-ties-to-even and
// go through binary32 (a half is always exactly representable as a float).
struct IEEEBinary16
{
    static const int32_t  NUM_TRAILING_BITS    = 10;
    static const int32_t  NUM_SIGNIFICAND_BITS = 11;
    static const int32_t  EXPONENT_BIAS        = 15;
    static const int32_t  MAX_BIASED_EXPONENT  = 31;
    static const uint16_t SIGN_MASK     = 0x8000;
    static const uint16_t NOT_SIGN_MASK = 0x7FFF;
    static const uint16_t TRAILING_MASK = 0x03FF;
    static const uint16_t EXPONENT_MASK = 0x7C00;
    static const uint16_t SUP_TRAILING  = 0x0400; // 1 << NUM_TRAILING_BITS
    static const uint16_t MIN_SUBNORMAL = 0x0001;
    static const uint16_t MAX_NORMAL    = 0x7BFF;
    static const uint16_t POS_INFINITY  = 0x7C00;
    static const uint16_t NEG_INFINITY  = 0xFC00;

    uint16_t encoding;

    IEEEBinary16() { encoding = 0; }
    IEEEBinary16(uint16_t inEncoding) { encoding = inEncoding; }
    IEEEBinary16(float inNumber);
    IEEEBinary16(double inNumber);

    operator float() const;
    operator double() const;

    uint16_t GetSign()     const { return (encoding & SIGN_MASK) >> 15; }
    uint16_t GetBiased()   const { return (encoding & EXPONENT_MASK) >> NUM_TRAILING_BITS; }
    uint16_t GetTrailing() const { return encoding & TRAILING_MASK; }

    bool operator==(IEEEBinary16 b) const { return (float)*this == (float)b; }
    bool operator!=(IEEEBinary16 b) const { return (float)*this != (float)b; }
    bool operator< (IEEEBinary16 b) const { return (float)*this <  (float)b; }
    bool operator<=(IEEEBinary16 b) const { return (float)*this <= (float)b; }
    bool operator> (IEEEBinary16 b) const { return (float)*this >  (float)b; }
    bool operator>=(IEEEBinary16 b) const { return (float)*this >= (float)b; }
};

#endif
