#define BUF_STANDALONE
#define BUF_IMPLEMENTATION
#include "uinteger.h"
#include "buf.h"
#include <string.h>

// UIntegerAP stores its 32-bit words in a core/buf.h stretchy buffer. m_size is
// kept explicitly (rather than read from buf_len) so uinteger.h's inline
// accessors need not pull in buf.h. NOTE: buf_resize does not zero new slots,
// unlike std::vector::resize that GTE relied on -- uint_mul compensates with
// set_all_zero(); every other ALU op writes all active words explicitly.

UIntegerAP::UIntegerAP(uint32_t number)
    : m_num_bits(0), m_size(0), m_words(nullptr)
{
    if (number > 0)
    {
        int32_t first = bit_leading_u32(number);
        int32_t last  = bit_trailing_u32(number);
        set_num_bits(first - last + 1);
        m_words[0] = (number >> last);
    }
}

UIntegerAP::UIntegerAP(uint64_t number)
    : m_num_bits(0), m_size(0), m_words(nullptr)
{
    if (number > 0)
    {
        int32_t first = bit_leading_u64(number);
        int32_t last  = bit_trailing_u64(number);
        number >>= last;
        set_num_bits(first - last + 1);
        m_words[0] = (uint32_t)(number & 0xFFFFFFFFull);
        if (m_size > 1)
            m_words[1] = (uint32_t)((number >> 32) & 0xFFFFFFFFull);
    }
}

UIntegerAP::UIntegerAP(const UIntegerAP& number)
    : m_num_bits(0), m_size(0), m_words(nullptr)
{
    *this = number;
}

UIntegerAP::UIntegerAP(UIntegerAP&& number)
    : m_num_bits(number.m_num_bits), m_size(number.m_size), m_words(number.m_words)
{
    number.m_num_bits = 0;
    number.m_size = 0;
    number.m_words = nullptr;
}

UIntegerAP& UIntegerAP::operator=(const UIntegerAP& number)
{
    if (this == &number)
        return *this;

    m_num_bits = number.m_num_bits;
    m_size = number.m_size;
    if (m_size > 0)
    {
        buf_resize(m_words, (size_t)m_size);
        memcpy(m_words, number.m_words, (size_t)m_size * sizeof(uint32_t));
    }
    else
    {
        buf_clear(m_words);
    }
    return *this;
}

UIntegerAP& UIntegerAP::operator=(UIntegerAP&& number)
{
    if (this == &number)
        return *this;

    buf_free(m_words);
    m_num_bits = number.m_num_bits;
    m_size = number.m_size;
    m_words = number.m_words;
    number.m_num_bits = 0;
    number.m_size = 0;
    number.m_words = nullptr;
    return *this;
}

UIntegerAP::~UIntegerAP()
{
    buf_free(m_words);
}

void UIntegerAP::set_num_bits(int32_t numBits)
{
    if (numBits > 0)
    {
        m_num_bits = numBits;
        m_size = 1 + (numBits - 1) / 32;
        buf_resize(m_words, (size_t)m_size);
    }
    else
    {
        m_num_bits = 0;
        m_size = 0;
        buf_clear(m_words);
    }
}

void UIntegerAP::set_all_zero()
{
    if (m_words && m_size > 0)
        memset(m_words, 0, (size_t)m_size * sizeof(uint32_t));
}
