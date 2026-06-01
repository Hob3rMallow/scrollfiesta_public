/*
 * tiffconf.h — Pre-configured for Windows x64 (MSVC).
 * Generated manually for vesuvius-c pipeline.
 * Only uncompressed TIFF I/O is needed.
 */

#ifndef _TIFFCONF_
#define _TIFFCONF_

#include <stddef.h>
#include <stdint.h>
#include <inttypes.h>

/* Signed types */
#define TIFF_INT8_T   int8_t
#define TIFF_INT16_T  int16_t
#define TIFF_INT32_T  int32_t
#define TIFF_INT64_T  int64_t

/* Unsigned types */
#define TIFF_UINT8_T  uint8_t
#define TIFF_UINT16_T uint16_t
#define TIFF_UINT32_T uint32_t
#define TIFF_UINT64_T uint64_t

/* Signed size type */
#ifdef _WIN64
#define TIFF_SSIZE_T  int64_t
#else
#define TIFF_SSIZE_T  int32_t
#endif

/* IEEE floating point */
#define HAVE_IEEEFP 1

/* Host fill order (x86 = LSB2MSB) */
#define HOST_FILLORDER FILLORDER_LSB2MSB

/* Little-endian (Intel/AMD) */
#define HOST_BIGENDIAN 0

/* Compression support — we only need uncompressed TIFF,
   but enable common ones for compatibility */
#define CCITT_SUPPORT    1
#define LOGLUV_SUPPORT   1
#define LZW_SUPPORT      1
#define NEXT_SUPPORT     1
#define PACKBITS_SUPPORT 1
#define THUNDER_SUPPORT  1

/* Strip chopping */
#define STRIPCHOP_DEFAULT TIFF_STRIPCHOP

/* SubIFD support */
#define SUBIFD_SUPPORT 1

/* Treat extra sample as alpha */
#define DEFAULT_EXTRASAMPLE_AS_ALPHA 1

/* MDI support */
#define MDI_SUPPORT 1

/* Backward compat macros */
#define COLORIMETRY_SUPPORT
#define YCBCR_SUPPORT
#define CMYK_SUPPORT
#define ICC_SUPPORT
#define PHOTOSHOP_SUPPORT
#define IPTC_SUPPORT

#endif /* _TIFFCONF_ */
