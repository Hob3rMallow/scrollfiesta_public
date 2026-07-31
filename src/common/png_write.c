#include "png_write.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- CRC-32 (PNG / zlib polynomial) --------------------------------------- */
static uint32_t g_crc[256];
static int g_crc_init = 0;
static void crc_make(void)
{
    for (uint32_t n = 0; n < 256; n++) {
        uint32_t c = n;
        for (int k = 0; k < 8; k++) c = (c & 1) ? (0xedb88320u ^ (c >> 1)) : (c >> 1);
        g_crc[n] = c;
    }
    g_crc_init = 1;
}
static uint32_t crc32_2(const uint8_t *a, size_t na, const uint8_t *b, size_t nb)
{
    if (!g_crc_init) crc_make();
    uint32_t c = 0xffffffffu;
    for (size_t i = 0; i < na; i++) c = g_crc[(c ^ a[i]) & 0xff] ^ (c >> 8);
    for (size_t i = 0; i < nb; i++) c = g_crc[(c ^ b[i]) & 0xff] ^ (c >> 8);
    return c ^ 0xffffffffu;
}

static uint32_t adler32(const uint8_t *d, size_t n)
{
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < n; i++) { a = (a + d[i]) % 65521u; b = (b + a) % 65521u; }
    return (b << 16) | a;
}

static void put_be32(FILE *f, uint32_t v)
{
    uint8_t b[4] = { (uint8_t)(v>>24), (uint8_t)(v>>16), (uint8_t)(v>>8), (uint8_t)v };
    fwrite(b, 1, 4, f);
}

static void put_chunk(FILE *f, const char *type, const uint8_t *data, size_t len)
{
    put_be32(f, (uint32_t)len);
    fwrite(type, 1, 4, f);
    if (len) fwrite(data, 1, len, f);
    put_be32(f, crc32_2((const uint8_t *)type, 4, data, len));
}

int PngWrite_rgb(const char *path, int w, int h, const uint8_t *rgb)
{
    if (!path || !rgb || w <= 0 || h <= 0) return -1;
    FILE *f = fopen(path, "wb");
    if (!f) return -1;

    static const uint8_t sig[8] = { 0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a };
    fwrite(sig, 1, 8, f);

    /* IHDR */
    uint8_t ihdr[13];
    ihdr[0]=(uint8_t)(w>>24); ihdr[1]=(uint8_t)(w>>16); ihdr[2]=(uint8_t)(w>>8); ihdr[3]=(uint8_t)w;
    ihdr[4]=(uint8_t)(h>>24); ihdr[5]=(uint8_t)(h>>16); ihdr[6]=(uint8_t)(h>>8); ihdr[7]=(uint8_t)h;
    ihdr[8]=8;   /* bit depth */
    ihdr[9]=2;   /* color type: truecolor RGB */
    ihdr[10]=0;  /* compression */
    ihdr[11]=0;  /* filter */
    ihdr[12]=0;  /* interlace */
    put_chunk(f, "IHDR", ihdr, 13);

    /* raw filtered scanlines: filter byte 0 + RGB row */
    size_t stride = (size_t)w*3 + 1;
    size_t rawlen = stride * (size_t)h;
    uint8_t *raw = (uint8_t *)malloc(rawlen);
    if (!raw) { fclose(f); return -1; }
    for (int y = 0; y < h; y++) {
        raw[(size_t)y*stride] = 0;
        memcpy(raw + (size_t)y*stride + 1, rgb + (size_t)y*(size_t)w*3, (size_t)w*3);
    }

    /* zlib stream: 2-byte header + stored DEFLATE blocks + adler32 */
    size_t nblocks = rawlen / 65535u + 1;
    size_t zcap = 2 + rawlen + 5*nblocks + 4;
    uint8_t *z = (uint8_t *)malloc(zcap);
    if (!z) { free(raw); fclose(f); return -1; }
    size_t zp = 0;
    z[zp++] = 0x78; z[zp++] = 0x01;            /* CMF, FLG */
    size_t off = 0;
    while (off < rawlen) {
        size_t blk = rawlen - off; if (blk > 65535u) blk = 65535u;
        int final = (off + blk >= rawlen) ? 1 : 0;
        z[zp++] = (uint8_t)final;              /* BFINAL + BTYPE(00) */
        z[zp++] = (uint8_t)(blk & 0xff); z[zp++] = (uint8_t)((blk >> 8) & 0xff);
        uint16_t nlen = (uint16_t)~(uint16_t)blk;
        z[zp++] = (uint8_t)(nlen & 0xff); z[zp++] = (uint8_t)((nlen >> 8) & 0xff);
        memcpy(z + zp, raw + off, blk); zp += blk; off += blk;
    }
    uint32_t ad = adler32(raw, rawlen);
    z[zp++]=(uint8_t)(ad>>24); z[zp++]=(uint8_t)(ad>>16); z[zp++]=(uint8_t)(ad>>8); z[zp++]=(uint8_t)ad;
    put_chunk(f, "IDAT", z, zp);

    put_chunk(f, "IEND", NULL, 0);

    free(z); free(raw); fclose(f);
    return 0;
}
