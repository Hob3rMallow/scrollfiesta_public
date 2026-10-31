/* mesh_chart_index -- per-chart index of a chart-contiguous VESMESH1 file (see mesh_chart_index.h). */
#include "mesh_chart_index.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _WIN32
#include <windows.h>
#endif

#include "ves_platform.h"

static const unsigned char mci_magic[8] = { 'V', 'M', 'C', 'H', 'I', 'D', 'X', '1' };
enum { MCI_VERSION = 1 };
#define MCI_BLOCK ((size_t)1 << 20)          /* vertices or faces per streamed block */
#define MCI_NONE UINT64_MAX

static int mci_seek(FILE *f, uint64_t off)
{
#ifdef _WIN32
    return off > (uint64_t)INT64_MAX ? -1 : _fseeki64(f, (__int64)off, SEEK_SET);
#else
    return fseeko(f, (off_t)off, SEEK_SET);
#endif
}

static int mci_read(FILE *f, void *buf, size_t n)
{
    return n == 0 || fread(buf, 1, n, f) == n ? 0 : -1;
}

static int mci_write(FILE *f, const void *buf, size_t n)
{
    return n == 0 || fwrite(buf, 1, n, f) == n ? 0 : -1;
}

static int mci_file_info(const char *path, uint64_t *bytes, int64_t *mtime)
{
#ifdef _WIN32
    struct _stat64 st;
    if (_stat64(path, &st) != 0) return -1;
#else
    struct stat st;
    if (stat(path, &st) != 0) return -1;
#endif
    if (bytes != NULL) *bytes = (uint64_t)st.st_size;
    if (mtime != NULL) *mtime = (int64_t)st.st_mtime;
    return 0;
}

static int mci_replace(const char *tmp, const char *dst)
{
#ifdef _WIN32
    return MoveFileExA(tmp, dst, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) ? 0 : -1;
#else
    return rename(tmp, dst);
#endif
}

int MeshChartIndex_path(const char *vmesh, char *out, size_t cap)
{
    int k = 0;
    if (vmesh == NULL || out == NULL || cap == 0) return -1;
    k = snprintf(out, cap, "%s.chartidx", vmesh);
    return k < 0 || (size_t)k >= cap ? -1 : 0;
}

static void mci_rec_init(MeshChartRec *r, int32_t chart, uint64_t v0)
{
    int k = 0;
    memset(r, 0, sizeof *r);
    r->chart = chart;
    r->v0 = v0;
    r->v1 = v0;
    r->f0 = MCI_NONE;
    r->f1 = MCI_NONE;
    for (k = 0; k < 3; k++) { r->lo[k] = INFINITY; r->hi[k] = -INFINITY; }
    for (k = 0; k < 2; k++) { r->uvlo[k] = INFINITY; r->uvhi[k] = -INFINITY; }
}

/* the record holding vertex v, from `hint` (faces are chart-ordered: the hint almost always holds) */
static size_t mci_find(const MeshChartRec *rec, size_t n, uint64_t v, size_t hint)
{
    size_t lo = 0, hi = n;
    if (hint < n && v >= rec[hint].v0 && v < rec[hint].v1) return hint;
    if (hint + 1 < n && v >= rec[hint + 1].v0 && v < rec[hint + 1].v1) return hint + 1;
    while (lo < hi) {                               /* the last record with v0 <= v */
        size_t mid = lo + (hi - lo) / 2;
        if (rec[mid].v0 <= v) lo = mid + 1;
        else hi = mid;
    }
    return lo == 0 ? n : lo - 1;
}

static int mci_write_index(const char *path, const MeshChartRec *rec, size_t n, const MeshBinHeader *h,
                           uint64_t chart_bytes, int64_t mtime)
{
    char tmp[4096];
    FILE *f = NULL;
    uint32_t version = MCI_VERSION, flags = h->flags;
    uint64_t nrec = (uint64_t)n;
    size_t i = 0;
    int ok = 0;
    /* the writer's own temporary (two close-ups may index the same mesh at once; the rename publishes either) */
    if (snprintf(tmp, sizeof tmp, "%s.partial.%d", path, ves_getpid()) >= (int)sizeof tmp) return -1;
    f = fopen(tmp, "wb");
    if (f == NULL) return -1;
    ok = mci_write(f, mci_magic, sizeof mci_magic) == 0 && mci_write(f, &version, sizeof version) == 0 &&
         mci_write(f, &flags, sizeof flags) == 0 && mci_write(f, &h->nv, sizeof h->nv) == 0 &&
         mci_write(f, &h->nf, sizeof h->nf) == 0 && mci_write(f, &h->file_bytes, sizeof h->file_bytes) == 0 &&
         mci_write(f, &chart_bytes, sizeof chart_bytes) == 0 && mci_write(f, &mtime, sizeof mtime) == 0 &&
         mci_write(f, &nrec, sizeof nrec) == 0;
    for (i = 0; ok && i < n; i++) {
        int32_t pad = 0;
        ok = mci_write(f, &rec[i].chart, sizeof rec[i].chart) == 0 && mci_write(f, &pad, sizeof pad) == 0 &&
             mci_write(f, &rec[i].v0, sizeof rec[i].v0) == 0 && mci_write(f, &rec[i].v1, sizeof rec[i].v1) == 0 &&
             mci_write(f, &rec[i].f0, sizeof rec[i].f0) == 0 && mci_write(f, &rec[i].f1, sizeof rec[i].f1) == 0 &&
             mci_write(f, rec[i].lo, sizeof rec[i].lo) == 0 && mci_write(f, rec[i].hi, sizeof rec[i].hi) == 0 &&
             mci_write(f, rec[i].uvlo, sizeof rec[i].uvlo) == 0 && mci_write(f, rec[i].uvhi, sizeof rec[i].uvhi) == 0;
    }
    if (fclose(f) != 0) ok = 0;
    if (!ok || mci_replace(tmp, path) != 0) { remove(tmp); return -1; }
    return 0;
}

/* The builder shared by the streaming and the in-memory index: the vertices with their charts, then the UVs, then the
 * faces, each fed in file order in any number of blocks.  A step returns 0, 1 (not chart-contiguous) or -1 (no memory,
 * bad input). */
typedef struct MciBuild {
    MeshChartRec *rec;
    size_t n, cap, r, cur, seen_n;
    uint8_t *seen;
    uint64_t nv, nuv, nf;           /* vertices, UVs and faces fed so far */
} MciBuild;

static void mci_build_free(MciBuild *s)
{
    free(s->rec);
    free(s->seen);
    memset(s, 0, sizeof *s);
}

static int mci_build_verts(MciBuild *s, const float *v, const int32_t *ch, size_t cnt)
{
    size_t k = 0;
    for (k = 0; k < cnt; k++) {
        uint64_t vi = s->nv + k;
        int32_t c = ch[k];
        int a = 0;
        MeshChartRec *q = NULL;
        if (s->n == 0 || s->rec[s->n - 1].chart != c) {
            if (s->n > 0) s->rec[s->n - 1].v1 = vi;
            if (c >= 0) {                          /* a chart id owns one run */
                if ((size_t)c >= s->seen_n) {
                    size_t nn = s->seen_n ? s->seen_n : 1024;
                    uint8_t *s2 = NULL;
                    while (nn <= (size_t)c) nn *= 2;
                    s2 = (uint8_t *)realloc(s->seen, nn);
                    if (s2 == NULL) return -1;
                    memset(s2 + s->seen_n, 0, nn - s->seen_n);
                    s->seen = s2;
                    s->seen_n = nn;
                }
                if (s->seen[c]) return 1;
                s->seen[c] = 1;
            }
            if (s->n == s->cap) {
                size_t nc = s->cap ? s->cap * 2 : 4096;
                MeshChartRec *r2 = (MeshChartRec *)realloc(s->rec, nc * sizeof *s->rec);
                if (r2 == NULL) return -1;
                s->rec = r2;
                s->cap = nc;
            }
            mci_rec_init(&s->rec[s->n++], c, vi);
        }
        q = &s->rec[s->n - 1];
        for (a = 0; a < 3; a++) {
            double x = (double)v[k * 3 + (size_t)a];
            if (x < q->lo[a]) q->lo[a] = x;
            if (x > q->hi[a]) q->hi[a] = x;
        }
    }
    s->nv += cnt;
    return 0;
}

/* after the last vertex: close the last run */
static int mci_build_verts_end(MciBuild *s)
{
    if (s->n == 0) return -1;
    s->rec[s->n - 1].v1 = s->nv;
    s->r = 0;
    return 0;
}

/* `cnt` UV pairs, float32 or float64 (f64), packed */
static int mci_build_uv(MciBuild *s, const void *uv, int f64, size_t cnt)
{
    const size_t es = f64 ? 2 * sizeof(double) : 2 * sizeof(float);
    const unsigned char *ub = (const unsigned char *)uv;
    size_t k = 0;
    if (s->nuv + cnt > s->nv) return -1;
    for (k = 0; k < cnt; k++) {
        double u = 0.0, v = 0.0;
        uint64_t vi = s->nuv + k;
        MeshChartRec *q = NULL;
        while (vi >= s->rec[s->r].v1) s->r++;
        q = &s->rec[s->r];
        if (f64) { memcpy(&u, ub + k * es, sizeof u); memcpy(&v, ub + k * es + sizeof u, sizeof v); }
        else {
            float fu = 0.0f, fv = 0.0f;
            memcpy(&fu, ub + k * es, sizeof fu); memcpy(&fv, ub + k * es + sizeof fu, sizeof fv);
            u = (double)fu; v = (double)fv;
        }
        if (u < q->uvlo[0]) q->uvlo[0] = u;
        if (u > q->uvhi[0]) q->uvhi[0] = u;
        if (v < q->uvlo[1]) q->uvlo[1] = v;
        if (v > q->uvhi[1]) q->uvhi[1] = v;
    }
    s->nuv += cnt;
    return 0;
}

/* faces -> each run's face range; a face spanning runs, or a run's faces in two places, is refused (1) */
static int mci_build_faces(MciBuild *s, const int32_t *fb, size_t cnt)
{
    size_t k = 0;
    for (k = 0; k < cnt; k++) {
        int32_t a = fb[k * 3], b = fb[k * 3 + 1], c = fb[k * 3 + 2];
        uint64_t f = s->nf + k;
        size_t q = 0;
        if (a < 0 || b < 0 || c < 0) return -1;
        q = mci_find(s->rec, s->n, (uint64_t)a, s->cur);
        if (q >= s->n || (uint64_t)b < s->rec[q].v0 || (uint64_t)b >= s->rec[q].v1 ||
            (uint64_t)c < s->rec[q].v0 || (uint64_t)c >= s->rec[q].v1) return 1;
        if (s->rec[q].f0 == MCI_NONE) { s->rec[q].f0 = f; s->rec[q].f1 = f + 1; }
        else if (s->rec[q].f1 == f) s->rec[q].f1 = f + 1;
        else return 1;
        s->cur = q;
    }
    s->nf += cnt;
    return 0;
}

static int mci_build_write(MciBuild *s, const char *out_idx, const MeshBinHeader *h, uint64_t chart_bytes, int64_t mtime)
{
    size_t r = 0;
    if (s->nv != h->nv || s->nf != h->nf || ((h->flags & MESH_BIN_HAS_UV) && s->nuv != h->nv)) return -1;
    for (r = 0; r < s->n; r++) if (s->rec[r].f0 == MCI_NONE) s->rec[r].f0 = s->rec[r].f1 = 0;
    return mci_write_index(out_idx, s->rec, s->n, h, chart_bytes, mtime) == 0 ? 0 : -1;
}

int MeshChartIndex_stamp(const char *vmesh, const char *chart_i32, MeshChartStamp *out)
{
    if (out == NULL) return -1;
    memset(out, 0, sizeof *out);
    if (vmesh == NULL || chart_i32 == NULL || MeshBin_header(vmesh, &out->mesh) != 0 ||
        mci_file_info(vmesh, NULL, &out->mtime) != 0 || mci_file_info(chart_i32, &out->chart_bytes, NULL) != 0 ||
        out->chart_bytes != out->mesh.nv * sizeof(int32_t))
        return -1;
    return 0;
}

int MeshChartIndex_build(const char *vmesh, const char *chart_i32, const char *out_idx)
{
    MeshChartStamp st;
    MciBuild s;
    FILE *fm = NULL, *fc = NULL;
    float *vb = NULL;
    int32_t *cb = NULL, *fb = NULL;
    unsigned char *ub = NULL;
    uint64_t i = 0;
    int rc = -1, f64 = 0;
    memset(&s, 0, sizeof s);
    if (out_idx == NULL || MeshChartIndex_stamp(vmesh, chart_i32, &st) != 0) return -1;
    f64 = (st.mesh.flags & MESH_BIN_UV_F64) != 0;
    fm = fopen(vmesh, "rb");
    fc = fopen(chart_i32, "rb");
    vb = (float *)malloc(MCI_BLOCK * 3 * sizeof *vb);
    cb = (int32_t *)malloc(MCI_BLOCK * sizeof *cb);
    fb = (int32_t *)malloc(MCI_BLOCK * 3 * sizeof *fb);
    ub = (unsigned char *)malloc(MCI_BLOCK * 2 * sizeof(double));
    if (fm == NULL || fc == NULL || vb == NULL || cb == NULL || fb == NULL || ub == NULL) goto out;
    /* pass 1: vertices and their charts -> runs and 3-D boxes */
    if (mci_seek(fm, st.mesh.vertices_offset) != 0) goto out;
    for (i = 0; i < st.mesh.nv;) {
        size_t cnt = (size_t)((st.mesh.nv - i) < MCI_BLOCK ? (st.mesh.nv - i) : MCI_BLOCK);
        if (mci_read(fm, vb, cnt * 3 * sizeof *vb) != 0 || mci_read(fc, cb, cnt * sizeof *cb) != 0) goto out;
        if ((rc = mci_build_verts(&s, vb, cb, cnt)) != 0) goto out;
        rc = -1;
        i += cnt;
    }
    if (mci_build_verts_end(&s) != 0) goto out;
    /* pass 2: UV boxes */
    if (st.mesh.flags & MESH_BIN_HAS_UV) {
        size_t es = f64 ? 2 * sizeof(double) : 2 * sizeof(float);
        if (mci_seek(fm, st.mesh.uv_offset) != 0) goto out;
        for (i = 0; i < st.mesh.nv;) {
            size_t cnt = (size_t)((st.mesh.nv - i) < MCI_BLOCK ? (st.mesh.nv - i) : MCI_BLOCK);
            if (mci_read(fm, ub, cnt * es) != 0 || mci_build_uv(&s, ub, f64, cnt) != 0) goto out;
            i += cnt;
        }
    }
    /* pass 3: faces */
    if (mci_seek(fm, st.mesh.faces_offset) != 0) goto out;
    for (i = 0; i < st.mesh.nf;) {
        size_t cnt = (size_t)((st.mesh.nf - i) < MCI_BLOCK ? (st.mesh.nf - i) : MCI_BLOCK);
        if (mci_read(fm, fb, cnt * 3 * sizeof *fb) != 0) goto out;
        if ((rc = mci_build_faces(&s, fb, cnt)) != 0) goto out;
        rc = -1;
        i += cnt;
    }
    rc = mci_build_write(&s, out_idx, &st.mesh, st.chart_bytes, st.mtime);
out:
    if (fm) fclose(fm);
    if (fc) fclose(fc);
    mci_build_free(&s);
    free(vb); free(cb); free(fb); free(ub);
    return rc;
}

int MeshChartIndex_build_mem(const MeshChartStamp *stamp, const char *vmesh, const char *chart_i32,
                             const MeshBinData *m, const int32_t *vchart, const char *out_idx)
{
    MeshChartStamp now;
    MciBuild s;
    int rc = -1;
    memset(&s, 0, sizeof s);
    if (stamp == NULL || m == NULL || vchart == NULL || out_idx == NULL) return -1;
    /* the files must still be the ones that were read: the same header numbers, size and mtime */
    if (MeshChartIndex_stamp(vmesh, chart_i32, &now) != 0 || memcmp(&now.mesh, &stamp->mesh, sizeof now.mesh) != 0 ||
        now.mtime != stamp->mtime || now.chart_bytes != stamp->chart_bytes || m->nv != stamp->mesh.nv ||
        m->nf != stamp->mesh.nf || ((stamp->mesh.flags & MESH_BIN_HAS_UV) && m->uv == NULL && m->uv64 == NULL))
        return -1;
    if ((rc = mci_build_verts(&s, m->verts, vchart, m->nv)) != 0) goto out;
    rc = -1;
    if (mci_build_verts_end(&s) != 0) goto out;
    if ((stamp->mesh.flags & MESH_BIN_HAS_UV) &&
        mci_build_uv(&s, m->uv64 ? (const void *)m->uv64 : (const void *)m->uv, m->uv64 != NULL, m->nv) != 0) goto out;
    if ((rc = mci_build_faces(&s, m->faces, m->nf)) != 0) goto out;
    rc = mci_build_write(&s, out_idx, &stamp->mesh, stamp->chart_bytes, stamp->mtime);
out:
    mci_build_free(&s);
    return rc;
}

int MeshChartIndex_load(const char *idx, const char *vmesh, MeshChartIndex *out)
{
    unsigned char magic[8];
    uint32_t version = 0, flags = 0;
    uint64_t nv = 0, nf = 0, bytes = 0, chart_bytes = 0, nrec = 0, k = 0;
    int64_t mtime = 0, now_mtime = 0;
    MeshBinHeader h;
    FILE *f = NULL;
    int ok = 0;
    if (out == NULL) return -1;
    memset(out, 0, sizeof *out);
    if (idx == NULL || vmesh == NULL || MeshBin_header(vmesh, &h) != 0 || mci_file_info(vmesh, NULL, &now_mtime) != 0)
        return -1;
    f = fopen(idx, "rb");
    if (f == NULL) return -1;
    ok = mci_read(f, magic, sizeof magic) == 0 && memcmp(magic, mci_magic, sizeof magic) == 0 &&
         mci_read(f, &version, sizeof version) == 0 && version == MCI_VERSION &&
         mci_read(f, &flags, sizeof flags) == 0 && mci_read(f, &nv, sizeof nv) == 0 &&
         mci_read(f, &nf, sizeof nf) == 0 && mci_read(f, &bytes, sizeof bytes) == 0 &&
         mci_read(f, &chart_bytes, sizeof chart_bytes) == 0 && mci_read(f, &mtime, sizeof mtime) == 0 &&
         mci_read(f, &nrec, sizeof nrec) == 0;
    /* stale: the mesh is not the one indexed */
    ok = ok && flags == h.flags && nv == h.nv && nf == h.nf && bytes == h.file_bytes && mtime == now_mtime &&
         nrec > 0 && nrec <= nv;
    if (ok) out->rec = (MeshChartRec *)calloc((size_t)nrec, sizeof *out->rec);
    ok = ok && out->rec != NULL;
    for (k = 0; ok && k < nrec; k++) {
        MeshChartRec *r = &out->rec[k];
        int32_t pad = 0;
        ok = mci_read(f, &r->chart, sizeof r->chart) == 0 && mci_read(f, &pad, sizeof pad) == 0 &&
             mci_read(f, &r->v0, sizeof r->v0) == 0 && mci_read(f, &r->v1, sizeof r->v1) == 0 &&
             mci_read(f, &r->f0, sizeof r->f0) == 0 && mci_read(f, &r->f1, sizeof r->f1) == 0 &&
             mci_read(f, r->lo, sizeof r->lo) == 0 && mci_read(f, r->hi, sizeof r->hi) == 0 &&
             mci_read(f, r->uvlo, sizeof r->uvlo) == 0 && mci_read(f, r->uvhi, sizeof r->uvhi) == 0 &&
             r->v0 < r->v1 && r->v1 <= nv && r->f0 <= r->f1 && r->f1 <= nf && (k == 0 || r->v0 == out->rec[k - 1].v1);
    }
    ok = ok && out->rec[nrec - 1].v1 == nv;
    fclose(f);
    if (!ok) { MeshChartIndex_free(out); return -1; }
    out->n = (size_t)nrec;
    out->mesh = h;
    out->chart_bytes = chart_bytes;
    return 0;
}

void MeshChartIndex_free(MeshChartIndex *ix)
{
    if (ix == NULL) return;
    free(ix->rec);
    memset(ix, 0, sizeof *ix);
}

/* Read one file section at a time. Adjacent selected ranges share a read; never read omitted charts. The old
 * chart-at-a-time reader sought from vertices to UVs to faces and back for every chart (minutes on a cold disk).
 * Face runs need not follow vertex order: only truly adjacent offsets are joined, preserving the compact order. */
/* Selected runs this close in the file are read through rather than sought past: on a mechanical disk one long read
 * beats a seek per run (a close-up pool's ~640,000 runs were ~70,000 seeks, minutes cold; 2026-10-30).  The gaps'
 * bytes go to scratch and are dropped, so the output holds exactly the selected runs. */
#define MCI_GAP_BYTES ((uint64_t)1 << 20)
#define MCI_SPAN_BYTES ((uint64_t)64 << 20)

static int mci_read_ranges(FILE *f, const MeshChartIndex *ix, const size_t *sel, size_t nsel,
                           int faces, uint64_t offset, size_t es, void *out)
{
    size_t s = 0, at = 0, scap = 0;
    unsigned char *scratch = NULL;
    int rc = 0;
    while (s < nsel && rc == 0) {
        const MeshChartRec *r = &ix->rec[sel[s]];
        uint64_t lo = faces ? r->f0 : r->v0, hi = faces ? r->f1 : r->v1;
        size_t end = s + 1, bytes = 0, k = 0;
        int gaps = 0;
        if (lo == hi) { s++; continue; }
        while (end < nsel) {
            const MeshChartRec *next = &ix->rec[sel[end]];
            uint64_t nlo = faces ? next->f0 : next->v0, nhi = faces ? next->f1 : next->v1;
            if (nlo == nhi) { end++; continue; }
            if (nlo < hi || (nlo - hi) * es > MCI_GAP_BYTES || (nhi - lo) * es > MCI_SPAN_BYTES) break;
            gaps |= nlo != hi;
            hi = nhi;
            end++;
        }
        bytes = (size_t)(hi - lo) * es;
        if (!gaps) {                   /* one contiguous span: straight into the output */
            if (mci_seek(f, offset + lo * es) != 0 || mci_read(f, (unsigned char *)out + at, bytes) != 0) rc = -1;
            at += bytes;
            s = end;
            continue;
        }
        if (bytes > scap) {
            unsigned char *grown = (unsigned char *)realloc(scratch, bytes);
            if (grown == NULL) { rc = -1; break; }
            scratch = grown; scap = bytes;
        }
        if (mci_seek(f, offset + lo * es) != 0 || mci_read(f, scratch, bytes) != 0) { rc = -1; break; }
        for (k = s; k < end; k++) {    /* the selected runs out of the span, in order */
            const MeshChartRec *q = &ix->rec[sel[k]];
            uint64_t qlo = faces ? q->f0 : q->v0, qhi = faces ? q->f1 : q->v1;
            size_t n = (size_t)(qhi - qlo) * es;
            if (n == 0) continue;
            memcpy((unsigned char *)out + at, scratch + (size_t)(qlo - lo) * es, n);
            at += n;
        }
        s = end;
    }
    free(scratch);
    return rc;
}

int MeshChartIndex_read(const MeshChartIndex *ix, const char *vmesh, const size_t *sel, size_t nsel,
                        MeshBinData *out, int32_t **vchart)
{
    FILE *f = NULL;
    uint64_t tv = 0, tf = 0;
    size_t s = 0, bv = 0, bf = 0;
    int has_uv = 0, f64 = 0, rc = -1;
    float *ufl = NULL;
    if (out == NULL || vchart == NULL) return -1;
    memset(out, 0, sizeof *out);
    *vchart = NULL;
    if (ix == NULL || ix->rec == NULL || vmesh == NULL || (nsel > 0 && sel == NULL)) return -1;
    has_uv = (ix->mesh.flags & MESH_BIN_HAS_UV) != 0;
    f64 = (ix->mesh.flags & MESH_BIN_UV_F64) != 0;
    for (s = 0; s < nsel; s++) {
        if (sel[s] >= ix->n || (s > 0 && sel[s] <= sel[s - 1])) return -1;     /* ascending, in range */
        tv += ix->rec[sel[s]].v1 - ix->rec[sel[s]].v0;
        tf += ix->rec[sel[s]].f1 - ix->rec[sel[s]].f0;
    }
    if (tv > (uint64_t)INT32_MAX || tv > SIZE_MAX / (2 * sizeof(double)) ||
        tf > SIZE_MAX / (3 * sizeof(int32_t))) return -1;
    out->verts = (float *)malloc((size_t)(tv ? tv : 1) * 3 * sizeof *out->verts);
    out->faces = (int32_t *)malloc((size_t)(tf ? tf : 1) * 3 * sizeof *out->faces);
    *vchart = (int32_t *)malloc((size_t)(tv ? tv : 1) * sizeof **vchart);
    if (has_uv) out->uv64 = (double *)malloc((size_t)(tv ? tv : 1) * 2 * sizeof *out->uv64);
    f = fopen(vmesh, "rb");
    if (out->verts == NULL || out->faces == NULL || *vchart == NULL || (has_uv && out->uv64 == NULL) || f == NULL) goto out;
    if (mci_read_ranges(f, ix, sel, nsel, 0, ix->mesh.vertices_offset, 3 * sizeof(float), out->verts) != 0) goto out;
    if (has_uv && f64 &&
        mci_read_ranges(f, ix, sel, nsel, 0, ix->mesh.uv_offset, 2 * sizeof(double), out->uv64) != 0) goto out;
    if (has_uv && !f64) {
        /* Float UVs widen exactly as before, in file order, using one bounded conversion buffer. */
        ufl = (float *)malloc(MCI_BLOCK * 2 * sizeof *ufl);
        if (ufl == NULL) goto out;
        for (s = 0; s < nsel; s++) {
            const MeshChartRec *r = &ix->rec[sel[s]];
            size_t nv = (size_t)(r->v1 - r->v0), at = 0;
            if (mci_seek(f, ix->mesh.uv_offset + r->v0 * 2 * sizeof(float)) != 0) goto out;
            while (at < nv) {
                size_t count = nv - at < MCI_BLOCK ? nv - at : MCI_BLOCK, k = 0;
                if (mci_read(f, ufl, count * 2 * sizeof *ufl) != 0) goto out;
                for (k = 0; k < count * 2; k++) out->uv64[(bv + at) * 2 + k] = (double)ufl[k];
                at += count;
            }
            bv += nv;
        }
        bv = 0;
    }
    if (mci_read_ranges(f, ix, sel, nsel, 1, ix->mesh.faces_offset, 3 * sizeof(int32_t), out->faces) != 0) goto out;
    for (s = 0; s < nsel; s++) {
        const MeshChartRec *r = &ix->rec[sel[s]];
        size_t nv = (size_t)(r->v1 - r->v0), nf = (size_t)(r->f1 - r->f0), k = 0;
        if (nf > 0) {
            for (k = 0; k < nf * 3; k++) {
                int64_t w = (int64_t)out->faces[bf * 3 + k] - (int64_t)r->v0;
                if (w < 0 || (uint64_t)w >= nv) goto out;
                out->faces[bf * 3 + k] = (int32_t)(w + (int64_t)bv);
            }
        }
        for (k = 0; k < nv; k++) (*vchart)[bv + k] = r->chart;
        bv += nv;
        bf += nf;
    }
    out->nv = bv;
    out->nf = bf;
    rc = 0;
out:
    if (f) fclose(f);
    free(ufl);
    if (rc != 0) { MeshBin_dispose(out); free(*vchart); *vchart = NULL; }
    return rc;
}

int MeshChartIndex_read_side(const MeshChartIndex *ix, const char *path, size_t es, const size_t *sel, size_t nsel,
                             void **out)
{
    FILE *f = NULL;
    uint64_t tv = 0, bytes = 0;
    size_t s = 0, b = 0;
    unsigned char *buf = NULL;
    if (out == NULL) return -1;
    *out = NULL;
    if (ix == NULL || path == NULL || es == 0 || mci_file_info(path, &bytes, NULL) != 0 || bytes != ix->mesh.nv * es)
        return -1;
    for (s = 0; s < nsel; s++) {
        if (sel[s] >= ix->n) return -1;
        tv += ix->rec[sel[s]].v1 - ix->rec[sel[s]].v0;
    }
    buf = (unsigned char *)malloc((size_t)(tv ? tv : 1) * es);
    f = fopen(path, "rb");
    if (buf == NULL || f == NULL) goto fail;
    for (s = 0; s < nsel; s++) {
        const MeshChartRec *r = &ix->rec[sel[s]];
        size_t nb = (size_t)(r->v1 - r->v0) * es;
        if (mci_seek(f, r->v0 * es) != 0 || mci_read(f, buf + b, nb) != 0) goto fail;
        b += nb;
    }
    fclose(f);
    *out = buf;
    return 0;
fail:
    if (f) fclose(f);
    free(buf);
    return -1;
}

/* ---------------------------------------------------------------------------------------------------- selftest */
static int mci_write_chart_file(const char *path, const int32_t *ch, size_t n)
{
    FILE *f = fopen(path, "wb");
    int ok = 0;
    if (f == NULL) return -1;
    ok = fwrite(ch, sizeof *ch, n, f) == n;
    if (fclose(f) != 0) ok = 0;
    return ok ? 0 : -1;
}

/* Exact whole-reader parity for every subset, all UV encodings, zero-face runs, and face runs in a different
 * order from the vertices. These are the cases in which coalescing may accidentally include or reorder data. */
static int mci_read_selftest(const char *dir)
{
    float v[36], uv32[24];
    double uv64[24];
    int32_t ch[12], faces[12] = {6, 7, 8, 0, 1, 2, 0, 2, 1, 3, 4, 5};
    char path[700], cp[700], ip[750];
    int mode = 0, fails = 0;
    size_t k = 0;
    for (k = 0; k < 36; k++) v[k] = (float)k * 0.25f;
    for (k = 0; k < 24; k++) { uv32[k] = (float)k * 0.125f; uv64[k] = (double)k + 0.123456789; }
    for (k = 0; k < 12; k++) ch[k] = (int32_t)(k / 3) + 11;
    snprintf(cp, sizeof cp, "%s/ranges_chart.i32", dir);
    if (mci_write_chart_file(cp, ch, 12) != 0) return 1;
    for (mode = 0; mode < 3; mode++) {
        MeshChartIndex ix;
        MeshBinData whole;
        unsigned int mask = 0;
        memset(&ix, 0, sizeof ix); memset(&whole, 0, sizeof whole);
        snprintf(path, sizeof path, "%s/ranges%d.vmesh", dir, mode);
        if ((mode == 2 ? MeshBin_write_uv64(path, v, 12, faces, 4, uv64) :
                         MeshBin_write(path, v, 12, faces, 4, mode == 1 ? uv32 : NULL)) != 0 ||
            MeshChartIndex_path(path, ip, sizeof ip) != 0 || MeshChartIndex_build(path, cp, ip) != 0 ||
            MeshChartIndex_load(ip, path, &ix) != 0 || MeshBin_read_precise_malloc(path, &whole) != 0) {
            MeshChartIndex_free(&ix); MeshBin_dispose(&whole); return fails + 1;
        }
        for (mask = 0; mask < 16; mask++) {
            size_t sel[4], ns = 0, bv = 0, bf = 0, s = 0;
            MeshBinData got;
            int32_t *vc = NULL;
            for (k = 0; k < 4; k++) if (mask & (1u << k)) sel[ns++] = k;
            if (MeshChartIndex_read(&ix, path, sel, ns, &got, &vc) != 0) { fails++; continue; }
            for (s = 0; s < ns; s++) {
                const MeshChartRec *r = &ix.rec[sel[s]];
                size_t nv = (size_t)(r->v1 - r->v0), nf = (size_t)(r->f1 - r->f0);
                if (memcmp(got.verts + bv * 3, whole.verts + (size_t)r->v0 * 3, nv * 3 * sizeof(float)) != 0) fails++;
                if (mode == 2 && memcmp(got.uv64 + bv * 2, whole.uv64 + (size_t)r->v0 * 2, nv * 2 * sizeof(double)) != 0) fails++;
                if (mode == 1) for (k = 0; k < nv * 2; k++) {
                    double expected = (double)whole.uv[(size_t)r->v0 * 2 + k];
                    if (memcmp(got.uv64 + bv * 2 + k, &expected, sizeof expected) != 0) fails++;
                }
                if (!mode && got.uv64 != NULL) fails++;
                for (k = 0; k < nv; k++) if (vc[bv + k] != r->chart) fails++;
                for (k = 0; k < nf * 3; k++)
                    if (got.faces[bf * 3 + k] != whole.faces[(size_t)r->f0 * 3 + k] - (int32_t)r->v0 + (int32_t)bv) fails++;
                bv += nv; bf += nf;
            }
            if (got.nv != bv || got.nf != bf) fails++;
            MeshBin_dispose(&got); free(vc);
        }
        {
            MeshBinData got;
            int32_t *vc = NULL;
            size_t bad[][2] = {{1, 0}, {1, 1}, {0, 4}};
            for (k = 0; k < 3; k++) {
                if (MeshChartIndex_read(&ix, path, bad[k], 2, &got, &vc) == 0 || got.verts || vc) fails++;
                MeshBin_dispose(&got); free(vc);
            }
            /* A short read must free all partially filled output. */
            { FILE *f = fopen(path, "wb"); if (f) fclose(f); else fails++; }
            { size_t one = 0;
              if (MeshChartIndex_read(&ix, path, &one, 1, &got, &vc) == 0 || got.verts || vc) fails++;
              MeshBin_dispose(&got); free(vc); }
        }
        MeshBin_dispose(&whole); MeshChartIndex_free(&ix);
    }
    fprintf(stderr, "  mesh_chart_index section reads: 48 subsets + invalid/truncated reads: %s\n", fails ? "FAIL" : "ok");
    return fails;
}

int MeshChartIndex_selftest(void)
{
    /* three charts (ids 5, 2, -1) of 4, 3, 3 vertices; faces chart-contiguous and self-contained */
    float verts[30] = { 0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0,
                        10, 10, 10,  11, 10, 10,  10, 11, 10,
                        20, 0, 5,  21, 0, 5,  20, 1, 5 };
    double uv[20] = { 0, 0,  1, 0,  1, 1,  0, 1,   100, 100,  101, 100,  100, 101,   50, 0,  51, 0,  50, 1 };
    int32_t faces[12] = { 0, 1, 2,  0, 2, 3,   4, 5, 6,   7, 8, 9 };
    int32_t charts[10] = { 5, 5, 5, 5,  2, 2, 2,  -1, -1, -1 };
    int32_t bad[10] = { 5, 5, 2, 2,  5, 5, 5,  -1, -1, -1 };        /* chart 5 in two runs */
    char dir[512], mesh[600], cpath[600], idx[700], bpath[600];
    MeshChartIndex ix;
    MeshBinData m;
    int32_t *vc = NULL;
    size_t sel[2] = { 0, 2 };
    int fails = 0;
    memset(&ix, 0, sizeof ix);
    memset(&m, 0, sizeof m);
    if (ves_selftest_dir(dir, sizeof dir, "mesh_chart_index") != 0) return 1;
    fails += mci_read_selftest(dir);
    snprintf(mesh, sizeof mesh, "%s/m.vmesh", dir);
    snprintf(cpath, sizeof cpath, "%s/m_chart.i32", dir);
    snprintf(bpath, sizeof bpath, "%s/bad_chart.i32", dir);
    if (MeshChartIndex_path(mesh, idx, sizeof idx) != 0 ||
        MeshBin_write_uv64(mesh, verts, 10, faces, 4, uv) != 0 || mci_write_chart_file(cpath, charts, 10) != 0 ||
        mci_write_chart_file(bpath, bad, 10) != 0) {
        fprintf(stderr, "  mesh_chart_index: cannot write the fixture\n");
        return 1;
    }
    if (MeshChartIndex_build(mesh, cpath, idx) != 0) { fprintf(stderr, "  mesh_chart_index: build failed\n"); fails++; }
    else if (MeshChartIndex_load(idx, mesh, &ix) != 0 || ix.n != 3) {
        fprintf(stderr, "  mesh_chart_index: load failed\n"); fails++;
    } else {
        const MeshChartRec *a = &ix.rec[0], *b = &ix.rec[1], *c = &ix.rec[2];
        if (a->chart != 5 || a->v0 != 0 || a->v1 != 4 || a->f0 != 0 || a->f1 != 2 || b->chart != 2 || b->v0 != 4 ||
            b->f0 != 2 || b->f1 != 3 || c->chart != -1 || c->f0 != 3 || c->f1 != 4 ||
            fabs(a->hi[0] - 1.0) > 1e-12 || fabs(b->lo[0] - 10.0) > 1e-12 || fabs(c->hi[0] - 21.0) > 1e-12 ||
            fabs(b->uvlo[0] - 100.0) > 1e-12 || fabs(a->uvhi[1] - 1.0) > 1e-12) {
            fprintf(stderr, "  mesh_chart_index: records wrong\n"); fails++;
        }
        /* charts 0 and 2: a compact mesh, faces renumbered, file order kept */
        if (MeshChartIndex_read(&ix, mesh, sel, 2, &m, &vc) != 0 || m.nv != 7 || m.nf != 3 || m.uv64 == NULL ||
            m.faces[0] != 0 || m.faces[5] != 3 || m.faces[6] != 4 || m.faces[7] != 5 || m.faces[8] != 6 ||
            fabs((double)m.verts[4 * 3] - 20.0) > 1e-12 || fabs(m.uv64[4 * 2] - 50.0) > 1e-12 ||
            vc[0] != 5 || vc[3] != 5 || vc[4] != -1 || vc[6] != -1) {
            fprintf(stderr, "  mesh_chart_index: paged read wrong\n"); fails++;
        }
        MeshBin_dispose(&m);
        free(vc);
        vc = NULL;
        {
            void *side = NULL;
            uint8_t attr[10] = { 1, 1, 1, 1, 2, 2, 2, 3, 3, 3 };
            char ap[600];
            snprintf(ap, sizeof ap, "%s/m_attr.u8", dir);
            {
                FILE *fa = fopen(ap, "wb");
                if (fa != NULL) { fwrite(attr, 1, 10, fa); fclose(fa); }
            }
            if (MeshChartIndex_read_side(&ix, ap, 1, sel, 2, &side) != 0 || ((uint8_t *)side)[0] != 1 ||
                ((uint8_t *)side)[4] != 3 || ((uint8_t *)side)[6] != 3) {
                fprintf(stderr, "  mesh_chart_index: side file read wrong\n"); fails++;
            }
            free(side);
        }
        MeshChartIndex_free(&ix);
    }
    {   /* the index from a whole read is the streamed index, byte for byte; a changed stamp writes nothing */
        MeshChartStamp st;
        MeshBinData whole;
        char imem[700];
        unsigned char *a = NULL, *b = NULL;
        size_t na = 0, nb = 0;
        memset(&whole, 0, sizeof whole);
        snprintf(imem, sizeof imem, "%s/m_mem.chartidx", dir);
        if (MeshChartIndex_stamp(mesh, cpath, &st) != 0 || MeshBin_read_precise_malloc(mesh, &whole) != 0 ||
            MeshChartIndex_build_mem(&st, mesh, cpath, &whole, charts, imem) != 0) {
            fprintf(stderr, "  mesh_chart_index: build from memory failed\n"); fails++;
        } else {
            FILE *fa = fopen(idx, "rb"), *fb = fopen(imem, "rb");
            if (fa && fb) {
                a = (unsigned char *)malloc(1 << 16); b = (unsigned char *)malloc(1 << 16);
                if (a && b) { na = fread(a, 1, 1 << 16, fa); nb = fread(b, 1, 1 << 16, fb); }
            }
            if (fa) fclose(fa);
            if (fb) fclose(fb);
            if (a == NULL || b == NULL || na == 0 || na != nb || memcmp(a, b, na) != 0) {
                fprintf(stderr, "  mesh_chart_index: the in-memory index differs from the streamed one\n"); fails++;
            }
            free(a); free(b);
            st.mtime += 1;                          /* the mesh "changed" after it was read: refused */
            if (MeshChartIndex_build_mem(&st, mesh, cpath, &whole, charts, imem) != -1) {
                fprintf(stderr, "  mesh_chart_index: a changed mesh was indexed from memory\n"); fails++;
            }
        }
        MeshBin_dispose(&whole);
    }
    if (MeshChartIndex_build(mesh, bpath, idx) != 1) {
        fprintf(stderr, "  mesh_chart_index: a chart in two runs was indexed\n"); fails++;
    }
    {   /* a face spanning two runs is refused; a rewritten mesh makes its old index stale */
        int32_t span[12] = { 0, 1, 2,  0, 2, 4,   4, 5, 6,   7, 8, 9 };
        char m2[600], i2[700];
        snprintf(m2, sizeof m2, "%s/span.vmesh", dir);
        if (MeshChartIndex_path(m2, i2, sizeof i2) != 0 || MeshBin_write_uv64(m2, verts, 10, span, 4, uv) != 0 ||
            MeshChartIndex_build(m2, cpath, i2) != 1) {
            fprintf(stderr, "  mesh_chart_index: a face spanning two charts was indexed\n"); fails++;
        }
        if (MeshBin_write_uv64(mesh, verts, 10, faces, 3, uv) != 0 || MeshChartIndex_load(idx, mesh, &ix) == 0) {
            fprintf(stderr, "  mesh_chart_index: a stale index was accepted\n"); fails++;
            MeshChartIndex_free(&ix);
        }
    }
    fprintf(stderr, "  mesh_chart_index selftest: %s\n", fails ? "FAIL" : "ok");
    if (!fails) (void)ves_selftest_clean(dir);
    return fails ? 1 : 0;
}
