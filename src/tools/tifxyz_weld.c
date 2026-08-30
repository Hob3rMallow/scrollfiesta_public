/*
 * tifxyz_weld.c -- weld z-slab tifxyz strips into full-height per-wrap segments.
 *
 * The full-scroll pipeline processes the band as independent z-slabs (1x21x21
 * strips), each exported to its own tifxyz atlas from ONE shared global solve
 * (scroll_unroll --z-range). Because the strips share the u lattice and cover
 * disjoint v (world z) bands, welding is a vertical stack: for each wrap k this
 * tool places every strip's k-segment at its absolute (u,v) into one
 * full-height canvas. Halo overlaps are coordinate-checked and averaged when
 * consistent; conflicts are hidden and fail the verified assembly budget. The
 * C form of python/scripts/weld_slab_atlas.py -- reassembly is the
 * scaling mechanism for the 241k-cube whole scroll.
 *
 * usage: tifxyz_weld [--coord-tol T] [--max-contested-fraction F]
 *                    [--lattice-tol L] <out_root> <strip_root...>
 *        tifxyz_weld --selftest
 *
 * out_root gets seg/welded_w###/{x,y,z,mask,provenance}.tif + meta.json and
 * atlas.json. Segments carry absolute origin_uv; a wrap's members are all the
 * "*_w<k>_*" segments across the input roots.
 */
#include "../common/ves_platform.h"

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef _WIN32
  #include <windows.h>
#else
  #include <dirent.h>
  #include <sys/stat.h>
#endif

#include "../common/arena.h"
#include "../common/tiff_io.h"

static int regular_file_exists(const char *path)
{
#ifdef _WIN32
    DWORD attributes = GetFileAttributesA(path);
    return attributes != INVALID_FILE_ATTRIBUTES &&
           !(attributes & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat status;
    return stat(path, &status) == 0 && S_ISREG(status.st_mode);
#endif
}

/* crude json key scan (same convention as the placed_index / tifxyz_render
 * readers): find "key" then the ':' after it. */
static const char *jfind(const char *s, const char *key)
{
    char pat[128];
    snprintf(pat, sizeof(pat), "\"%s\"", key);
    const char *p = strstr(s, pat);
    if (p == NULL) return NULL;
    p = strchr(p + strlen(pat), ':');
    return p != NULL ? p + 1 : NULL;
}

typedef struct {
    char   dir[1024];        /* full segment dir (…/seg/<name>) */
    int    k;                /* wrap index */
    double ou, ov, du, dv;   /* absolute origin (vox) + grid step */
    int    W, H;
} SegRef;

typedef struct {
    double coord_tol, max_contested_fraction, lattice_tol;
} WeldParams;

/* Wrap k from "…_w[+-]?<digits>_…". INT_MIN means no parse. */
static int parse_wrap_k(const char *name)
{
    const char *p = name;
    while ((p = strstr(p, "_w")) != NULL) {
        const char *number = p + 2;
        char *end = NULL;
        errno = 0;
        long value = strtol(number, &end, 10);
        if (errno == 0 && end != number && (*end == '_' || *end == '\0') &&
            value > INT_MIN && value <= INT_MAX)
            return (int)value;
        p += 2;
    }
    return INT_MIN;
}

static void wrap_name(char *out, size_t capacity, int k)
{
    if (k < 0) snprintf(out, capacity, "-%03d", -k);
    else snprintf(out, capacity, "%03d", k);
}

/* read origin_uv/du/dv/flip from meta.json + W/H from mask.tif. 0 on success. */
static int read_ref(Arena_T scr, const char *seg_dir, SegRef *r)
{
    char path[1200];
    snprintf(path, sizeof(path), "%s/meta.json", seg_dir);
    FILE *f = fopen(path, "rb");
    if (f == NULL) return -1;
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    buf[n] = '\0';
    const char *p;
    r->du = r->dv = r->ou = r->ov = 0.0;
    int flip_u = 0, flip_v = 0;
    p = jfind(buf, "du");
    if (p == NULL || sscanf(p, " %lf", &r->du) != 1) goto invalid_meta;
    p = jfind(buf, "dv");
    if (p == NULL || sscanf(p, " %lf", &r->dv) != 1) goto invalid_meta;
    p = jfind(buf, "origin_uv");
    if (p == NULL || sscanf(p, " [ %lf , %lf", &r->ou, &r->ov) != 2)
        goto invalid_meta;
    if ((p = jfind(buf, "flip_u")) != NULL && sscanf(p, " %d", &flip_u) != 1)
        goto invalid_meta;
    if ((p = jfind(buf, "flip_v")) != NULL && sscanf(p, " %d", &flip_v) != 1)
        goto invalid_meta;
    if (flip_u || flip_v) {
        fprintf(stderr, "tifxyz_weld: flip_u/flip_v unsupported (%s)\n", seg_dir);
        return -1;
    }
    if (!isfinite(r->du) || !isfinite(r->dv) ||
        !isfinite(r->ou) || !isfinite(r->ov) || r->du <= 0.0 || r->dv <= 0.0)
        goto invalid_meta;
    snprintf(path, sizeof(path), "%s/mask.tif", seg_dir);
    uint8_t *m = NULL; int D = 0, H = 0, W = 0;
    if (TiffIO_load(scr, path, &m, &D, &H, &W) != 0 || D != 1 || H < 1 || W < 1)
        return -1;
    r->W = W; r->H = H;
    snprintf(r->dir, sizeof(r->dir), "%s", seg_dir);
    return 0;
invalid_meta:
    fprintf(stderr, "tifxyz_weld: invalid lattice metadata (%s)\n", seg_dir);
    return -1;
}

static int append_ref(SegRef **refs, size_t *n, size_t *cap, const SegRef *r)
{
    if (*n == *cap) {
        size_t next_cap = *cap ? *cap * 2 : 256;
        if (next_cap < *cap || next_cap > SIZE_MAX / sizeof(SegRef)) return -1;
        SegRef *next = (SegRef *)realloc(*refs, next_cap * sizeof(SegRef));
        if (next == NULL) return -1;
        *refs = next;
        *cap = next_cap;
    }
    (*refs)[(*n)++] = *r;
    return 0;
}

/* Enumerate <root>/seg/<name>/ dirs. A recognized but malformed segment makes
 * the whole input fail; verified assembly must never silently omit a slab. */
static int scan_segs(Arena_T scr, const char *root,
                     SegRef **refs, size_t *n, size_t *cap)
{
    char segroot[1100];
    snprintf(segroot, sizeof(segroot), "%s/seg", root);
#ifdef _WIN32
    char glob[1200]; snprintf(glob, sizeof(glob), "%s\\*", segroot);
    WIN32_FIND_DATAA fd; HANDLE h = FindFirstFileA(glob, &fd);
    if (h == INVALID_HANDLE_VALUE) {
        fprintf(stderr, "tifxyz_weld: cannot enumerate %s\n", segroot);
        return -1;
    }
    int rc = 0;
    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) continue;
        if (fd.cFileName[0] == '.') continue;
        int k = parse_wrap_k(fd.cFileName);
        if (k == INT_MIN) continue;
        char sd[1200]; snprintf(sd, sizeof(sd), "%s/%s", segroot, fd.cFileName);
        Arena_free(scr);
        SegRef r; r.k = k;
        if (read_ref(scr, sd, &r) != 0) {
            fprintf(stderr, "tifxyz_weld: invalid segment %s\n", sd);
            rc = -1; break;
        }
        if (append_ref(refs, n, cap, &r) != 0) { rc = -1; break; }
    } while (FindNextFileA(h, &fd));
    FindClose(h);
    return rc;
#else
    DIR *d = opendir(segroot);
    if (d == NULL) {
        fprintf(stderr, "tifxyz_weld: cannot enumerate %s\n", segroot);
        return -1;
    }
    int rc = 0;
    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
        if (e->d_name[0] == '.') continue;
        int k = parse_wrap_k(e->d_name);
        if (k == INT_MIN) continue;
        char sd[1200]; snprintf(sd, sizeof(sd), "%s/%s", segroot, e->d_name);
        Arena_free(scr);
        SegRef r; r.k = k;
        if (read_ref(scr, sd, &r) != 0) {
            fprintf(stderr, "tifxyz_weld: invalid segment %s\n", sd);
            rc = -1; break;
        }
        if (append_ref(refs, n, cap, &r) != 0) { rc = -1; break; }
    }
    closedir(d);
    return rc;
#endif
}

static int segref_cmp(const void *va, const void *vb)
{
    const SegRef *a = (const SegRef *)va, *b = (const SegRef *)vb;
    if (a->k != b->k) return a->k < b->k ? -1 : 1;
    return strcmp(a->dir, b->dir);
}

typedef struct {
    int k, W, H;
    size_t nslab;
    double ou, ov, du, dv;
    uint64_t valid, overlap, contested;
} Piece;

/* weld all wraps found under the strip roots into out_root. 0 on success. */
static int weld_run(const char *out_root, char **roots, int nroot,
                    const WeldParams *params,
                    size_t *out_nwrap, uint64_t *out_valid)
{
    Arena_T scr = Arena_new();
    SegRef *refs = NULL; size_t nseg = 0, cap = 0;
    for (int i = 0; i < nroot; i++) {
        if (scan_segs(scr, roots[i], &refs, &nseg, &cap) != 0) {
            fprintf(stderr, "tifxyz_weld: source scan failed (%s)\n", roots[i]);
            free(refs); Arena_free(scr); return 1;
        }
    }
    if (nseg == 0) {
        fprintf(stderr, "tifxyz_weld: no segments under the %d root(s)\n", nroot);
        free(refs); Arena_free(scr); return 1;
    }
    qsort(refs, nseg, sizeof(*refs), segref_cmp);
    int kmin = INT_MAX, kmax = INT_MIN;
    for (size_t i = 0; i < nseg; i++) {
        if (refs[i].k < kmin) kmin = refs[i].k;
        if (refs[i].k > kmax) kmax = refs[i].k;
    }
    fprintf(stderr, "tifxyz_weld: %zu strip segments, wraps k[%d,%d]\n",
            nseg, kmin, kmax);

    Piece *pieces = (Piece *)malloc(nseg * sizeof(Piece));
    if (pieces == NULL) {
        fprintf(stderr, "tifxyz_weld: out of memory allocating wrap table\n");
        free(refs); Arena_free(scr); return 1;
    }
    size_t npiece = 0;
    uint64_t total_valid = 0, total_contested = 0;
    int failed = 0;

    size_t group_end = 0;
    for (size_t group = 0; group < nseg; group = group_end) {
        int k = refs[group].k;
        group_end = group + 1;
        while (group_end < nseg && refs[group_end].k == k) group_end++;
        /* members = strip segments for this wrap */
        double u0 = 1e300, u1 = -1e300, v0 = 1e300, v1 = -1e300;
        double du = 0.0, dv = 0.0;
        size_t nmem = 0;
        int lattice_bad = 0, input_bad = 0;
        for (size_t i = group; i < group_end; i++) {
            SegRef *s = &refs[i];
            if (nmem == 0) { du = s->du; dv = s->dv; }
            else if (fabs(s->du-du) > params->lattice_tol*fmax(1.0, fabs(du)) ||
                     fabs(s->dv-dv) > params->lattice_tol*fmax(1.0, fabs(dv))) {
                fprintf(stderr, "  w%03d: incompatible lattice step in %s\n", k, s->dir);
                lattice_bad = 1;
            }
            if (s->ou < u0) u0 = s->ou;
            if (s->ou + s->W * s->du > u1) u1 = s->ou + s->W * s->du;
            if (s->ov < v0) v0 = s->ov;
            if (s->ov + s->H * s->dv > v1) v1 = s->ov + s->H * s->dv;
            nmem++;
        }
        if (!isfinite(u0) || !isfinite(u1) || !isfinite(v0) || !isfinite(v1) ||
            group_end - group > INT_MAX) {
            fprintf(stderr, "  w%03d: non-finite extent or too many members\n", k);
            failed = 1; continue;
        }
        for (size_t i = group; i < group_end && !lattice_bad; i++) {
            double col = (refs[i].ou-u0)/du, row = (refs[i].ov-v0)/dv;
            if (!isfinite(col) || !isfinite(row) ||
                fabs(col-floor(col+0.5)) > params->lattice_tol ||
                fabs(row-floor(row+0.5)) > params->lattice_tol) {
                fprintf(stderr, "  w%03d: off-lattice origin in %s\n", k, refs[i].dir);
                lattice_bad = 1;
            }
        }
        if (lattice_bad) { failed = 1; continue; }
        double Wd = floor((u1 - u0) / du + 0.5);
        double Hd = floor((v1 - v0) / dv + 0.5);
        if (!isfinite(Wd) || !isfinite(Hd) || Wd < 1.0 || Hd < 1.0 ||
            Wd > INT_MAX || Hd > INT_MAX || Wd * Hd > 950000000.0) {
            fprintf(stderr, "  w%03d: invalid canvas (%.0fx%.0f px)\n", k, Wd, Hd);
            failed = 1;
            continue;
        }
        int Wk = (int)Wd, Hk = (int)Hd;
        size_t np = (size_t)Wk * (size_t)Hk;
        float *X = (float *)malloc(np * sizeof(float));
        float *Y = (float *)malloc(np * sizeof(float));
        float *Z = (float *)malloc(np * sizeof(float));
        uint8_t *M = (uint8_t *)calloc(np, 1);
        uint8_t *P = (uint8_t *)calloc(np, 1);
        uint32_t *cover = (uint32_t *)calloc(np, sizeof(uint32_t));
        if (!X || !Y || !Z || !M || !P || !cover) {
            free(X);free(Y);free(Z);free(M);free(P);free(cover);
            free(pieces); free(refs); Arena_free(scr); return 1; }
        for (size_t i = 0; i < np; i++) { X[i] = -1.0f; Y[i] = -1.0f; Z[i] = -1.0f; }

        for (size_t i = group; i < group_end; i++) {
            SegRef *s = &refs[i];
            Arena_free(scr);
            char p[1300]; float *sx=NULL,*sy=NULL,*sz=NULL; int sw=0,sh=0,w2,h2;
            snprintf(p, sizeof p, "%s/x.tif", s->dir);
            if (TiffIO_load_float2d(scr, p, &sx, &sw, &sh) != 0) {
                fprintf(stderr, "  w%03d: unreadable x.tif in %s\n", k, s->dir);
                input_bad = 1; continue;
            }
            snprintf(p, sizeof p, "%s/y.tif", s->dir);
            if (TiffIO_load_float2d(scr, p, &sy, &w2, &h2) != 0 || w2!=sw || h2!=sh) {
                fprintf(stderr, "  w%03d: invalid y.tif in %s\n", k, s->dir);
                input_bad = 1; continue;
            }
            snprintf(p, sizeof p, "%s/z.tif", s->dir);
            if (TiffIO_load_float2d(scr, p, &sz, &w2, &h2) != 0 || w2!=sw || h2!=sh) {
                fprintf(stderr, "  w%03d: invalid z.tif in %s\n", k, s->dir);
                input_bad = 1; continue;
            }
            uint8_t *sm=NULL,*sp=NULL; int D3=0,hh=0,ww=0;
            snprintf(p, sizeof p, "%s/mask.tif", s->dir);
            if (TiffIO_load(scr, p, &sm, &D3, &hh, &ww) != 0 ||
                D3!=1 || ww!=sw || hh!=sh) {
                fprintf(stderr, "  w%03d: invalid mask.tif in %s\n", k, s->dir);
                input_bad = 1; continue;
            }
            snprintf(p, sizeof p, "%s/provenance.tif", s->dir);
            if (regular_file_exists(p)) {
                if (TiffIO_load(scr, p, &sp, &D3, &hh, &ww) != 0 ||
                    D3 != 1 || ww != sw || hh != sh) {
                    fprintf(stderr, "  w%03d: invalid provenance.tif in %s\n", k, s->dir);
                    input_bad = 1; continue;
                }
            }
            int col0 = (int)floor((s->ou - u0) / du + 0.5);
            int row0 = (int)floor((s->ov - v0) / dv + 0.5);
            for (int r = 0; r < sh; r++) {
                int dr = r + row0;
                if (dr < 0 || dr >= Hk) continue;
                for (int c = 0; c < sw; c++) {
                    size_t si = (size_t)r * (size_t)sw + (size_t)c;
                    int dc = c + col0;
                    if (dc < 0 || dc >= Wk) continue;
                    size_t di = (size_t)dr * (size_t)Wk + (size_t)dc;
                    if (sp != NULL) {
                        if (sp[si] > 3 ||
                            (sp[si] == 0 && sm[si] == 255) ||
                            (sp[si] == 3 && sm[si] == 255) ||
                            ((sp[si] == 1 || sp[si] == 2) && sm[si] < 255)) {
                            fprintf(stderr, "  w%03d: inconsistent provenance/mask in %s\n",
                                    k, s->dir);
                            input_bad = 1; continue;
                        }
                        if (sp[si] == 3) {
                            X[di]=Y[di]=Z[di]=-1.0f; M[di]=0; P[di]=3; cover[di]=0;
                            continue;
                        }
                    }
                    if (sm[si] < 255) continue;
                    if (!isfinite(sx[si]) || !isfinite(sy[si]) ||
                        !isfinite(sz[si]) || sx[si] < 0.0f ||
                        sy[si] < 0.0f || sz[si] < 0.0f) {
                        fprintf(stderr, "  w%03d: mask exposes invalid xyz in %s\n", k, s->dir);
                        input_bad = 1; continue;
                    }
                    if (P[di] == 3) continue; /* a conflict is never resurrected */
                    if (M[di] == 0) {
                        X[di] = sx[si]; Y[di] = sy[si]; Z[di] = sz[si];
                        M[di] = 255; P[di] = 1; cover[di] = 1;
                    } else {
                        double dx = (double)X[di]-sx[si];
                        double dy = (double)Y[di]-sy[si];
                        double dz = (double)Z[di]-sz[si];
                        double difference = sqrt(dx*dx+dy*dy+dz*dz);
                        if (difference <= params->coord_tol) {
                            double count = (double)cover[di];
                            X[di] = (float)((count*X[di]+sx[si])/(count+1.0));
                            Y[di] = (float)((count*Y[di]+sy[si])/(count+1.0));
                            Z[di] = (float)((count*Z[di]+sz[si])/(count+1.0));
                            if (cover[di] < UINT32_MAX) cover[di]++;
                            P[di] = 2;
                        } else {
                            X[di]=Y[di]=Z[di]=-1.0f; M[di]=0; P[di]=3; cover[di]=0;
                        }
                    }
                }
            }
        }

        if (input_bad) {
            fprintf(stderr, "  w%03d: input validation failed\n", k);
            free(X);free(Y);free(Z);free(M);free(P);free(cover);
            failed = 1; continue;
        }

        uint64_t nv = 0, no = 0, nc = 0;
        double blo[3] = { 1e300,1e300,1e300 }, bhi[3] = { -1e300,-1e300,-1e300 };
        for (size_t i = 0; i < np; i++) {
            if (P[i] == 3) { nc++; continue; }
            if (M[i] < 255) continue;
            nv++;
            if (P[i] == 2) no++;
            if (X[i] < blo[0]) blo[0] = X[i]; if (X[i] > bhi[0]) bhi[0] = X[i];
            if (Y[i] < blo[1]) blo[1] = Y[i]; if (Y[i] > bhi[1]) bhi[1] = Y[i];
            if (Z[i] < blo[2]) blo[2] = Z[i]; if (Z[i] > bhi[2]) bhi[2] = Z[i];
        }
        if (nv == 0) {
            total_contested += nc;
            fprintf(stderr, "  w%03d: no verified pixels (%" PRIu64 " contested)\n", k, nc);
            free(X);free(Y);free(Z);free(M);free(P);free(cover);
            failed = 1; continue;
        }
        double contested_fraction = (double)nc / (double)(nv + nc);
        if (contested_fraction > params->max_contested_fraction + 1e-15) {
            fprintf(stderr,
                    "  w%03d: contested fraction %.8f exceeds %.8f (%" PRIu64 " cells)\n",
                    k, contested_fraction, params->max_contested_fraction, nc);
            total_contested += nc;
            free(X);free(Y);free(Z);free(M);free(P);free(cover);
            failed = 1; continue;
        }

        char od[1200], p[1300], wkey[32];
        wrap_name(wkey, sizeof(wkey), k);
        snprintf(od, sizeof od, "%s/seg/welded_w%s", out_root, wkey);
        snprintf(p, sizeof p, "%s/x.tif", od);
        int wr = 0;
        if (ves_ensure_parent_dir(p) != 0) {
            fprintf(stderr, "mkdir %s failed\n", od); wr = 1;
        }
        wr |= TiffIO_save_float2d(p, X, Wk, Hk);
        snprintf(p, sizeof p, "%s/y.tif", od); wr |= TiffIO_save_float2d(p, Y, Wk, Hk);
        snprintf(p, sizeof p, "%s/z.tif", od); wr |= TiffIO_save_float2d(p, Z, Wk, Hk);
        snprintf(p, sizeof p, "%s/mask.tif", od); wr |= TiffIO_save(p, M, 1, Hk, Wk);
        snprintf(p, sizeof p, "%s/provenance.tif", od); wr |= TiffIO_save(p, P, 1, Hk, Wk);
        snprintf(p, sizeof p, "%s/meta.json", od);
        FILE *mf = fopen(p, "w");
        if (mf) {
            fprintf(mf,
                "{\"scale\": [%.17g, %.17g], \"type\": \"seg\", \"format\": \"tifxyz\", "
                "\"uuid\": \"welded_w%s\", \"source\": \"scrollfiesta-welded\", "
                "\"scrollfiesta\": {\"du\": %.17g, \"dv\": %.17g, \"flip_u\": 0, \"flip_v\": 0, "
                "\"origin_uv\": [%.17g, %.17g], \"wrap_k\": %d, \"n_slabs\": %zu, "
                "\"valid_px\": %" PRIu64 ", \"overlap_px\": %" PRIu64
                ", \"contested_px\": %" PRIu64 ", "
                "\"coordinate_tolerance\": %.17g, \"verified\": true, "
                "\"provenance\": \"provenance.tif: 0=empty 1=single 2=overlap 3=contested\"}, "
                "\"bbox\": [[%.9g, %.9g, %.9g], [%.9g, %.9g, %.9g]]}\n",
                1.0/du, 1.0/dv, wkey, du, dv, u0, v0, k, nmem, nv, no, nc,
                params->coord_tol,
                blo[0], blo[1], blo[2], bhi[0], bhi[1], bhi[2]);
            if (fclose(mf) != 0) wr = 1;
        } else wr = 1;
        if (wr) { fprintf(stderr, "  w%03d: a band write failed\n", k); failed = 1; }

        pieces[npiece].k = k; pieces[npiece].W = Wk; pieces[npiece].H = Hk;
        pieces[npiece].ou = u0; pieces[npiece].ov = v0;
        pieces[npiece].du = du; pieces[npiece].dv = dv;
        pieces[npiece].valid = nv; pieces[npiece].overlap = no;
        pieces[npiece].contested = nc; pieces[npiece].nslab = nmem;
        npiece++; total_valid += nv; total_contested += nc;
        free(X); free(Y); free(Z); free(M); free(P); free(cover);
    }

    /* atlas.json */
    char ap[1200]; snprintf(ap, sizeof ap, "%s/atlas.json", out_root);
    if (ves_ensure_parent_dir(ap) != 0) failed = 1;
    FILE *af = fopen(ap, "w");
    if (af) {
        double vlo = 1e300, vhi = -1e300;
        for (size_t i = 0; i < npiece; i++) {
            if (pieces[i].ov < vlo) vlo = pieces[i].ov;
            double end_v = pieces[i].ov + pieces[i].H * pieces[i].dv;
            if (end_v > vhi) vhi = end_v;
        }
        fprintf(af, "{\n  \"format\": \"scrollfiesta-atlas-welded\", \"version\": 2,\n"
                    "  \"prefix\": \"welded\", \"n_pieces\": %zu,\n"
                    "  \"v_full\": [%.17g, %.17g], \"total_valid_px\": %" PRIu64 ", "
                    "\"total_contested_px\": %" PRIu64 ",\n"
                    "  \"verification\":{\"verified\":%s,\"lattice_compatible\":true,"
                    "\"coordinate_tolerance\":%.17g,\"max_contested_fraction\":%.17g},\n"
                    "  \"pieces\": [\n",
                npiece, npiece ? vlo : 0.0, npiece ? vhi : 0.0,
                total_valid, total_contested, failed ? "false" : "true",
                params->coord_tol, params->max_contested_fraction);
        for (size_t i = 0; i < npiece; i++) {
            char wkey[32]; wrap_name(wkey, sizeof(wkey), pieces[i].k);
            fprintf(af, "    {\"uuid\": \"welded_w%s\", \"k\": %d, \"n_slabs\": %zu, "
                        "\"W\": %d, \"H\": %d, \"origin_uv\": [%.17g, %.17g], "
                        "\"du\": %.17g, \"dv\": %.17g, \"valid_px\": %" PRIu64
                        ", \"overlap_px\": %" PRIu64 ", \"contested_px\": %" PRIu64 "}%s\n",
                    wkey, pieces[i].k, pieces[i].nslab, pieces[i].W, pieces[i].H,
                    pieces[i].ou, pieces[i].ov, pieces[i].du, pieces[i].dv,
                    pieces[i].valid, pieces[i].overlap, pieces[i].contested,
                    i + 1 < npiece ? "," : "");
        }
        fprintf(af, "  ]\n}\n");
        if (fclose(af) != 0) failed = 1;
    } else failed = 1;

    char cp[1200]; snprintf(cp, sizeof cp, "%s/assembly_certificate.json", out_root);
    if (ves_ensure_parent_dir(cp) != 0) failed = 1;
    FILE *certificate = fopen(cp, "w");
    if (certificate != NULL) {
        fprintf(certificate,
                "{\n  \"format\":\"scrollfiesta-column-assembly-certificate-v1\",\n"
                "  \"verified\":%s, \"source_roots\":%d, \"segments\":%zu,\n"
                "  \"valid_px\":%" PRIu64 ", \"contested_px\":%" PRIu64 ",\n"
                "  \"coordinate_tolerance\":%.17g, \"max_contested_fraction\":%.17g,\n"
                "  \"checks\":[\"shared_uv_lattice\",\"masked_xyz_validity\","
                "\"overlap_coordinate_agreement\",\"contested_cells_hidden\"]\n}\n",
                failed ? "false" : "true", nroot, npiece, total_valid,
                total_contested, params->coord_tol,
                params->max_contested_fraction);
        if (fclose(certificate) != 0) failed = 1;
    } else failed = 1;
    fprintf(stderr, "tifxyz_weld: %zu welded wraps, %" PRIu64
                    " valid, %" PRIu64 " contested px -> %s%s\n",
            npiece, total_valid, total_contested, out_root,
            failed ? " (FAILED VERIFICATION)" : "");
    if (out_nwrap) *out_nwrap = npiece;
    if (out_valid) *out_valid = total_valid;
    free(pieces); free(refs); Arena_free(scr);
    return failed ? 1 : 0;
}

/* ============================================================================
 * Self-test: two synthetic z-strips (disjoint v bands, shared u) each with one
 * w005 segment -> weld -> one full-height welded_w005 with exact pixel parity.
 * ==========================================================================*/
static int tw_write_strip(const char *root, double ou, double ov, int W, int H,
                          float xbase)
{
    char od[512], p[640];
    snprintf(od, sizeof od, "%s/seg/tw_w005_z%05d", root, (int)ov);
    snprintf(p, sizeof p, "%s/x.tif", od);
    if (ves_ensure_parent_dir(p) != 0) return -1;
    size_t np = (size_t)W * (size_t)H;
    float *X = (float *)malloc(np*sizeof(float)), *Y = (float *)malloc(np*sizeof(float)),
          *Z = (float *)malloc(np*sizeof(float));
    uint8_t *M = (uint8_t *)malloc(np), *P = (uint8_t *)malloc(np);
    if (!X || !Y || !Z || !M || !P) {
        free(X);free(Y);free(Z);free(M);free(P); return -1;
    }
    for (size_t i = 0; i < np; i++) {
        X[i] = xbase + (float)i; Y[i] = 20.0f; Z[i] = (float)(ov) + (float)(i / (size_t)W);
        M[i] = 255; P[i] = 1;
    }
    int rc = 0;
    rc |= TiffIO_save_float2d(p, X, W, H);
    snprintf(p, sizeof p, "%s/y.tif", od); rc |= TiffIO_save_float2d(p, Y, W, H);
    snprintf(p, sizeof p, "%s/z.tif", od); rc |= TiffIO_save_float2d(p, Z, W, H);
    snprintf(p, sizeof p, "%s/mask.tif", od); rc |= TiffIO_save(p, M, 1, H, W);
    snprintf(p, sizeof p, "%s/provenance.tif", od); rc |= TiffIO_save(p, P, 1, H, W);
    snprintf(p, sizeof p, "%s/meta.json", od);
    FILE *mf = fopen(p, "w");
    if (mf) { fprintf(mf, "{\"scale\":[1,1],\"scrollfiesta\":{\"du\":1.0,\"dv\":1.0,"
                          "\"flip_u\":0,\"flip_v\":0,\"origin_uv\":[%.1f,%.1f]}}\n", ou, ov);
              fclose(mf); } else rc = -1;
    free(X);free(Y);free(Z);free(M);free(P);
    return rc;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
        const char *base = "output/_selftest_tifxyz_weld";
        char s0[256], s1[256], out[256];
        snprintf(s0, sizeof s0, "%s/strip0", base);
        snprintf(s1, sizeof s1, "%s/strip1", base);
        snprintf(out, sizeof out, "%s/welded", base);
        int fails = 0;
        WeldParams params = {0.01, 0.0, 1e-6};
        /* strip0: u[100,104) v[0,3);  strip1: same u, v[10,13) */
        if (tw_write_strip(s0, 100.0, 0.0, 4, 3, 1000.0f) != 0) { fprintf(stderr,"FAIL write s0\n"); fails++; }
        if (tw_write_strip(s1, 100.0, 10.0, 4, 3, 2000.0f) != 0) { fprintf(stderr,"FAIL write s1\n"); fails++; }
        char *roots[2] = { s0, s1 };
        size_t nw = 0; uint64_t nv = 0;
        if (weld_run(out, roots, 2, &params, &nw, &nv) != 0) { fprintf(stderr,"FAIL weld rc\n"); fails++; }
        if (nw != 1) { fprintf(stderr, "FAIL: expected 1 welded wrap, got %zu\n", nw); fails++; }
        if (nv != 24) { fprintf(stderr, "FAIL: expected 24 valid px (2x12), got %" PRIu64 "\n", nv); fails++; }
        /* readback: welded_w005 must be 4 wide x 13 tall (v 0..12), strip0 rows
         * 0-2, gap rows 3-9 empty, strip1 rows 10-12. */
        Arena_T ar = Arena_new();
        char p[512]; uint8_t *M = NULL; int D=0,H=0,W=0;
        snprintf(p, sizeof p, "%s/seg/welded_w005/mask.tif", out);
        if (TiffIO_load(ar, p, &M, &D, &H, &W) != 0) { fprintf(stderr,"FAIL readback mask\n"); fails++; }
        else {
            if (W != 4 || H != 13) { fprintf(stderr, "FAIL: welded %dx%d, want 4x13\n", W, H); fails++; }
            int gap_empty = 1, s0_full = 1, s1_full = 1;
            for (int r = 0; r < H && r < 13; r++)
                for (int c = 0; c < 4; c++) {
                    uint8_t v = M[r*W+c];
                    if (r <= 2 && v != 255) s0_full = 0;
                    else if (r >= 3 && r <= 9 && v != 0) gap_empty = 0;
                    else if (r >= 10 && v != 255) s1_full = 0;
                }
            if (!s0_full) { fprintf(stderr, "FAIL: strip0 rows not full\n"); fails++; }
            if (!gap_empty) { fprintf(stderr, "FAIL: gap rows not empty\n"); fails++; }
            if (!s1_full) { fprintf(stderr, "FAIL: strip1 rows not full\n"); fails++; }
        }
        /* x value at strip1 row10 col0 must be strip1's xbase (2000), proving
         * absolute-v placement (not overwritten / mis-rowed). */
        float *X = NULL;
        snprintf(p, sizeof p, "%s/seg/welded_w005/x.tif", out);
        if (TiffIO_load_float2d(ar, p, &X, &W, &H) == 0 && H == 13) {
            if (fabs((double)X[0] - 1000.0) > 1e-3) { fprintf(stderr,"FAIL: row0 x!=1000\n"); fails++; }
            if (fabs((double)X[10*W] - 2000.0) > 1e-3) { fprintf(stderr,"FAIL: row10 x!=2000\n"); fails++; }
        } else { fprintf(stderr, "FAIL readback x\n"); fails++; }
        Arena_dispose(&ar);

        /* Exact halo overlap is averaged and labelled provenance=2. */
        char s2[256], s3[256], overlap_out[256];
        snprintf(s2, sizeof s2, "%s/overlap0", base);
        snprintf(s3, sizeof s3, "%s/overlap1", base);
        snprintf(overlap_out, sizeof overlap_out, "%s/overlap_welded", base);
        if (tw_write_strip(s2, 0.0, 0.0, 4, 3, 3000.0f) != 0 ||
            tw_write_strip(s3, 0.0, 0.0, 4, 3, 3000.0f) != 0) {
            fprintf(stderr, "FAIL write overlap strips\n"); fails++;
        } else {
            char *overlap_roots[2] = {s2, s3}; nw = 0; nv = 0;
            if (weld_run(overlap_out, overlap_roots, 2, &params, &nw, &nv) != 0 ||
                nw != 1 || nv != 12) {
                fprintf(stderr, "FAIL verified overlap weld\n"); fails++;
            } else {
                Arena_T oa = Arena_new(); uint8_t *provenance = NULL;
                int od=0,oh=0,ow=0;
                snprintf(p, sizeof p, "%s/seg/welded_w005/provenance.tif", overlap_out);
                if (TiffIO_load(oa, p, &provenance, &od, &oh, &ow) != 0) {
                    fprintf(stderr, "FAIL overlap provenance read\n"); fails++;
                } else for (int i = 0; i < ow*oh; i++) if (provenance[i] != 2) {
                    fprintf(stderr, "FAIL overlap provenance !=2\n"); fails++; break;
                }
                Arena_dispose(&oa);
            }
        }

        /* A geometrically conflicting overlap is hidden and fails closed. */
        char s4[256], conflict_out[256];
        snprintf(s4, sizeof s4, "%s/conflict", base);
        snprintf(conflict_out, sizeof conflict_out, "%s/conflict_welded", base);
        if (tw_write_strip(s4, 0.0, 0.0, 4, 3, 9000.0f) != 0) {
            fprintf(stderr, "FAIL write conflict strip\n"); fails++;
        } else {
            char *conflict_roots[2] = {s2, s4};
            if (weld_run(conflict_out, conflict_roots, 2, &params, NULL, NULL) == 0) {
                fprintf(stderr, "FAIL conflicting overlap was accepted\n"); fails++;
            }
        }
        fprintf(stderr, "[tifxyz_weld selftest] %s (%d failures)\n",
                fails ? "FAILED" : "PASSED", fails);
        return fails ? 1 : 0;
    }
    WeldParams params = {0.01, 0.0, 1e-6};
    int first = 1;
    while (first < argc && argv[first][0] == '-') {
        if (!strcmp(argv[first], "--coord-tol") && first+1 < argc)
            params.coord_tol = atof(argv[++first]);
        else if (!strcmp(argv[first], "--max-contested-fraction") && first+1 < argc)
            params.max_contested_fraction = atof(argv[++first]);
        else if (!strcmp(argv[first], "--lattice-tol") && first+1 < argc)
            params.lattice_tol = atof(argv[++first]);
        else { fprintf(stderr, "tifxyz_weld: unknown/incomplete option %s\n", argv[first]); return 2; }
        first++;
    }
    if (argc - first < 2 || !isfinite(params.coord_tol) || params.coord_tol < 0.0 ||
        !isfinite(params.max_contested_fraction) || params.max_contested_fraction < 0.0 ||
        params.max_contested_fraction > 1.0 || !isfinite(params.lattice_tol) ||
        params.lattice_tol <= 0.0) {
        fprintf(stderr, "usage: tifxyz_weld [--coord-tol T] [--max-contested-fraction F]\n"
                        "                    [--lattice-tol L] <out_root> <strip_root...>\n"
                        "       tifxyz_weld --selftest\n"
                        "  each strip_root holds seg/<name>/{x,y,z,mask}.tif + meta.json;\n"
                        "  same-wrap segments (*_w<k>_*) are stacked in v onto the shared\n"
                        "  u lattice into out_root/seg/welded_w###/ + atlas.json\n");
        return 2;
    }
    return weld_run(argv[first], &argv[first+1], argc - first - 1,
                    &params, NULL, NULL);
}
