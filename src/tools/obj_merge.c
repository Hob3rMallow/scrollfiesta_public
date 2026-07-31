/* obj_merge.c -- plain-concatenate OBJ meshes into one, renumbering face indices
 * by the running vertex count. No offset, no weld, no dedup -- the inputs are
 * already in a shared (world) coordinate frame. Vertex colours, if present, are
 * dropped (ObjIO reads/writes plain xyz). Used to bucket gate-classified sheets
 * into one "good" / one "bad" OBJ.
 *
 * Usage:
 *   obj_merge <out.obj> <in1.obj> [in2.obj ...]
 *   obj_merge --selftest
 */
#include "../common/ves_platform.h"
#include "../common/arena.h"
#include "../common/obj_io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* Merge `n` OBJ files (paths) into out. Returns 0 on success. *out_nv/*out_nf
 * receive the merged counts. Inputs that fail to read are skipped (warned). */
static int merge_files(const char *out, char **paths, int n,
                       size_t *out_nv, size_t *out_nf, int *out_nread)
{
    float   *V = NULL; int32_t *F = NULL;
    size_t   nv = 0, nf = 0, capv = 0, capf = 0;
    int      nread = 0;
    for (int a = 0; a < n; a++) {
        Arena_T ar = Arena_new();
        float *v = NULL; int32_t *f = NULL; size_t lnv = 0, lnf = 0;
        if (ObjIO_read(ar, paths[a], &v, &lnv, &f, &lnf) != 0) {
            fprintf(stderr, "obj_merge: skip (unreadable) %s\n", paths[a]);
            Arena_dispose(&ar); continue;
        }
        if (nv + lnv > capv) { capv = (nv + lnv) * 2 + 1024;
            V = (float*)realloc(V, capv * 3 * sizeof(float)); }
        memcpy(V + nv * 3, v, lnv * 3 * sizeof(float));
        if (nf + lnf > capf) { capf = (nf + lnf) * 2 + 1024;
            F = (int32_t*)realloc(F, capf * 3 * sizeof(int32_t)); }
        for (size_t i = 0; i < lnf; i++) {
            F[(nf + i) * 3 + 0] = f[i * 3 + 0] + (int32_t)nv;
            F[(nf + i) * 3 + 1] = f[i * 3 + 1] + (int32_t)nv;
            F[(nf + i) * 3 + 2] = f[i * 3 + 2] + (int32_t)nv;
        }
        nv += lnv; nf += lnf; nread++;
        Arena_dispose(&ar);
    }
    ves_ensure_parent_dir(out);
    int rc = ObjIO_write(out, V, nv, F, nf);
    free(V); free(F);
    if (out_nv) *out_nv = nv;
    if (out_nf) *out_nf = nf;
    if (out_nread) *out_nread = nread;
    return rc;
}

static int run_selftest(void){
    int fail = 0;
    const char *pa = "output/_objmerge_test_a.obj";
    const char *pb = "output/_objmerge_test_b.obj";
    const char *po = "output/_objmerge_test_out.obj";
    float va[9] = {0,0,0, 1,0,0, 0,1,0};      int32_t fa[3] = {0,1,2};
    float vb[9] = {5,5,5, 6,5,5, 5,6,5};      int32_t fb[3] = {0,1,2};
    ves_ensure_parent_dir(pa);
    ObjIO_write(pa, va, 3, fa, 1);
    ObjIO_write(pb, vb, 3, fb, 1);
    char *paths[2]; paths[0] = (char*)pa; paths[1] = (char*)pb;
    size_t nv=0, nf=0; int nread=0;
    merge_files(po, paths, 2, &nv, &nf, &nread);
    /* read back and verify second triangle's indices were offset by 3 */
    Arena_T ar = Arena_new();
    float *V=NULL; int32_t *F=NULL; size_t rnv=0, rnf=0;
    int rc = ObjIO_read(ar, po, &V, &rnv, &F, &rnf);
    int f1 = (rc != 0 || rnv != 6 || rnf != 2);
    int f2 = !f1 && !(F[3]==3 && F[4]==4 && F[5]==5);   /* 2nd tri offset */
    fail = f1 || f2;
    fprintf(stderr,"[selftest] merge 2 tris -> %zu v %zu f (want 6/2); 2nd tri=(%d,%d,%d) want(3,4,5) -> %s\n",
            rnv, rnf, f1?-1:F[3], f1?-1:F[4], f1?-1:F[5], fail?"FAIL":"ok");
    Arena_dispose(&ar);
    remove(pa); remove(pb); remove(po);
    fprintf(stderr,"=== obj_merge selftest %s ===\n", fail?"FAILED":"PASSED");
    return fail?1:0;
}

int main(int argc, char **argv){
    if(argc>=2 && !strcmp(argv[1],"--selftest")) return run_selftest();
    if(argc<3){
        fprintf(stderr,"Usage: %s <out.obj> <in1.obj> [in2.obj ...]\n"
                       "       %s --selftest\n", argv[0], argv[0]);
        return 1;
    }
    size_t nv=0, nf=0; int nread=0;
    int rc = merge_files(argv[1], &argv[2], argc-2, &nv, &nf, &nread);
    fprintf(stderr,"obj_merge: %d/%d inputs -> %zu v %zu f -> %s\n",
            nread, argc-2, nv, nf, argv[1]);
    return rc?1:0;
}
