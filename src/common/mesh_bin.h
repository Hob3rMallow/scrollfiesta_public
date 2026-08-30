#ifndef MESH_BIN_INCLUDED
#define MESH_BIN_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "arena.h"

/* VESMESH1 is the authoritative, fixed-width mesh handoff used inside the
 * pipeline.  OBJ siblings remain interchange/debug dumps and are never a
 * fallback for these readers.  Layout is little-endian:
 *
 *   64-byte header
 *   float32 vertices[nv][3]
 *   optional float32 uv[nv][2]
 *   int32 faces[nf][3]
 */
enum { MESH_BIN_HAS_UV = 1u };

typedef struct MeshBinData {
    float *verts;       /* [nv*3] */
    float *uv;          /* NULL or [nv*2] */
    int32_t *faces;     /* [nf*3] */
    size_t nv, nf;
} MeshBinData;

/* Replace a final .obj or .vmesh suffix with .vmesh; otherwise append it. */
int MeshBin_companion_path(const char *path, char *out, size_t out_cap);

/* Cheap O(1) completeness check for resume/build-system scans.  Validates the
 * fixed header, offsets, flags, dimensions, overflow, and exact file length
 * without allocating or reading the payload.  MeshBin_write publishes by an
 * atomic rename, so a valid final-layout file is a completed write. */
int MeshBin_looks_complete(const char *path);

/* Atomic binary write.  uv == NULL writes a geometry-only container. */
int MeshBin_write(const char *path, const float *verts, size_t nv,
                  const int32_t *faces, size_t nf, const float *uv);

/* Strict readers: version/endian/layout, exact EOF, and every face index are
 * validated.  They do not inspect or fall back to an OBJ sibling. */
int MeshBin_read_arena(Arena_T arena, const char *path, MeshBinData *out);
int MeshBin_read_malloc(const char *path, MeshBinData *out);
void MeshBin_dispose(MeshBinData *mesh); /* malloc-reader results only */

int MeshBin_selftest(void);

#endif
