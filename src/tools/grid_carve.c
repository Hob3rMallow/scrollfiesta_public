/* Native input extraction only. Geometry remains in the existing C tools. */
#include "../common/zarr_grid.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    long bbox[6];
    if (argc == 2 && strcmp(argv[1],"--selftest") == 0)
        return ZarrGrid_selftest() != 0;
    if (argc != 9 && argc != 10) goto usage;
    if (argc == 10 && strcmp(argv[9],"binary") != 0) goto usage;
    for (int k = 0; k < 6; k++) {
        char *end = NULL;
        errno = 0;
        bbox[k] = strtol(argv[3+k],&end,10);
        if (errno || end == argv[3+k] || *end) goto usage;
    }
    return ZarrGrid_export(argv[1],argv[2],bbox,argc == 10) != 0;
usage:
    fprintf(stderr,"usage: grid_carve <local_zarr_array_dir> <new_cube_dir> "
                   "z0 z1 y0 y1 x0 x1 [binary]\n"
                   "       grid_carve --selftest\n"
                   "Bounds are 128-aligned, half-open original voxel coordinates.\n"
                   "Every source chunk must exist; every output TIFF is read back.\n");
    return 2;
}
