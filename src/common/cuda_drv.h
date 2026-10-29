#ifndef CUDA_DRV_INCLUDED
#define CUDA_DRV_INCLUDED

#include <stddef.h>
#include <stdint.h>

/* ============================================================================
 * cuda_drv.h -- the CUDA driver API and NVRTC, loaded at run time.
 *
 * No CUDA toolkit is needed to build or run: nvcuda.dll comes with the
 * NVIDIA driver, and the runtime compiler (nvrtc64_*.dll with its
 * nvrtc-builtins64_*.dll) is loaded from the executable's own directory.
 * Kernels are CUDA C source strings compiled once per process for the
 * device's own architecture (a cubin for sm_XY).
 *
 * The device is ordinal 0 of the devices the driver lets this process see,
 * so CUDA_VISIBLE_DEVICES (with CUDA_DEVICE_ORDER=PCI_BUS_ID) picks the card,
 * as for any other CUDA program.
 * ==========================================================================*/

typedef unsigned long long CudaDrv_ptr;      /* CUdeviceptr */

typedef struct CudaDrv *CudaDrv;

/* Load the driver and NVRTC, create the device's primary context, compile
 * `source` (options: no FMA contraction, IEEE division and square root) and
 * fetch the named kernels.  On failure returns NULL and writes a reason. */
CudaDrv CudaDrv_open(const char *source, const char *const *kernel_names, int n_kernels,
                     char *why, size_t why_len);
void    CudaDrv_close(CudaDrv cu);

/* "<name> (sm_XY, pci <bus id>)" of the device in use. */
const char *CudaDrv_device(CudaDrv cu);

int  CudaDrv_alloc(CudaDrv cu, CudaDrv_ptr *p, size_t bytes);
void CudaDrv_free(CudaDrv cu, CudaDrv_ptr p);
int  CudaDrv_upload(CudaDrv cu, CudaDrv_ptr dst, const void *src, size_t bytes);
int  CudaDrv_download(CudaDrv cu, void *dst, CudaDrv_ptr src, size_t bytes);
int  CudaDrv_memset8(CudaDrv cu, CudaDrv_ptr dst, uint8_t value, size_t bytes);

/* Launch kernel `k` (index into kernel_names) on a grid of blocks; params is
 * the usual array of pointers to each argument.  Synchronizes and reports
 * any launch or execution error. */
int CudaDrv_launch(CudaDrv cu, int k, unsigned gx, unsigned gy, unsigned bx, unsigned by, void **params);

/* The last error text (driver, NVRTC or loader). */
const char *CudaDrv_error(CudaDrv cu);

int CudaDrv_selftest(void);

#endif
