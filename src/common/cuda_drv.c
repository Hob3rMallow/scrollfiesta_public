/* cuda_drv.c -- see cuda_drv.h. */
#include "common/cuda_drv.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
typedef HMODULE Lib;
static Lib lib_open(const char *path) { return LoadLibraryExA(path, NULL, LOAD_WITH_ALTERED_SEARCH_PATH); }
static void *lib_sym(Lib l, const char *name) { return (void *)GetProcAddress(l, name); }
static void lib_close(Lib l) { if (l) FreeLibrary(l); }
static const char *DRIVER_LIBS[] = {"nvcuda.dll", NULL};
static const char *NVRTC_LIBS[] = {"nvrtc64_130_0.dll", "nvrtc64_120_0.dll", "nvrtc64_112_0.dll", NULL};
static void exe_dir(char *out, size_t n)
{
    DWORD len = GetModuleFileNameA(NULL, out, (DWORD)n);
    out[len < n ? len : n - 1] = '\0';
    char *s = strrchr(out, '\\');
    if (s) s[1] = '\0'; else out[0] = '\0';
}
#define CUDA_API __stdcall
#else
#include <dlfcn.h>
#include <unistd.h>
typedef void *Lib;
static Lib lib_open(const char *path) { return dlopen(path, RTLD_NOW | RTLD_LOCAL); }
static void *lib_sym(Lib l, const char *name) { return dlsym(l, name); }
static void lib_close(Lib l) { if (l) dlclose(l); }
static const char *DRIVER_LIBS[] = {"libcuda.so.1", "libcuda.so", NULL};
static const char *NVRTC_LIBS[] = {"libnvrtc.so.13", "libnvrtc.so.12", "libnvrtc.so", NULL};
static void exe_dir(char *out, size_t n)
{
    ssize_t len = readlink("/proc/self/exe", out, n - 1);
    out[len > 0 ? (size_t)len : 0] = '\0';
    char *s = strrchr(out, '/');
    if (s) s[1] = '\0'; else out[0] = '\0';
}
#define CUDA_API
#endif

/* ---- the slice of the driver API and NVRTC this file uses ---------------- */
typedef int CUresult;
typedef int CUdevice;
typedef struct CUctx_st *CUcontext;
typedef struct CUmod_st *CUmodule;
typedef struct CUfunc_st *CUfunction;
typedef struct CUstream_st *CUstream;
typedef int nvrtcResult;
typedef struct _nvrtcProgram *nvrtcProgram;

enum { CU_ATTR_CC_MAJOR = 75, CU_ATTR_CC_MINOR = 76 };

typedef struct {
    CUresult (CUDA_API *Init)(unsigned);
    CUresult (CUDA_API *DeviceGetCount)(int *);
    CUresult (CUDA_API *DeviceGet)(CUdevice *, int);
    CUresult (CUDA_API *DeviceGetName)(char *, int, CUdevice);
    CUresult (CUDA_API *DeviceGetAttribute)(int *, int, CUdevice);
    CUresult (CUDA_API *DeviceGetPCIBusId)(char *, int, CUdevice);
    CUresult (CUDA_API *PrimaryCtxRetain)(CUcontext *, CUdevice);
    CUresult (CUDA_API *PrimaryCtxRelease)(CUdevice);
    CUresult (CUDA_API *CtxSetCurrent)(CUcontext);
    CUresult (CUDA_API *CtxSynchronize)(void);
    CUresult (CUDA_API *ModuleLoadData)(CUmodule *, const void *);
    CUresult (CUDA_API *ModuleUnload)(CUmodule);
    CUresult (CUDA_API *ModuleGetFunction)(CUfunction *, CUmodule, const char *);
    CUresult (CUDA_API *MemAlloc)(CudaDrv_ptr *, size_t);
    CUresult (CUDA_API *MemFree)(CudaDrv_ptr);
    CUresult (CUDA_API *MemcpyHtoD)(CudaDrv_ptr, const void *, size_t);
    CUresult (CUDA_API *MemcpyDtoH)(void *, CudaDrv_ptr, size_t);
    CUresult (CUDA_API *MemsetD8)(CudaDrv_ptr, unsigned char, size_t);
    CUresult (CUDA_API *MemGetInfo)(size_t *, size_t *);
    CUresult (CUDA_API *LaunchKernel)(CUfunction, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned,
                                      unsigned, CUstream, void **, void **);
    CUresult (CUDA_API *GetErrorString)(CUresult, const char **);
} DriverApi;

typedef struct {
    nvrtcResult (*CreateProgram)(nvrtcProgram *, const char *, const char *, int, const char *const *,
                                 const char *const *);
    nvrtcResult (*CompileProgram)(nvrtcProgram, int, const char *const *);
    nvrtcResult (*GetProgramLogSize)(nvrtcProgram, size_t *);
    nvrtcResult (*GetProgramLog)(nvrtcProgram, char *);
    nvrtcResult (*GetCUBINSize)(nvrtcProgram, size_t *);
    nvrtcResult (*GetCUBIN)(nvrtcProgram, char *);
    nvrtcResult (*GetPTXSize)(nvrtcProgram, size_t *);
    nvrtcResult (*GetPTX)(nvrtcProgram, char *);
    nvrtcResult (*DestroyProgram)(nvrtcProgram *);
    const char *(*GetErrorString)(nvrtcResult);
} NvrtcApi;

#define MAX_KERNELS 16

struct CudaDrv {
    Lib        driver_lib, nvrtc_lib;
    DriverApi  d;
    NvrtcApi   n;
    CUdevice   dev;
    CUcontext  ctx;
    CUmodule   mod;
    CUfunction fn[MAX_KERNELS];
    int        n_fn;
    char       device[256];
    char       err[2048];
};

static int fail(CudaDrv cu, const char *what, CUresult rc)
{
    const char *s = NULL;
    if (cu->d.GetErrorString) cu->d.GetErrorString(rc, &s);
    snprintf(cu->err, sizeof cu->err, "%s: CUDA error %d (%s)", what, rc, s ? s : "?");
    return -1;
}

#define SYM(lib, field, name)                                                                     \
    do {                                                                                          \
        *(void **)&(field) = lib_sym((lib), (name));                                              \
        if (!(field)) { snprintf(cu->err, sizeof cu->err, "missing symbol %s", (name)); return -1; } \
    } while (0)

static int load_driver(CudaDrv cu)
{
    for (int i = 0; DRIVER_LIBS[i] && !cu->driver_lib; i++) cu->driver_lib = lib_open(DRIVER_LIBS[i]);
    if (!cu->driver_lib) { snprintf(cu->err, sizeof cu->err, "the NVIDIA driver library is not loadable"); return -1; }
    Lib l = cu->driver_lib;
    SYM(l, cu->d.Init, "cuInit");
    SYM(l, cu->d.DeviceGetCount, "cuDeviceGetCount");
    SYM(l, cu->d.DeviceGet, "cuDeviceGet");
    SYM(l, cu->d.DeviceGetName, "cuDeviceGetName");
    SYM(l, cu->d.DeviceGetAttribute, "cuDeviceGetAttribute");
    SYM(l, cu->d.DeviceGetPCIBusId, "cuDeviceGetPCIBusId");
    SYM(l, cu->d.PrimaryCtxRetain, "cuDevicePrimaryCtxRetain");
    SYM(l, cu->d.PrimaryCtxRelease, "cuDevicePrimaryCtxRelease_v2");
    SYM(l, cu->d.CtxSetCurrent, "cuCtxSetCurrent");
    SYM(l, cu->d.CtxSynchronize, "cuCtxSynchronize");
    SYM(l, cu->d.ModuleLoadData, "cuModuleLoadData");
    SYM(l, cu->d.ModuleUnload, "cuModuleUnload");
    SYM(l, cu->d.ModuleGetFunction, "cuModuleGetFunction");
    SYM(l, cu->d.MemAlloc, "cuMemAlloc_v2");
    SYM(l, cu->d.MemFree, "cuMemFree_v2");
    SYM(l, cu->d.MemcpyHtoD, "cuMemcpyHtoD_v2");
    SYM(l, cu->d.MemcpyDtoH, "cuMemcpyDtoH_v2");
    SYM(l, cu->d.MemsetD8, "cuMemsetD8_v2");
    SYM(l, cu->d.MemGetInfo, "cuMemGetInfo_v2");
    SYM(l, cu->d.LaunchKernel, "cuLaunchKernel");
    SYM(l, cu->d.GetErrorString, "cuGetErrorString");
    return 0;
}

static int load_nvrtc(CudaDrv cu)
{
    char dir[1024], path[1200];
    exe_dir(dir, sizeof dir);
    for (int i = 0; NVRTC_LIBS[i] && !cu->nvrtc_lib; i++) {
        snprintf(path, sizeof path, "%s%s", dir, NVRTC_LIBS[i]);
        cu->nvrtc_lib = lib_open(path);
    }
    if (!cu->nvrtc_lib) {
        snprintf(cu->err, sizeof cu->err, "no NVRTC library beside the executable (%s%s and its nvrtc-builtins)", dir,
                 NVRTC_LIBS[0]);
        return -1;
    }
    Lib l = cu->nvrtc_lib;
    SYM(l, cu->n.CreateProgram, "nvrtcCreateProgram");
    SYM(l, cu->n.CompileProgram, "nvrtcCompileProgram");
    SYM(l, cu->n.GetProgramLogSize, "nvrtcGetProgramLogSize");
    SYM(l, cu->n.GetProgramLog, "nvrtcGetProgramLog");
    SYM(l, cu->n.GetCUBINSize, "nvrtcGetCUBINSize");
    SYM(l, cu->n.GetCUBIN, "nvrtcGetCUBIN");
    SYM(l, cu->n.GetPTXSize, "nvrtcGetPTXSize");
    SYM(l, cu->n.GetPTX, "nvrtcGetPTX");
    SYM(l, cu->n.DestroyProgram, "nvrtcDestroyProgram");
    SYM(l, cu->n.GetErrorString, "nvrtcGetErrorString");
    return 0;
}

/* Compile for the device: a cubin for sm_XY; if this NVRTC cannot target it, PTX for compute_XY (driver JIT). */
static int compile(CudaDrv cu, const char *source, int major, int minor, char **image)
{
    char arch[64];
    for (int pass = 0; pass < 2; pass++) {
        snprintf(arch, sizeof arch, "--gpu-architecture=%s_%d%d", pass == 0 ? "sm" : "compute", major, minor);
        /* exact IEEE single precision: no multiply-add contraction, correctly rounded / and sqrt, denormals kept */
        const char *opts[] = {arch, "--fmad=false", "--prec-div=true", "--prec-sqrt=true", "--ftz=false",
                              "--std=c++17"};
        nvrtcProgram prog = NULL;
        nvrtcResult rc = cu->n.CreateProgram(&prog, source, "vcr_kernels.cu", 0, NULL, NULL);
        if (rc != 0) { snprintf(cu->err, sizeof cu->err, "nvrtcCreateProgram: %s", cu->n.GetErrorString(rc)); return -1; }
        rc = cu->n.CompileProgram(prog, (int)(sizeof opts / sizeof opts[0]), opts);
        if (rc != 0) {
            size_t ls = 0;
            cu->n.GetProgramLogSize(prog, &ls);
            char *log = (char *)calloc(ls + 1, 1);
            if (log) cu->n.GetProgramLog(prog, log);
            snprintf(cu->err, sizeof cu->err, "NVRTC %s (%s): %.1800s", cu->n.GetErrorString(rc), arch, log ? log : "");
            free(log);
            cu->n.DestroyProgram(&prog);
            continue;
        }
        size_t sz = 0;
        rc = pass == 0 ? cu->n.GetCUBINSize(prog, &sz) : cu->n.GetPTXSize(prog, &sz);
        *image = (char *)malloc(sz + 1);
        if (rc == 0 && *image)
            rc = pass == 0 ? cu->n.GetCUBIN(prog, *image) : cu->n.GetPTX(prog, *image);
        cu->n.DestroyProgram(&prog);
        if (rc == 0 && *image) return 0;
        free(*image);
        *image = NULL;
    }
    return -1;
}

CudaDrv CudaDrv_open(const char *source, const char *const *kernel_names, int n_kernels, char *why, size_t why_len)
{
    CudaDrv cu = (CudaDrv)calloc(1, sizeof *cu);
    if (!cu) { snprintf(why, why_len, "out of memory"); return NULL; }
    CUresult rc = 0;
    int count = 0, major = 0, minor = 0;
    char *image = NULL;
    if (n_kernels > MAX_KERNELS) { snprintf(cu->err, sizeof cu->err, "too many kernels"); goto bad; }
    if (load_driver(cu) != 0 || load_nvrtc(cu) != 0) goto bad;
    if ((rc = cu->d.Init(0)) != 0) { fail(cu, "cuInit", rc); goto bad; }
    if ((rc = cu->d.DeviceGetCount(&count)) != 0 || count < 1) {
        snprintf(cu->err, sizeof cu->err, "no CUDA device visible (CUDA_VISIBLE_DEVICES?)");
        goto bad;
    }
    if ((rc = cu->d.DeviceGet(&cu->dev, 0)) != 0) { fail(cu, "cuDeviceGet", rc); goto bad; }
    cu->d.DeviceGetAttribute(&major, CU_ATTR_CC_MAJOR, cu->dev);
    cu->d.DeviceGetAttribute(&minor, CU_ATTR_CC_MINOR, cu->dev);
    {
        char name[160] = "?", bus[64] = "?";
        cu->d.DeviceGetName(name, (int)sizeof name, cu->dev);
        cu->d.DeviceGetPCIBusId(bus, (int)sizeof bus, cu->dev);
        snprintf(cu->device, sizeof cu->device, "%s (sm_%d%d, pci %s)", name, major, minor, bus);
    }
    if ((rc = cu->d.PrimaryCtxRetain(&cu->ctx, cu->dev)) != 0) { fail(cu, "cuDevicePrimaryCtxRetain", rc); goto bad; }
    if ((rc = cu->d.CtxSetCurrent(cu->ctx)) != 0) { fail(cu, "cuCtxSetCurrent", rc); goto bad; }
    if (compile(cu, source, major, minor, &image) != 0) goto bad;
    rc = cu->d.ModuleLoadData(&cu->mod, image);
    free(image);
    image = NULL;
    if (rc != 0) { fail(cu, "cuModuleLoadData", rc); goto bad; }
    for (int i = 0; i < n_kernels; i++) {
        if ((rc = cu->d.ModuleGetFunction(&cu->fn[i], cu->mod, kernel_names[i])) != 0) {
            fail(cu, kernel_names[i], rc);
            goto bad;
        }
    }
    cu->n_fn = n_kernels;
    return cu;
bad:
    snprintf(why, why_len, "%s", cu->err);
    free(image);
    CudaDrv_close(cu);
    return NULL;
}

void CudaDrv_close(CudaDrv cu)
{
    if (!cu) return;
    if (cu->mod) cu->d.ModuleUnload(cu->mod);
    if (cu->ctx) cu->d.PrimaryCtxRelease(cu->dev);
    lib_close(cu->nvrtc_lib);
    lib_close(cu->driver_lib);
    free(cu);
}

const char *CudaDrv_device(CudaDrv cu) { return cu->device; }
const char *CudaDrv_error(CudaDrv cu) { return cu->err; }

int CudaDrv_alloc(CudaDrv cu, CudaDrv_ptr *p, size_t bytes)
{
    CUresult rc = cu->d.MemAlloc(p, bytes ? bytes : 1);
    if (rc) { *p = 0; return fail(cu, "cuMemAlloc", rc); }
    return 0;
}

void CudaDrv_free(CudaDrv cu, CudaDrv_ptr p) { if (p) cu->d.MemFree(p); }

int CudaDrv_upload(CudaDrv cu, CudaDrv_ptr dst, const void *src, size_t bytes)
{
    CUresult rc = bytes ? cu->d.MemcpyHtoD(dst, src, bytes) : 0;
    return rc ? fail(cu, "cuMemcpyHtoD", rc) : 0;
}

int CudaDrv_download(CudaDrv cu, void *dst, CudaDrv_ptr src, size_t bytes)
{
    CUresult rc = bytes ? cu->d.MemcpyDtoH(dst, src, bytes) : 0;
    return rc ? fail(cu, "cuMemcpyDtoH", rc) : 0;
}

int CudaDrv_memset8(CudaDrv cu, CudaDrv_ptr dst, uint8_t value, size_t bytes)
{
    CUresult rc = bytes ? cu->d.MemsetD8(dst, value, bytes) : 0;
    return rc ? fail(cu, "cuMemsetD8", rc) : 0;
}

int CudaDrv_launch(CudaDrv cu, int k, unsigned gx, unsigned gy, unsigned bx, unsigned by, void **params)
{
    if (k < 0 || k >= cu->n_fn) { snprintf(cu->err, sizeof cu->err, "kernel %d not loaded", k); return -1; }
    if (gx == 0 || gy == 0) return 0;
    CUresult rc = cu->d.LaunchKernel(cu->fn[k], gx, gy, 1, bx, by, 1, 0, NULL, params, NULL);
    if (rc) return fail(cu, "cuLaunchKernel", rc);
    rc = cu->d.CtxSynchronize();
    return rc ? fail(cu, "kernel execution", rc) : 0;
}

/* ---- selftest: compile, launch, read back ---------------------------------- */
static const char *TEST_SRC =
    "extern \"C\" __global__ void t_axpy(const float *x, float *y, float a, int n)\n"
    "{ int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) y[i] = __fadd_rn(__fmul_rn(a, x[i]), y[i]); }\n";

int CudaDrv_selftest(void)
{
    char why[2048];
    const char *names[] = {"t_axpy"};
    CudaDrv cu = CudaDrv_open(TEST_SRC, names, 1, why, sizeof why);
    if (!cu) { fprintf(stderr, "cuda_drv selftest: SKIP (no device: %s)\n", why); return 0; }
    enum { N = 1000 };
    float x[N], y[N], out[N];
    for (int i = 0; i < N; i++) { x[i] = (float)i * 0.25f; y[i] = 1.0f + (float)i; }
    CudaDrv_ptr dx = 0, dy = 0;
    int ok = CudaDrv_alloc(cu, &dx, sizeof x) == 0 && CudaDrv_alloc(cu, &dy, sizeof y) == 0
          && CudaDrv_upload(cu, dx, x, sizeof x) == 0 && CudaDrv_upload(cu, dy, y, sizeof y) == 0;
    float a = 3.0f;
    int n = N;
    void *params[] = {&dx, &dy, &a, &n};
    ok = ok && CudaDrv_launch(cu, 0, (N + 255) / 256, 1, 256, 1, params) == 0
            && CudaDrv_download(cu, out, dy, sizeof out) == 0;
    for (int i = 0; ok && i < N; i++)
        if (out[i] != a * x[i] + y[i]) ok = 0;
    if (!ok) fprintf(stderr, "cuda_drv selftest: FAIL (%s)\n", CudaDrv_error(cu));
    else printf("cuda_drv selftest: PASS on %s\n", CudaDrv_device(cu));
    CudaDrv_free(cu, dx);
    CudaDrv_free(cu, dy);
    CudaDrv_close(cu);
    return ok ? 0 : 1;
}
