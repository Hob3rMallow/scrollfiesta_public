#![deny(unsafe_op_in_unsafe_fn)]

mod grid;
mod oracle;

#[cfg(any(
    feature = "cubecl-cpu",
    feature = "cubecl-hip",
    feature = "cubecl-cuda",
    feature = "cubecl-wgpu"
))]
mod cube_backend;

pub use grid::{Backend, MlsArena};
pub use oracle::{ProjectionResult, project as project_rust_cpu};

#[cfg(feature = "cubecl-cpu")]
pub fn project_cubecl_cpu(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, String> {
    cube_backend::project_cpu(arena, vertices, radius)
}

#[cfg(feature = "cubecl-hip")]
pub fn project_cubecl_hip(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, String> {
    cube_backend::project_hip(arena, vertices, radius)
}

#[cfg(feature = "cubecl-cuda")]
pub fn project_cubecl_cuda(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, String> {
    cube_backend::project_cuda(arena, vertices, radius)
}

#[cfg(feature = "cubecl-wgpu")]
pub fn project_cubecl_wgpu(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, String> {
    cube_backend::project_wgpu(arena, vertices, radius)
}

use libc::{c_char, c_void, size_t};
use std::cell::RefCell;
use std::ffi::{CStr, CString};
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::ptr;
use std::slice;

thread_local! {
    static LAST_ERROR: RefCell<CString> = RefCell::new(CString::new("no error").expect("literal contains no NUL"));
}

fn set_last_error(message: impl AsRef<str>) {
    let clean = message.as_ref().replace('\0', " ");
    LAST_ERROR.with(|slot| {
        *slot.borrow_mut() = CString::new(clean).expect("NUL bytes were removed");
    });
}

fn panic_message(payload: Box<dyn std::any::Any + Send>) -> String {
    if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_owned()
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else {
        "unknown Rust panic".to_owned()
    }
}

fn requested_backend(arena_backend: Backend) -> Result<Backend, String> {
    if arena_backend != Backend::Auto {
        return Ok(arena_backend);
    }
    if let Ok(value) = std::env::var("MLS_BACKEND") {
        return match value.trim().to_ascii_lowercase().as_str() {
            "auto" => Ok(Backend::Auto),
            "rust" | "rust-cpu" | "cpu" => Ok(Backend::RustCpu),
            "cube-cpu" | "cubecl-cpu" => Ok(Backend::CubeCpu),
            "hip" | "cube-hip" | "cubecl-hip" => Ok(Backend::CubeHip),
            "cuda" | "cube-cuda" | "cubecl-cuda" => Ok(Backend::CubeCuda),
            "wgpu" | "cube-wgpu" | "cubecl-wgpu" => Ok(Backend::CubeWgpu),
            other => Err(format!("unknown MLS_BACKEND={other:?}")),
        };
    }
    Ok(Backend::Auto)
}

fn auto_backend() -> Backend {
    #[cfg(feature = "cubecl-hip")]
    {
        return Backend::CubeHip;
    }
    #[cfg(all(not(feature = "cubecl-hip"), feature = "cubecl-cuda"))]
    {
        return Backend::CubeCuda;
    }
    // cubecl-wgpu is deliberately NOT auto-selected: its output is known
    // incorrect on both the SPIR-V and WGSL code-generation paths (see
    // WGPU_STATUS.md). It is reachable only via an explicit MLS_BACKEND=wgpu /
    // Backend::CubeWgpu request, kept for debugging and upstream regression use.
    #[cfg(all(
        not(feature = "cubecl-hip"),
        not(feature = "cubecl-cuda"),
        feature = "cubecl-cpu"
    ))]
    {
        return Backend::CubeCpu;
    }
    #[allow(unreachable_code)]
    Backend::RustCpu
}

fn dispatch_with_normals(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
    compute_normals: bool,
) -> Result<ProjectionResult, String> {
    let mut backend = requested_backend(arena.backend)?;
    if backend == Backend::Auto {
        backend = auto_backend();
    }
    match backend {
        Backend::Auto => unreachable!("auto backend must be resolved"),
        Backend::RustCpu => oracle::project_with_normals(arena, vertices, radius, compute_normals)
            .map_err(str::to_owned),
        Backend::CubeCpu => {
            #[cfg(feature = "cubecl-cpu")]
            {
                cube_backend::project_cpu_with_normals(arena, vertices, radius, compute_normals)
            }
            #[cfg(not(feature = "cubecl-cpu"))]
            {
                Err("CubeCL CPU backend was not compiled; build with --features cubecl-cpu".into())
            }
        }
        Backend::CubeHip => {
            #[cfg(feature = "cubecl-hip")]
            {
                cube_backend::project_hip_with_normals(arena, vertices, radius, compute_normals)
            }
            #[cfg(not(feature = "cubecl-hip"))]
            {
                Err("CubeCL HIP backend was not compiled; build with --features cubecl-hip".into())
            }
        }
        Backend::CubeCuda => {
            #[cfg(feature = "cubecl-cuda")]
            {
                cube_backend::project_cuda_with_normals(arena, vertices, radius, compute_normals)
            }
            #[cfg(not(feature = "cubecl-cuda"))]
            {
                Err(
                    "CubeCL CUDA backend was not compiled; build with --features cubecl-cuda"
                        .into(),
                )
            }
        }
        Backend::CubeWgpu => {
            #[cfg(feature = "cubecl-wgpu")]
            {
                cube_backend::project_wgpu_with_normals(arena, vertices, radius, compute_normals)
            }
            #[cfg(not(feature = "cubecl-wgpu"))]
            {
                Err(
                    "CubeCL WGPU backend was not compiled; build with --features cubecl-wgpu"
                        .into(),
                )
            }
        }
    }
}

fn failure_outputs(
    vertices: &[f32],
    float_count: usize,
    out_vertices: *mut f32,
    out_normals: *mut f32,
) {
    if !out_vertices.is_null() && float_count != 0 {
        // SAFETY: the FFI contract requires output buffers to hold float_count floats.
        unsafe {
            if vertices.len() == float_count {
                ptr::copy(vertices.as_ptr(), out_vertices, float_count);
            } else {
                ptr::write_bytes(out_vertices, 0, float_count);
            }
        };
    }
    if !out_normals.is_null() && float_count != 0 {
        // SAFETY: the FFI contract requires output buffers to hold float_count floats.
        unsafe { ptr::write_bytes(out_normals, 0, float_count) };
    }
}

fn scrollfiesta_origin_xyz(
    points_xyz: &[[f32; 3]],
    radius: f32,
    cell_origin_zyx: [f32; 3],
) -> [f32; 3] {
    let mut minimum = points_xyz[0];
    for point in &points_xyz[1..] {
        for axis in 0..3 {
            minimum[axis] = minimum[axis].min(point[axis]);
        }
    }
    let additive_xyz = [cell_origin_zyx[2], cell_origin_zyx[1], cell_origin_zyx[0]];
    let mut origin = [0.0_f32; 3];
    for axis in 0..3 {
        let phase_origin = -additive_xyz[axis];
        let cell_shift = ((minimum[axis] - phase_origin) / radius).floor() - 1.0;
        origin[axis] = phase_origin + cell_shift * radius;
    }
    origin
}

/// Perform one production-equivalent ScrollFiesta projection pass.
///
/// The input and output use ScrollFiesta's z-y-x coordinate order. The support
/// grid is rebuilt from `vertices_zyx`, which is both the query and support
/// cloud for this pass. Caller-owned 5/20-pass loops must call this function
/// again with the previous output.
pub fn project_scrollfiesta_one_pass(
    vertices_zyx: &[f32],
    radius: f32,
    cell_origin_zyx: [f32; 3],
    backend: Backend,
    compute_normals: bool,
) -> Result<ProjectionResult, String> {
    if vertices_zyx.is_empty() || !vertices_zyx.len().is_multiple_of(3) {
        return Err("vertices must contain one or more zyx triples".into());
    }
    if !vertices_zyx.iter().all(|value| value.is_finite()) {
        return Err("vertices contain a non-finite value".into());
    }
    if !radius.is_finite() || radius <= 0.0 {
        return Err("radius must be finite and > 0".into());
    }
    if !cell_origin_zyx.iter().all(|value| value.is_finite()) {
        return Err("cell_origin contains a non-finite value".into());
    }

    let points_xyz: Vec<[f32; 3]> = vertices_zyx
        .chunks_exact(3)
        .map(|point| [point[2], point[1], point[0]])
        .collect();
    let origin_xyz = scrollfiesta_origin_xyz(&points_xyz, radius, cell_origin_zyx);
    let mut arena = MlsArena::from_points(&points_xyz, radius, Some(origin_xyz))
        .map_err(|error| error.to_string())?;
    arena.backend = backend;
    let vertices_xyz: Vec<f32> = points_xyz.iter().flatten().copied().collect();
    let projected = dispatch_with_normals(&arena, &vertices_xyz, radius, compute_normals)?;

    let mut vertices = Vec::with_capacity(projected.vertices.len());
    for point in projected.vertices.chunks_exact(3) {
        vertices.extend_from_slice(&[point[2], point[1], point[0]]);
    }
    let mut normals = Vec::with_capacity(projected.normals.len());
    for normal in projected.normals.chunks_exact(3) {
        normals.extend_from_slice(&[normal[2], normal[1], normal[0]]);
    }
    Ok(ProjectionResult { vertices, normals })
}

/// Returns a thread-local human-readable error string. The pointer remains valid
/// until the next library call on the same thread.
#[unsafe(no_mangle)]
pub extern "C" fn MLS_cubecl_last_error() -> *const c_char {
    LAST_ERROR.with(|slot| slot.borrow().as_ptr())
}

/// Bitmask: bit 0 Rust CPU, bit 1 CubeCL CPU, bit 2 CubeCL HIP, bit 3 CubeCL CUDA, bit 4 CubeCL WGPU.
#[unsafe(no_mangle)]
pub extern "C" fn MLS_cubecl_available_backends() -> u32 {
    let mask = 1_u32;
    #[cfg(feature = "cubecl-cpu")]
    let mask = mask | (1 << 1);
    #[cfg(feature = "cubecl-hip")]
    let mask = mask | (1 << 2);
    #[cfg(feature = "cubecl-cuda")]
    let mask = mask | (1 << 3);
    #[cfg(feature = "cubecl-wgpu")]
    let mask = mask | (1 << 4);
    mask
}

/// Build a Rust-owned counting-sort arena from xyz point triples.
///
/// # Safety
///
/// `points_xyz` must reference `npoints * 3` readable floats. When non-null,
/// `cell_origin_or_null` must reference three readable floats.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn MLS_cubecl_arena_create_from_points(
    points_xyz: *const f32,
    npoints: size_t,
    cell_size: f32,
    cell_origin_or_null: *const f32,
) -> *mut c_void {
    let result = catch_unwind(AssertUnwindSafe(|| -> Result<*mut c_void, String> {
        if points_xyz.is_null() || npoints == 0 {
            return Err("points_xyz is null or npoints is zero".into());
        }
        let float_count = npoints.checked_mul(3).ok_or("point count overflow")?;
        // SAFETY: caller promises npoints xyz triples.
        let flat = unsafe { slice::from_raw_parts(points_xyz, float_count) };
        let points: Vec<[f32; 3]> = flat.chunks_exact(3).map(|p| [p[0], p[1], p[2]]).collect();
        let origin = if cell_origin_or_null.is_null() {
            None
        } else {
            // SAFETY: non-null origin points to three floats by contract.
            let o = unsafe { slice::from_raw_parts(cell_origin_or_null, 3) };
            Some([o[0], o[1], o[2]])
        };
        let arena = MlsArena::from_points(&points, cell_size, origin).map_err(|e| e.to_string())?;
        Ok(Box::into_raw(Box::new(arena)).cast())
    }));

    match result {
        Ok(Ok(pointer)) => {
            set_last_error("no error");
            pointer
        }
        Ok(Err(error)) => {
            set_last_error(error);
            ptr::null_mut()
        }
        Err(payload) => {
            set_last_error(format!(
                "panic while creating arena: {}",
                panic_message(payload)
            ));
            ptr::null_mut()
        }
    }
}

/// Copy an existing counting-sort cell grid into a Rust-owned arena for the
/// explicit static-support API. ScrollFiesta's legacy path does not use this;
/// it rebuilds the grid from current vertices on every call.
///
/// # Safety
///
/// Every pointer must reference the number of readable elements specified by
/// its associated count. `dims_xyz` and `cell_origin` must each reference
/// three readable values.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn MLS_cubecl_arena_create_from_grid(
    points_xyz: *const f32,
    npoints: size_t,
    cell_offsets: *const u32,
    cell_offsets_len: size_t,
    point_indices: *const u32,
    point_indices_len: size_t,
    dims_xyz: *const u32,
    cell_origin: *const f32,
    cell_size: f32,
) -> *mut c_void {
    let result = catch_unwind(AssertUnwindSafe(|| -> Result<*mut c_void, String> {
        if points_xyz.is_null()
            || cell_offsets.is_null()
            || point_indices.is_null()
            || dims_xyz.is_null()
            || cell_origin.is_null()
            || npoints == 0
        {
            return Err("one or more required grid pointers are null".into());
        }
        let point_float_count = npoints.checked_mul(3).ok_or("point count overflow")?;
        // SAFETY: lengths and pointers are supplied by the caller under the C contract.
        let points = unsafe { slice::from_raw_parts(points_xyz, point_float_count) }.to_vec();
        let offsets = unsafe { slice::from_raw_parts(cell_offsets, cell_offsets_len) }.to_vec();
        let indices = unsafe { slice::from_raw_parts(point_indices, point_indices_len) }.to_vec();
        let dims_slice = unsafe { slice::from_raw_parts(dims_xyz, 3) };
        let origin_slice = unsafe { slice::from_raw_parts(cell_origin, 3) };
        let arena = MlsArena::from_grid(
            points,
            offsets,
            indices,
            [dims_slice[0], dims_slice[1], dims_slice[2]],
            [origin_slice[0], origin_slice[1], origin_slice[2]],
            cell_size,
        )
        .map_err(|e| e.to_string())?;
        Ok(Box::into_raw(Box::new(arena)).cast())
    }));

    match result {
        Ok(Ok(pointer)) => {
            set_last_error("no error");
            pointer
        }
        Ok(Err(error)) => {
            set_last_error(error);
            ptr::null_mut()
        }
        Err(payload) => {
            set_last_error(format!(
                "panic while importing arena: {}",
                panic_message(payload)
            ));
            ptr::null_mut()
        }
    }
}

/// Destroy an arena returned by one of this library's constructors.
///
/// # Safety
///
/// `arena` must be null or a live, uniquely owned arena pointer returned by
/// this library. A non-null pointer may be destroyed exactly once.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn MLS_cubecl_arena_destroy(arena: *mut c_void) {
    if arena.is_null() {
        return;
    }
    // SAFETY: pointer must originate from one of the arena constructors and be destroyed once.
    unsafe { drop(Box::from_raw(arena.cast::<MlsArena>())) };
}

/// Select the backend used by subsequent projections on an arena.
///
/// # Safety
///
/// `arena` must be a live, uniquely borrowed arena pointer returned by this
/// library for the duration of the call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn MLS_cubecl_arena_set_backend(arena: *mut c_void, backend: u32) -> i32 {
    let result = catch_unwind(AssertUnwindSafe(|| -> Result<(), String> {
        if arena.is_null() {
            return Err("arena pointer is null".into());
        }
        let backend = Backend::try_from(backend).map_err(|e| e.to_string())?;
        // SAFETY: caller supplied an arena created by this library.
        let arena = unsafe { &mut *arena.cast::<MlsArena>() };
        if !arena.validate_tag() {
            return Err("arena magic/version mismatch".into());
        }
        arena.backend = backend;
        Ok(())
    }));

    match result {
        Ok(Ok(())) => {
            set_last_error("no error");
            0
        }
        Ok(Err(error)) => {
            set_last_error(error);
            -1
        }
        Err(payload) => {
            set_last_error(format!(
                "panic while setting backend: {}",
                panic_message(payload)
            ));
            -1
        }
    }
}

/// Project xyz query vertices against an explicitly supplied Rust-owned arena.
/// This is the static-support API for standalone benchmarks and imported grids.
///
/// # Safety
///
/// `arena` must be a live arena returned by this library. The other non-null
/// pointers must reference `nv * 3` floats. The input may alias either output;
/// the two output ranges must not overlap each other.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn MLS_cubecl_project_with_arena(
    arena: *mut c_void,
    verts_xyz: *const f32,
    nv: size_t,
    radius_vox: f32,
    out_verts_xyz: *mut f32,
    out_normals_xyz: *mut f32,
) {
    let float_count = match nv.checked_mul(3) {
        Some(value) => value,
        None => {
            set_last_error("vertex count overflow");
            return;
        }
    };
    if nv == 0 {
        set_last_error("no error");
        return;
    }
    let vertices = if verts_xyz.is_null() {
        &[]
    } else {
        // SAFETY: caller promises nv xyz triples.
        unsafe { slice::from_raw_parts(verts_xyz, float_count) }
    };
    let result = catch_unwind(AssertUnwindSafe(|| -> Result<ProjectionResult, String> {
        if arena.is_null() || verts_xyz.is_null() || out_verts_xyz.is_null() {
            return Err("arena, input, or vertex output pointer is null".into());
        }
        // SAFETY: caller supplied a live arena created by this library.
        let arena = unsafe { &*arena.cast::<MlsArena>() };
        if !arena.validate_tag() {
            return Err("arena magic/version mismatch".into());
        }
        dispatch_with_normals(arena, vertices, radius_vox, !out_normals_xyz.is_null())
    }));
    match result {
        Ok(Ok(projected)) => {
            // SAFETY: output lengths are nv*3 under the C contract. `copy`
            // permits the in-place aliasing supported by ScrollFiesta callers.
            unsafe {
                ptr::copy(projected.vertices.as_ptr(), out_verts_xyz, float_count);
                if !out_normals_xyz.is_null() {
                    ptr::copy(projected.normals.as_ptr(), out_normals_xyz, float_count);
                }
            }
            set_last_error("no error");
        }
        Ok(Err(error)) => {
            failure_outputs(vertices, float_count, out_verts_xyz, out_normals_xyz);
            set_last_error(error);
        }
        Err(payload) => {
            failure_outputs(vertices, float_count, out_verts_xyz, out_normals_xyz);
            set_last_error(format!(
                "panic in static arena projection: {}",
                panic_message(payload)
            ));
        }
    }
}

/// Exact one-pass legacy symbol and parameter order expected by ScrollFiesta.
///
/// `arena` is ScrollFiesta scratch allocation state and is intentionally not
/// dereferenced. `verts` is both the z-y-x query cloud and the support cloud;
/// the counting-sort grid is rebuilt on every call. The caller retains its
/// existing moving-support 5/20-pass loop.
///
/// # Safety
///
/// For `nv > 0`, `verts` and `out_verts` must reference `nv * 3` floats.
/// `out_normals` may be null. `cell_origin`, when non-null, references three
/// floats. The input may alias either output; the two output ranges must not
/// overlap each other.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn MLS_project_verts(
    arena: *mut c_void,
    verts: *const f32,
    nv: size_t,
    radius_vox: f32,
    cell_origin: *const f32,
    out_verts: *mut f32,
    out_normals: *mut f32,
) {
    let _ = arena;
    let float_count = match nv.checked_mul(3) {
        Some(value) => value,
        None => {
            set_last_error("vertex count overflow");
            return;
        }
    };
    if nv == 0 {
        set_last_error("no error");
        return;
    }
    let vertices = if verts.is_null() {
        &[]
    } else {
        // SAFETY: caller promises nv zyx triples.
        unsafe { slice::from_raw_parts(verts, float_count) }
    };
    let result = catch_unwind(AssertUnwindSafe(|| -> Result<ProjectionResult, String> {
        if verts.is_null() || out_verts.is_null() {
            return Err("input or vertex output pointer is null".into());
        }
        if nv > u32::MAX as usize {
            return Err("CubeCL uses u32 indexing; nv exceeds u32::MAX".into());
        }
        let origin = if cell_origin.is_null() {
            [0.0_f32; 3]
        } else {
            // SAFETY: non-null origin points to three floats by contract.
            let supplied = unsafe { slice::from_raw_parts(cell_origin, 3) };
            [supplied[0], supplied[1], supplied[2]]
        };
        project_scrollfiesta_one_pass(
            vertices,
            radius_vox,
            origin,
            Backend::Auto,
            !out_normals.is_null(),
        )
    }));
    match result {
        Ok(Ok(projected)) => {
            // SAFETY: output lengths are nv*3 under the legacy ABI contract.
            unsafe {
                ptr::copy(projected.vertices.as_ptr(), out_verts, float_count);
                if !out_normals.is_null() {
                    ptr::copy(projected.normals.as_ptr(), out_normals, float_count);
                }
            }
            set_last_error("no error");
        }
        Ok(Err(error)) => {
            failure_outputs(vertices, float_count, out_verts, out_normals);
            set_last_error(error);
        }
        Err(payload) => {
            failure_outputs(vertices, float_count, out_verts, out_normals);
            set_last_error(format!(
                "panic in MLS_project_verts: {}",
                panic_message(payload)
            ));
        }
    }
}

/// Optional helper for diagnostics from C/C++.
#[unsafe(no_mangle)]
pub extern "C" fn MLS_cubecl_backend_name(backend: u32) -> *const c_char {
    static AUTO: &[u8] = b"auto\0";
    static RUST_CPU: &[u8] = b"rust-cpu\0";
    static CUBE_CPU: &[u8] = b"cubecl-cpu\0";
    static CUBE_HIP: &[u8] = b"cubecl-hip\0";
    static CUBE_CUDA: &[u8] = b"cubecl-cuda\0";
    static CUBE_WGPU: &[u8] = b"cubecl-wgpu\0";
    static UNKNOWN: &[u8] = b"unknown\0";
    match backend {
        0 => AUTO.as_ptr().cast(),
        1 => RUST_CPU.as_ptr().cast(),
        2 => CUBE_CPU.as_ptr().cast(),
        3 => CUBE_HIP.as_ptr().cast(),
        4 => CUBE_CUDA.as_ptr().cast(),
        5 => CUBE_WGPU.as_ptr().cast(),
        _ => UNKNOWN.as_ptr().cast(),
    }
}

/// Parse backend names in scripts without duplicating the enum table.
///
/// # Safety
///
/// `name` must point to a readable NUL-terminated C string for the duration of
/// the call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn MLS_cubecl_backend_from_name(name: *const c_char) -> i32 {
    if name.is_null() {
        set_last_error("backend name is null");
        return -1;
    }
    // SAFETY: caller supplies a NUL-terminated C string.
    let text = match unsafe { CStr::from_ptr(name) }.to_str() {
        Ok(value) => value.trim().to_ascii_lowercase(),
        Err(error) => {
            set_last_error(format!("backend name is not UTF-8: {error}"));
            return -1;
        }
    };
    match text.as_str() {
        "auto" => 0,
        "rust" | "rust-cpu" | "cpu" => 1,
        "cube-cpu" | "cubecl-cpu" => 2,
        "hip" | "cube-hip" | "cubecl-hip" => 3,
        "cuda" | "cube-cuda" | "cubecl-cuda" => 4,
        "wgpu" | "cube-wgpu" | "cubecl-wgpu" => 5,
        _ => {
            set_last_error(format!("unknown backend name {text:?}"));
            -1
        }
    }
}
