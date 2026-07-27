//! One CubeCL kernel shared by the CPU and HIP runtimes.
//!
//! The kernel performs exactly one MLS projection. ScrollFiesta owns the
//! moving-support 5/20-pass loop and rebuilds the grid between calls.

use crate::grid::MlsArena;
use crate::oracle::ProjectionResult;
use cubecl::prelude::*;

const CUBE_DIM: u32 = 256;

#[cube]
fn compensated_add(value: f32, sum: &mut f32, compensation: &mut f32) {
    let corrected = value - *compensation;
    let next = *sum + corrected;
    *compensation = (next - *sum) - corrected;
    *sum = next;
}

#[cube(launch_unchecked)]
fn mls_iteration(
    points: &Array<f32>,
    cell_offsets: &Array<u32>,
    point_indices: &Array<u32>,
    input_vertices: &Array<f32>,
    output_vertices: &mut Array<f32>,
    output_normals: &mut Array<f32>,
    meta_f32: &Array<f32>,
    meta_u32: &Array<u32>,
) {
    let vertex_index = ABSOLUTE_POS;
    let nvertices = meta_u32[3usize] as usize;

    if vertex_index < nvertices {
        let vertex_base = vertex_index * 3usize;
        let px = input_vertices[vertex_base];
        let py = input_vertices[vertex_base + 1usize];
        let pz = input_vertices[vertex_base + 2usize];

        let origin_x = meta_f32[0usize];
        let origin_y = meta_f32[1usize];
        let origin_z = meta_f32[2usize];
        let cell_size = meta_f32[3usize];
        let radius = meta_f32[4usize];
        let radius2 = radius * radius;

        let dim_x = meta_u32[0usize] as usize;
        let dim_y = meta_u32[1usize] as usize;
        let dim_z = meta_u32[2usize] as usize;
        let reach = meta_u32[4usize] as usize;

        let mut raw_x = (px - origin_x) / cell_size;
        let mut raw_y = (py - origin_y) / cell_size;
        let mut raw_z = (pz - origin_z) / cell_size;
        if raw_x < 0.0 {
            raw_x = 0.0;
        }
        if raw_y < 0.0 {
            raw_y = 0.0;
        }
        if raw_z < 0.0 {
            raw_z = 0.0;
        }
        let max_x = f32::cast_from(dim_x - 1usize);
        let max_y = f32::cast_from(dim_y - 1usize);
        let max_z = f32::cast_from(dim_z - 1usize);
        if raw_x > max_x {
            raw_x = max_x;
        }
        if raw_y > max_y {
            raw_y = max_y;
        }
        if raw_z > max_z {
            raw_z = max_z;
        }
        let cx = usize::cast_from(raw_x.floor());
        let cy = usize::cast_from(raw_y.floor());
        let cz = usize::cast_from(raw_z.floor());

        let mut x0 = 0usize;
        let mut y0 = 0usize;
        let mut z0 = 0usize;
        if cx > reach {
            x0 = cx - reach;
        }
        if cy > reach {
            y0 = cy - reach;
        }
        if cz > reach {
            z0 = cz - reach;
        }
        let x1 = if reach >= dim_x - 1usize - cx {
            dim_x - 1usize
        } else {
            cx + reach
        };
        let y1 = if reach >= dim_y - 1usize - cy {
            dim_y - 1usize
        } else {
            cy + reach
        };
        let z1 = if reach >= dim_z - 1usize - cz {
            dim_z - 1usize
        } else {
            cz + reach
        };

        let mut sum_w = 0.0f32;
        let mut centroid_x = 0.0f32;
        let mut centroid_y = 0.0f32;
        let mut centroid_z = 0.0f32;
        let mut compensation_w = 0.0f32;
        let mut compensation_x = 0.0f32;
        let mut compensation_y = 0.0f32;
        let mut compensation_z = 0.0f32;
        let mut neighbour_count: u32 = 0;

        let mut z = z0;
        while z <= z1 {
            let mut y = y0;
            while y <= y1 {
                let mut x = x0;
                while x <= x1 {
                    let cell = x + dim_x * (y + dim_y * z);
                    let mut slot = cell_offsets[cell] as usize;
                    let end = cell_offsets[cell + 1usize] as usize;
                    while slot < end {
                        let point_index = point_indices[slot] as usize;
                        let point_base = point_index * 3usize;
                        let qx = points[point_base];
                        let qy = points[point_base + 1usize];
                        let qz = points[point_base + 2usize];
                        let dx = qx - px;
                        let dy = qy - py;
                        let dz = qz - pz;
                        let distance2 = dx * dx + dy * dy + dz * dz;
                        if distance2 < radius2 {
                            let ratio = distance2.sqrt() / radius;
                            let one_minus = 1.0 - ratio;
                            let square = one_minus * one_minus;
                            let weight = square * square * (4.0 * ratio + 1.0);
                            compensated_add(weight, &mut sum_w, &mut compensation_w);
                            compensated_add(
                                weight * (qx - px),
                                &mut centroid_x,
                                &mut compensation_x,
                            );
                            compensated_add(
                                weight * (qy - py),
                                &mut centroid_y,
                                &mut compensation_y,
                            );
                            compensated_add(
                                weight * (qz - pz),
                                &mut centroid_z,
                                &mut compensation_z,
                            );
                            neighbour_count += 1u32;
                        }
                        slot += 1usize;
                    }
                    x += 1usize;
                }
                y += 1usize;
            }
            z += 1usize;
        }

        if neighbour_count >= 4u32 && sum_w > 1.0e-12 {
            centroid_x = px + centroid_x / sum_w;
            centroid_y = py + centroid_y / sum_w;
            centroid_z = pz + centroid_z / sum_w;

            if meta_u32[5usize] == 0u32 {
                output_vertices[vertex_base] = centroid_x;
                output_vertices[vertex_base + 1usize] = centroid_y;
                output_vertices[vertex_base + 2usize] = centroid_z;
                output_normals[vertex_base] = 0.0;
                output_normals[vertex_base + 1usize] = 0.0;
                output_normals[vertex_base + 2usize] = 0.0;
                terminate!();
            }

            let mut c00 = 0.0f32;
            let mut c01 = 0.0f32;
            let mut c02 = 0.0f32;
            let mut c11 = 0.0f32;
            let mut c12 = 0.0f32;
            let mut c22 = 0.0f32;
            let mut cc00 = 0.0f32;
            let mut cc01 = 0.0f32;
            let mut cc02 = 0.0f32;
            let mut cc11 = 0.0f32;
            let mut cc12 = 0.0f32;
            let mut cc22 = 0.0f32;

            let mut zc = z0;
            while zc <= z1 {
                let mut yc = y0;
                while yc <= y1 {
                    let mut xc = x0;
                    while xc <= x1 {
                        let cell = xc + dim_x * (yc + dim_y * zc);
                        let mut slot = cell_offsets[cell] as usize;
                        let end = cell_offsets[cell + 1usize] as usize;
                        while slot < end {
                            let point_index = point_indices[slot] as usize;
                            let point_base = point_index * 3usize;
                            let qx = points[point_base];
                            let qy = points[point_base + 1usize];
                            let qz = points[point_base + 2usize];
                            let dxp = qx - px;
                            let dyp = qy - py;
                            let dzp = qz - pz;
                            let distance2 = dxp * dxp + dyp * dyp + dzp * dzp;
                            if distance2 < radius2 {
                                let ratio = distance2.sqrt() / radius;
                                let one_minus = 1.0 - ratio;
                                let square = one_minus * one_minus;
                                let weight = square * square * (4.0 * ratio + 1.0);
                                let dx = qx - centroid_x;
                                let dy = qy - centroid_y;
                                let dz = qz - centroid_z;
                                compensated_add(weight * dx * dx, &mut c00, &mut cc00);
                                compensated_add(weight * dx * dy, &mut c01, &mut cc01);
                                compensated_add(weight * dx * dz, &mut c02, &mut cc02);
                                compensated_add(weight * dy * dy, &mut c11, &mut cc11);
                                compensated_add(weight * dy * dz, &mut c12, &mut cc12);
                                compensated_add(weight * dz * dz, &mut c22, &mut cc22);
                            }
                            slot += 1usize;
                        }
                        xc += 1usize;
                    }
                    yc += 1usize;
                }
                zc += 1usize;
            }

            let inv_w = 1.0 / sum_w;
            c00 *= inv_w;
            c01 *= inv_w;
            c02 *= inv_w;
            c11 *= inv_w;
            c12 *= inv_w;
            c22 *= inv_w;

            // Scalar Jacobi state avoids backend-sensitive dynamic indexing in
            // the tiny 3x3 eigensolve. The basis is ScrollFiesta's native z-y-x.
            let mut a00 = f64::cast_from(c22);
            let mut a01 = f64::cast_from(c12);
            let mut a02 = f64::cast_from(c02);
            let mut a11 = f64::cast_from(c11);
            let mut a12 = f64::cast_from(c01);
            let mut a22 = f64::cast_from(c00);
            let mut v00 = 1.0f64;
            let mut v01 = 0.0f64;
            let mut v02 = 0.0f64;
            let mut v10 = 0.0f64;
            let mut v11 = 1.0f64;
            let mut v12 = 0.0f64;
            let mut v20 = 0.0f64;
            let mut v21 = 0.0f64;
            let mut v22 = 1.0f64;
            let mut iteration: u32 = 0;
            while iteration < 100u32 {
                if a01.abs() >= a02.abs() && a01.abs() >= a12.abs() {
                    if a01.abs() > 1.0e-15 {
                        let angle = 0.5f64 * (2.0f64 * a01).atan2(a11 - a00);
                        let c = angle.cos();
                        let s = angle.sin();
                        let old00 = a00;
                        let old01 = a01;
                        let old02 = a02;
                        let old11 = a11;
                        let old12 = a12;
                        a00 = c * c * old00 - 2.0 * s * c * old01 + s * s * old11;
                        a11 = s * s * old00 + 2.0 * s * c * old01 + c * c * old11;
                        a01 = 0.0;
                        a02 = c * old02 - s * old12;
                        a12 = s * old02 + c * old12;
                        let t00 = v00;
                        let t01 = v01;
                        let t10 = v10;
                        let t11 = v11;
                        let t20 = v20;
                        let t21 = v21;
                        v00 = c * t00 - s * t01;
                        v01 = s * t00 + c * t01;
                        v10 = c * t10 - s * t11;
                        v11 = s * t10 + c * t11;
                        v20 = c * t20 - s * t21;
                        v21 = s * t20 + c * t21;
                    }
                } else if a02.abs() >= a12.abs() {
                    if a02.abs() > 1.0e-15 {
                        let angle = 0.5f64 * (2.0f64 * a02).atan2(a22 - a00);
                        let c = angle.cos();
                        let s = angle.sin();
                        let old00 = a00;
                        let old01 = a01;
                        let old02 = a02;
                        let old12 = a12;
                        let old22 = a22;
                        a00 = c * c * old00 - 2.0 * s * c * old02 + s * s * old22;
                        a22 = s * s * old00 + 2.0 * s * c * old02 + c * c * old22;
                        a02 = 0.0;
                        a01 = c * old01 - s * old12;
                        a12 = s * old01 + c * old12;
                        let t00 = v00;
                        let t02 = v02;
                        let t10 = v10;
                        let t12 = v12;
                        let t20 = v20;
                        let t22 = v22;
                        v00 = c * t00 - s * t02;
                        v02 = s * t00 + c * t02;
                        v10 = c * t10 - s * t12;
                        v12 = s * t10 + c * t12;
                        v20 = c * t20 - s * t22;
                        v22 = s * t20 + c * t22;
                    }
                } else if a12.abs() > 1.0e-15 {
                    let angle = 0.5f64 * (2.0f64 * a12).atan2(a22 - a11);
                    let c = angle.cos();
                    let s = angle.sin();
                    let old01 = a01;
                    let old02 = a02;
                    let old11 = a11;
                    let old12 = a12;
                    let old22 = a22;
                    a11 = c * c * old11 - 2.0 * s * c * old12 + s * s * old22;
                    a22 = s * s * old11 + 2.0 * s * c * old12 + c * c * old22;
                    a12 = 0.0;
                    a01 = c * old01 - s * old02;
                    a02 = s * old01 + c * old02;
                    let t01 = v01;
                    let t02 = v02;
                    let t11 = v11;
                    let t12 = v12;
                    let t21 = v21;
                    let t22 = v22;
                    v01 = c * t01 - s * t02;
                    v02 = s * t01 + c * t02;
                    v11 = c * t11 - s * t12;
                    v12 = s * t11 + c * t12;
                    v21 = c * t21 - s * t22;
                    v22 = s * t21 + c * t22;
                }
                iteration += 1u32;
            }

            let mut nz = v00;
            let mut ny = v10;
            let mut nx = v20;
            if a11 < a00 && a11 <= a22 {
                nz = v01;
                ny = v11;
                nx = v21;
            } else if a22 < a00 && a22 < a11 {
                nz = v02;
                ny = v12;
                nx = v22;
            }
            let length2 = nx * nx + ny * ny + nz * nz;

            if length2 > 1.0e-15 {
                let inv_length = 1.0 / length2.sqrt();
                nx *= inv_length;
                ny *= inv_length;
                nz *= inv_length;
                let sign_dot = nx + ny + nz;
                if sign_dot < -1.0e-9 {
                    nx = -nx;
                    ny = -ny;
                    nz = -nz;
                } else if sign_dot < 1.0e-9 {
                    let sign = if ny.abs() > nx.abs() {
                        if nz.abs() > ny.abs() { nz } else { ny }
                    } else if nz.abs() > nx.abs() {
                        nz
                    } else {
                        nx
                    };
                    if sign < 0.0 {
                        nx = -nx;
                        ny = -ny;
                        nz = -nz;
                    }
                }
                let normal_x = f32::cast_from(nx);
                let normal_y = f32::cast_from(ny);
                let normal_z = f32::cast_from(nz);
                let signed_distance = (px - centroid_x) * normal_x
                    + (py - centroid_y) * normal_y
                    + (pz - centroid_z) * normal_z;
                output_vertices[vertex_base] = px - signed_distance * normal_x;
                output_vertices[vertex_base + 1usize] = py - signed_distance * normal_y;
                output_vertices[vertex_base + 2usize] = pz - signed_distance * normal_z;
                output_normals[vertex_base] = normal_x;
                output_normals[vertex_base + 1usize] = normal_y;
                output_normals[vertex_base + 2usize] = normal_z;
            } else {
                output_vertices[vertex_base] = px;
                output_vertices[vertex_base + 1usize] = py;
                output_vertices[vertex_base + 2usize] = pz;
                output_normals[vertex_base] = 0.0;
                output_normals[vertex_base + 1usize] = 0.0;
                output_normals[vertex_base + 2usize] = 0.0;
            }
        } else {
            output_vertices[vertex_base] = px;
            output_vertices[vertex_base + 1usize] = py;
            output_vertices[vertex_base + 2usize] = pz;
            output_normals[vertex_base] = 0.0;
            output_normals[vertex_base + 1usize] = 0.0;
            output_normals[vertex_base + 2usize] = 0.0;
        }
    }
}

pub fn project_with_runtime<R: Runtime>(
    device: &R::Device,
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
    compute_normals: bool,
) -> Result<ProjectionResult, String> {
    if vertices.is_empty() || !vertices.len().is_multiple_of(3) {
        return Err("vertices must contain one or more xyz triples".into());
    }
    if !radius.is_finite() || radius <= 0.0 {
        return Err("radius must be finite and > 0".into());
    }

    let client = R::client(device);
    let nvertices = vertices.len() / 3;
    let vertex_bytes = core::mem::size_of_val(vertices);
    let reach = (radius / arena.cell_size).ceil().max(1.0) as u32;

    let points = client.create_from_slice(f32::as_bytes(&arena.points));
    let offsets = client.create_from_slice(u32::as_bytes(&arena.cell_offsets));
    let indices = client.create_from_slice(u32::as_bytes(&arena.point_indices));
    let current = client.create_from_slice(f32::as_bytes(vertices));
    let next = client.empty(vertex_bytes);
    let normals = client.empty(vertex_bytes);
    let meta_f32_values = [
        arena.origin[0],
        arena.origin[1],
        arena.origin[2],
        arena.cell_size,
        radius,
    ];
    let meta_u32_values = [
        arena.dims[0],
        arena.dims[1],
        arena.dims[2],
        nvertices as u32,
        reach,
        u32::from(compute_normals),
    ];
    let meta_f32 = client.create_from_slice(f32::as_bytes(&meta_f32_values));
    let meta_u32 = client.create_from_slice(u32::as_bytes(&meta_u32_values));

    let cube_count_x = (nvertices as u32).div_ceil(CUBE_DIM);

    unsafe {
        mls_iteration::launch_unchecked::<R>(
            &client,
            CubeCount::Static(cube_count_x, 1, 1),
            CubeDim::new_1d(CUBE_DIM),
            ArrayArg::from_raw_parts(points, arena.points.len()),
            ArrayArg::from_raw_parts(offsets, arena.cell_offsets.len()),
            ArrayArg::from_raw_parts(indices, arena.point_indices.len()),
            ArrayArg::from_raw_parts(current, vertices.len()),
            ArrayArg::from_raw_parts(next.clone(), vertices.len()),
            ArrayArg::from_raw_parts(normals.clone(), vertices.len()),
            ArrayArg::from_raw_parts(meta_f32, meta_f32_values.len()),
            ArrayArg::from_raw_parts(meta_u32, meta_u32_values.len()),
        );
    }

    let vertex_bytes = client
        .read_one(next)
        .map_err(|error| format!("CubeCL vertex readback failed: {error}"))?;
    let normal_bytes = client
        .read_one(normals)
        .map_err(|error| format!("CubeCL normal readback failed: {error}"))?;

    Ok(ProjectionResult {
        vertices: f32::from_bytes(&vertex_bytes).to_vec(),
        normals: f32::from_bytes(&normal_bytes).to_vec(),
    })
}

#[cfg(feature = "cubecl-cpu")]
pub fn project_cpu(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, String> {
    use cubecl::cpu::{CpuDevice, CpuRuntime};
    project_with_runtime::<CpuRuntime>(&CpuDevice, arena, vertices, radius, true)
}

#[cfg(feature = "cubecl-hip")]
pub fn project_hip(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, String> {
    use cubecl::hip::{AmdDevice, HipRuntime};
    project_with_runtime::<HipRuntime>(&AmdDevice::default(), arena, vertices, radius, true)
}

#[cfg(feature = "cubecl-cuda")]
pub fn project_cuda(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, String> {
    use cubecl::cuda::{CudaDevice, CudaRuntime};
    project_with_runtime::<CudaRuntime>(&CudaDevice::default(), arena, vertices, radius, true)
}

#[cfg(feature = "cubecl-wgpu")]
pub fn project_wgpu(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, String> {
    use cubecl::wgpu::{WgpuDevice, WgpuRuntime};
    project_with_runtime::<WgpuRuntime>(&WgpuDevice::default(), arena, vertices, radius, true)
}

#[cfg(feature = "cubecl-cpu")]
pub(crate) fn project_cpu_with_normals(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
    compute_normals: bool,
) -> Result<ProjectionResult, String> {
    use cubecl::cpu::{CpuDevice, CpuRuntime};
    project_with_runtime::<CpuRuntime>(&CpuDevice, arena, vertices, radius, compute_normals)
}

#[cfg(feature = "cubecl-hip")]
pub(crate) fn project_hip_with_normals(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
    compute_normals: bool,
) -> Result<ProjectionResult, String> {
    use cubecl::hip::{AmdDevice, HipRuntime};
    project_with_runtime::<HipRuntime>(
        &AmdDevice::default(),
        arena,
        vertices,
        radius,
        compute_normals,
    )
}

#[cfg(feature = "cubecl-cuda")]
pub(crate) fn project_cuda_with_normals(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
    compute_normals: bool,
) -> Result<ProjectionResult, String> {
    use cubecl::cuda::{CudaDevice, CudaRuntime};
    project_with_runtime::<CudaRuntime>(
        &CudaDevice::default(),
        arena,
        vertices,
        radius,
        compute_normals,
    )
}

#[cfg(feature = "cubecl-wgpu")]
pub(crate) fn project_wgpu_with_normals(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
    compute_normals: bool,
) -> Result<ProjectionResult, String> {
    use cubecl::wgpu::{WgpuDevice, WgpuRuntime};
    project_with_runtime::<WgpuRuntime>(
        &WgpuDevice::default(),
        arena,
        vertices,
        radius,
        compute_normals,
    )
}
