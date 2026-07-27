use crate::grid::{MlsArena, linear_cell, point_cell_coords};

const EPS_WEIGHT: f32 = 1.0e-12;
const EPS_EIGEN: f32 = 1.0e-12;

#[derive(Debug, Clone)]
pub struct ProjectionResult {
    pub vertices: Vec<f32>,
    pub normals: Vec<f32>,
}

#[inline]
fn wendland_c2(distance_squared: f32, radius: f32) -> f32 {
    let q = distance_squared.sqrt() / radius;
    if q >= 1.0 {
        0.0
    } else {
        let one_minus_q = 1.0 - q;
        let square = one_minus_q * one_minus_q;
        square * square * (4.0 * q + 1.0)
    }
}

#[inline]
fn orient_normal(mut n: [f32; 3]) -> [f32; 3] {
    let length2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    if length2 <= EPS_EIGEN || !length2.is_finite() {
        return [0.0, 0.0, 0.0];
    }
    let inv = length2.sqrt().recip();
    n[0] *= inv;
    n[1] *= inv;
    n[2] *= inv;
    let sign_dot = n[0] + n[1] + n[2];
    if sign_dot < -1.0e-9 {
        n = [-n[0], -n[1], -n[2]];
    } else if sign_dot < 1.0e-9 {
        let largest = if n[1].abs() > n[0].abs() { 1 } else { 0 };
        let largest = if n[2].abs() > n[largest].abs() {
            2
        } else {
            largest
        };
        if n[largest] < 0.0 {
            n = [-n[0], -n[1], -n[2]];
        }
    }
    n
}

pub fn smallest_eigenvector_symmetric(covariance: [[f32; 3]; 3]) -> [f32; 3] {
    // Solve in ScrollFiesta's native z-y-x basis, then map the eigenvector
    // back to xyz. This keeps tied-eigenspace behaviour aligned with upstream.
    let mut a = [
        [
            f64::from(covariance[2][2]),
            f64::from(covariance[2][1]),
            f64::from(covariance[2][0]),
        ],
        [
            f64::from(covariance[1][2]),
            f64::from(covariance[1][1]),
            f64::from(covariance[1][0]),
        ],
        [
            f64::from(covariance[0][2]),
            f64::from(covariance[0][1]),
            f64::from(covariance[0][0]),
        ],
    ];
    let mut v = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..100 {
        let mut p = 0;
        let mut q = 1;
        let mut max_off = a[0][1].abs();
        if a[0][2].abs() > max_off {
            max_off = a[0][2].abs();
            q = 2;
        }
        if a[1][2].abs() > max_off {
            max_off = a[1][2].abs();
            p = 1;
            q = 2;
        }
        if max_off < 1.0e-15 {
            break;
        }
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        let angle = 0.5 * (2.0 * apq).atan2(aqq - app);
        let c = angle.cos();
        let s = angle.sin();
        let mut next = a;
        next[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        next[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        next[p][q] = 0.0;
        next[q][p] = 0.0;
        for row in 0..3 {
            if row != p && row != q {
                next[row][p] = c * a[row][p] - s * a[row][q];
                next[p][row] = next[row][p];
                next[row][q] = s * a[row][p] + c * a[row][q];
                next[q][row] = next[row][q];
            }
        }
        a = next;
        for row in &mut v {
            let vp = row[p];
            let vq = row[q];
            row[p] = c * vp - s * vq;
            row[q] = s * vp + c * vq;
        }
    }
    let mut smallest = 0;
    if a[1][1] < a[smallest][smallest] {
        smallest = 1;
    }
    if a[2][2] < a[smallest][smallest] {
        smallest = 2;
    }
    orient_normal([
        v[2][smallest] as f32,
        v[1][smallest] as f32,
        v[0][smallest] as f32,
    ])
}

#[inline]
fn compensated_add(value: f32, sum: &mut f32, compensation: &mut f32) {
    let corrected = value - *compensation;
    let next = *sum + corrected;
    *compensation = (next - *sum) - corrected;
    *sum = next;
}

#[inline]
fn each_neighbour<F>(arena: &MlsArena, p: [f32; 3], radius: f32, mut visit: F)
where
    F: FnMut([f32; 3], f32),
{
    let c = point_cell_coords(p, arena.origin, arena.cell_size, arena.dims);
    let reach = (radius / arena.cell_size).ceil().max(1.0) as u32;
    let x0 = c[0].saturating_sub(reach);
    let y0 = c[1].saturating_sub(reach);
    let z0 = c[2].saturating_sub(reach);
    let x1 = c[0].saturating_add(reach).min(arena.dims[0] - 1);
    let y1 = c[1].saturating_add(reach).min(arena.dims[1] - 1);
    let z1 = c[2].saturating_add(reach).min(arena.dims[2] - 1);
    let r2 = radius * radius;

    for z in z0..=z1 {
        for y in y0..=y1 {
            for x in x0..=x1 {
                let cell = linear_cell(x, y, z, arena.dims);
                let begin = arena.cell_offsets[cell] as usize;
                let end = arena.cell_offsets[cell + 1] as usize;
                for slot in begin..end {
                    let point_index = arena.point_indices[slot] as usize;
                    let base = point_index * 3;
                    let q = [
                        arena.points[base],
                        arena.points[base + 1],
                        arena.points[base + 2],
                    ];
                    let dx = q[0] - p[0];
                    let dy = q[1] - p[1];
                    let dz = q[2] - p[2];
                    let d2 = dx * dx + dy * dy + dz * dz;
                    if d2 < r2 {
                        let weight = wendland_c2(d2, radius);
                        if weight > 0.0 {
                            visit(q, weight);
                        }
                    }
                }
            }
        }
    }
}

fn project_vertex_once(arena: &MlsArena, p: [f32; 3], radius: f32) -> ([f32; 3], [f32; 3]) {
    let mut sum_w = 0.0_f32;
    let mut centroid = [0.0_f32; 3];
    let mut compensation_w = 0.0_f32;
    let mut compensation = [0.0_f32; 3];
    let mut count = 0_u32;
    each_neighbour(arena, p, radius, |q, w| {
        compensated_add(w, &mut sum_w, &mut compensation_w);
        for axis in 0..3 {
            compensated_add(
                w * (q[axis] - p[axis]),
                &mut centroid[axis],
                &mut compensation[axis],
            );
        }
        count += 1;
    });

    if count < 4 || sum_w <= EPS_WEIGHT || !sum_w.is_finite() {
        return (p, [0.0, 0.0, 0.0]);
    }
    for axis in 0..3 {
        centroid[axis] = p[axis] + centroid[axis] / sum_w;
    }

    let mut covariance = [[0.0_f32; 3]; 3];
    let mut covariance_compensation = [[0.0_f32; 3]; 3];
    each_neighbour(arena, p, radius, |q, w| {
        let d = [q[0] - centroid[0], q[1] - centroid[1], q[2] - centroid[2]];
        for row in 0..3 {
            for column in row..3 {
                compensated_add(
                    w * d[row] * d[column],
                    &mut covariance[row][column],
                    &mut covariance_compensation[row][column],
                );
            }
        }
    });
    covariance[1][0] = covariance[0][1];
    covariance[2][0] = covariance[0][2];
    covariance[2][1] = covariance[1][2];

    let inv_w = sum_w.recip();
    for row in &mut covariance {
        for value in row {
            *value *= inv_w;
        }
    }

    let normal = smallest_eigenvector_symmetric(covariance);
    if normal == [0.0, 0.0, 0.0] {
        return (p, normal);
    }
    let signed_distance = (p[0] - centroid[0]) * normal[0]
        + (p[1] - centroid[1]) * normal[1]
        + (p[2] - centroid[2]) * normal[2];
    let projected = [
        p[0] - signed_distance * normal[0],
        p[1] - signed_distance * normal[1],
        p[2] - signed_distance * normal[2],
    ];
    (projected, normal)
}

pub fn project(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, &'static str> {
    if !vertices.len().is_multiple_of(3) {
        return Err("vertices must be flat xyz triples");
    }
    if !radius.is_finite() || radius <= 0.0 {
        return Err("radius must be finite and > 0");
    }
    if !vertices.iter().all(|x| x.is_finite()) {
        return Err("vertices contain a non-finite value");
    }

    project_with_normals(arena, vertices, radius, true)
}

fn project_vertex_centroid(arena: &MlsArena, p: [f32; 3], radius: f32) -> [f32; 3] {
    let mut sum_w = 0.0_f32;
    let mut centroid = [0.0_f32; 3];
    let mut compensation_w = 0.0_f32;
    let mut compensation = [0.0_f32; 3];
    let mut count = 0_u32;
    each_neighbour(arena, p, radius, |support, weight| {
        compensated_add(weight, &mut sum_w, &mut compensation_w);
        for axis in 0..3 {
            compensated_add(
                weight * (support[axis] - p[axis]),
                &mut centroid[axis],
                &mut compensation[axis],
            );
        }
        count += 1;
    });
    if count < 4 || sum_w <= EPS_WEIGHT || !sum_w.is_finite() {
        return p;
    }
    for axis in 0..3 {
        centroid[axis] = p[axis] + centroid[axis] / sum_w;
    }
    centroid
}

pub(crate) fn project_with_normals(
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
    compute_normals: bool,
) -> Result<ProjectionResult, &'static str> {
    if vertices.is_empty() || !vertices.len().is_multiple_of(3) {
        return Err("vertices must contain one or more flat xyz triples");
    }
    if !radius.is_finite() || radius <= 0.0 {
        return Err("radius must be finite and > 0");
    }
    if !vertices.iter().all(|x| x.is_finite()) {
        return Err("vertices contain a non-finite value");
    }

    let mut next = vec![0.0_f32; vertices.len()];
    let mut normals = vec![0.0_f32; vertices.len()];

    #[cfg(feature = "rust-cpu")]
    {
        use rayon::prelude::*;
        next.par_chunks_mut(3)
            .zip(normals.par_chunks_mut(3))
            .enumerate()
            .for_each(|(index, (out, normal_out))| {
                let base = index * 3;
                let p = [vertices[base], vertices[base + 1], vertices[base + 2]];
                let (q, n) = project_vertex_once(arena, p, radius);
                if compute_normals {
                    out.copy_from_slice(&q);
                    normal_out.copy_from_slice(&n);
                } else {
                    out.copy_from_slice(&project_vertex_centroid(arena, p, radius));
                }
            });
    }
    #[cfg(not(feature = "rust-cpu"))]
    {
        for index in 0..(vertices.len() / 3) {
            let base = index * 3;
            let p = [vertices[base], vertices[base + 1], vertices[base + 2]];
            if compute_normals {
                let (q, n) = project_vertex_once(arena, p, radius);
                next[base..base + 3].copy_from_slice(&q);
                normals[base..base + 3].copy_from_slice(&n);
            } else {
                next[base..base + 3].copy_from_slice(&project_vertex_centroid(arena, p, radius));
            }
        }
    }

    Ok(ProjectionResult {
        vertices: next,
        normals,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::MlsArena;

    #[test]
    fn diagonal_covariance_returns_smallest_axis() {
        let n =
            smallest_eigenvector_symmetric([[4.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.25]]);
        assert!(n[2] > 0.9999, "{n:?}");
    }

    #[test]
    fn plane_projection_moves_vertices_to_plane() {
        let mut points = Vec::new();
        for y in -6..=6 {
            for x in -6..=6 {
                points.push([x as f32, y as f32, 0.0]);
            }
        }
        let arena = MlsArena::from_points(&points, 4.0, None).unwrap();
        let result = project(&arena, &[0.25, -0.5, 1.75], 4.0).unwrap();
        assert!(result.vertices[2].abs() < 1.0e-5, "{:?}", result.vertices);
        assert!(result.normals[2] > 0.9999, "{:?}", result.normals);
    }

    #[test]
    fn public_project_is_exactly_one_fixed_support_pass() {
        let mut points = Vec::new();
        for y in -4..=4 {
            for x in -4..=4 {
                let z = 0.03 * (x * x) as f32 + 0.01 * y as f32;
                points.push([x as f32, y as f32, z]);
            }
        }
        let arena = MlsArena::from_points(&points, 4.0, None).unwrap();
        let input = [0.4, -0.2, 1.2];
        let expected = project_vertex_once(&arena, input, 4.0);
        let actual = project(&arena, &input, 4.0).unwrap();
        assert_eq!(actual.vertices, expected.0);
        assert_eq!(actual.normals, expected.1);
    }
}
