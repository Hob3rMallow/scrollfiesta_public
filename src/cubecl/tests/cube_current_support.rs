#![cfg(feature = "cubecl-cpu")]

use herculaneum_mls_cubecl::{MlsArena, project_cubecl_cpu, project_rust_cpu};

#[test]
fn cube_kernel_matches_rust_when_queries_are_the_support_cloud() {
    let mut points = Vec::new();
    let mut vertices = Vec::new();
    for y in -15..=15 {
        for x in -15..=15 {
            let xf = x as f32 * 0.7;
            let yf = y as f32 * 0.7;
            let z = 0.2 * (xf * 0.31).sin() + 0.15 * (yf * 0.23).cos();
            points.push([xf, yf, z]);
            vertices.extend_from_slice(&[xf, yf, z]);
        }
    }
    let arena = MlsArena::from_points(&points, 4.0, None).unwrap();
    let reference = project_rust_cpu(&arena, &vertices, 4.0).unwrap();
    let candidate = project_cubecl_cpu(&arena, &vertices, 4.0).unwrap();
    let max_vertex = reference
        .vertices
        .iter()
        .zip(candidate.vertices)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f32, f32::max);
    let max_normal = reference
        .normals
        .iter()
        .zip(candidate.normals)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f32, f32::max);
    assert!(max_vertex < 1.0e-5, "max vertex error {max_vertex}");
    assert!(max_normal < 1.0e-5, "max normal error {max_normal}");
}
