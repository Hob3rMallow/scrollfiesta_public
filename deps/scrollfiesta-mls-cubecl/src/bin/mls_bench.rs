use herculaneum_mls_cubecl::{MlsArena, ProjectionResult, project_rust_cpu};
use std::time::{Duration, Instant};

fn synthetic_plane(side: i32, spacing: f32) -> (Vec<[f32; 3]>, Vec<f32>) {
    let mut points = Vec::new();
    let mut vertices = Vec::new();
    let half = side / 2;
    for y in -half..=half {
        for x in -half..=half {
            let xf = x as f32 * spacing;
            let yf = y as f32 * spacing;
            // Mild deterministic corrugation avoids an unrealistically exact
            // covariance while retaining a known, smooth surface.
            let z = 0.015 * (xf * 0.17).sin() * (yf * 0.13).cos();
            points.push([xf, yf, z]);
            if x % 2 == 0 && y % 2 == 0 {
                vertices.extend_from_slice(&[xf + 0.1, yf - 0.1, z + 1.5]);
            }
        }
    }
    (points, vertices)
}

fn compare(reference: &ProjectionResult, candidate: &ProjectionResult) -> (f32, f32) {
    assert_eq!(reference.vertices.len(), candidate.vertices.len());
    assert_eq!(reference.normals.len(), candidate.normals.len());
    let max_vertex = reference
        .vertices
        .iter()
        .zip(&candidate.vertices)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f32, f32::max);
    let max_normal = reference
        .normals
        .iter()
        .zip(&candidate.normals)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f32, f32::max);
    (max_vertex, max_normal)
}

fn run_backend(
    backend: &str,
    arena: &MlsArena,
    vertices: &[f32],
    radius: f32,
) -> Result<ProjectionResult, Box<dyn std::error::Error>> {
    match backend {
        "rust" | "rust-cpu" => project_rust_cpu(arena, vertices, radius)
            .map_err(|error| Box::<dyn std::error::Error>::from(std::io::Error::other(error))),
        "cube-cpu" | "cubecl-cpu" => {
            #[cfg(feature = "cubecl-cpu")]
            {
                herculaneum_mls_cubecl::project_cubecl_cpu(arena, vertices, radius).map_err(
                    |error| Box::<dyn std::error::Error>::from(std::io::Error::other(error)),
                )
            }
            #[cfg(not(feature = "cubecl-cpu"))]
            {
                Err("binary was not built with --features cubecl-cpu".into())
            }
        }
        "hip" | "cube-hip" | "cubecl-hip" => {
            #[cfg(feature = "cubecl-hip")]
            {
                herculaneum_mls_cubecl::project_cubecl_hip(arena, vertices, radius).map_err(
                    |error| Box::<dyn std::error::Error>::from(std::io::Error::other(error)),
                )
            }
            #[cfg(not(feature = "cubecl-hip"))]
            {
                Err("binary was not built with --features cubecl-hip".into())
            }
        }
        "cuda" | "cube-cuda" | "cubecl-cuda" => {
            #[cfg(feature = "cubecl-cuda")]
            {
                herculaneum_mls_cubecl::project_cubecl_cuda(arena, vertices, radius).map_err(
                    |error| Box::<dyn std::error::Error>::from(std::io::Error::other(error)),
                )
            }
            #[cfg(not(feature = "cubecl-cuda"))]
            {
                Err("binary was not built with --features cubecl-cuda".into())
            }
        }
        other => Err(format!("unknown backend {other:?}").into()),
    }
}

fn timed<T>(f: impl FnOnce() -> T) -> (T, Duration) {
    let start = Instant::now();
    let value = f();
    (value, start.elapsed())
}

fn env_f32(name: &str, default: f32) -> f32 {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let backend = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "rust-cpu".to_owned());
    let side = std::env::var("MLS_BENCH_SIDE")
        .ok()
        .and_then(|value| value.parse::<i32>().ok())
        .unwrap_or(127);
    if side < 7 || side % 2 == 0 {
        return Err("MLS_BENCH_SIDE must be an odd integer >= 7".into());
    }
    let radius = env_f32("MLS_RADIUS", 12.0);
    let max_vertex_allowed = env_f32("MLS_MAX_VERTEX_ERROR", 0.002);
    let max_normal_allowed = env_f32("MLS_MAX_NORMAL_ERROR", 0.005);

    let (points, vertices) = synthetic_plane(side, 1.0);
    let arena = MlsArena::from_points(&points, radius, None)?;
    eprintln!(
        "surface_points={} vertices={} grid={:?} radius={radius}",
        arena.npoints(),
        vertices.len() / 3,
        arena.dims
    );

    let (reference_result, reference_duration) =
        timed(|| project_rust_cpu(&arena, &vertices, radius).map_err(std::io::Error::other));
    let reference = reference_result?;

    let (cold_result, cold_duration) = timed(|| run_backend(&backend, &arena, &vertices, radius));
    let cold = cold_result?;
    let (warm_result, warm_duration) = timed(|| run_backend(&backend, &arena, &vertices, radius));
    let warm = warm_result?;

    let (cold_vertex_error, cold_normal_error) = compare(&reference, &cold);
    let (max_vertex_error, max_normal_error) = compare(&reference, &warm);

    println!("backend={backend}");
    println!(
        "reference_rust_cpu_ms={:.3}",
        reference_duration.as_secs_f64() * 1000.0
    );
    println!(
        "candidate_cold_ms={:.3}",
        cold_duration.as_secs_f64() * 1000.0
    );
    println!(
        "candidate_warm_ms={:.3}",
        warm_duration.as_secs_f64() * 1000.0
    );
    println!("cold_max_vertex_abs_error={cold_vertex_error:.9}");
    println!("cold_max_normal_abs_error={cold_normal_error:.9}");
    println!("max_vertex_abs_error={max_vertex_error:.9}");
    println!("max_normal_abs_error={max_normal_error:.9}");
    println!("first_vertex={:?}", &warm.vertices[..3]);
    println!("first_normal={:?}", &warm.normals[..3]);

    if !max_vertex_error.is_finite() || max_vertex_error > max_vertex_allowed {
        return Err(format!(
            "vertex error {max_vertex_error} exceeds MLS_MAX_VERTEX_ERROR={max_vertex_allowed}"
        )
        .into());
    }
    if !max_normal_error.is_finite() || max_normal_error > max_normal_allowed {
        return Err(format!(
            "normal error {max_normal_error} exceeds MLS_MAX_NORMAL_ERROR={max_normal_allowed}"
        )
        .into());
    }

    println!("correctness_gate=PASS");
    Ok(())
}
