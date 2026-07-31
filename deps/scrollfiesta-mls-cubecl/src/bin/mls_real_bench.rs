use herculaneum_mls_cubecl::{Backend, ProjectionResult, project_scrollfiesta_one_pass};
use std::fs;
use std::path::Path;
use std::time::Instant;

fn read_f32(path: &str) -> Result<Vec<f32>, Box<dyn std::error::Error>> {
    let bytes = fs::read(path)?;
    if bytes.len() % size_of::<f32>() != 0 {
        return Err(format!("{path} length is not a multiple of four bytes").into());
    }
    Ok(bytes
        .chunks_exact(4)
        .map(|chunk| f32::from_le_bytes(chunk.try_into().expect("four-byte chunk")))
        .collect())
}

fn write_f32(path: impl AsRef<Path>, values: &[f32]) -> std::io::Result<()> {
    let mut bytes = Vec::with_capacity(size_of_val(values));
    for value in values {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    fs::write(path, bytes)
}

fn parse_backend(value: &str) -> Result<Backend, String> {
    match value {
        "rust" | "rust-cpu" => Ok(Backend::RustCpu),
        "cube-cpu" | "cubecl-cpu" => Ok(Backend::CubeCpu),
        "hip" | "cube-hip" | "cubecl-hip" => Ok(Backend::CubeHip),
        "cuda" | "cube-cuda" | "cubecl-cuda" => Ok(Backend::CubeCuda),
        _ => Err(format!("unknown backend {value:?}")),
    }
}

fn run_passes(
    input: &[f32],
    radius: f32,
    origin: [f32; 3],
    passes: usize,
    backend: Backend,
    compute_normals: bool,
) -> Result<ProjectionResult, String> {
    let mut current = input.to_vec();
    let mut normals = vec![0.0; input.len()];
    for _ in 0..passes {
        let projected =
            project_scrollfiesta_one_pass(&current, radius, origin, backend, compute_normals)?;
        current = projected.vertices;
        normals = projected.normals;
    }
    Ok(ProjectionResult {
        vertices: current,
        normals,
    })
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    let middle = values.len() / 2;
    if values.len().is_multiple_of(2) {
        (values[middle - 1] + values[middle]) * 0.5
    } else {
        values[middle]
    }
}

fn compare(label: &str, positions: &[f32], normals: &[f32], candidate: &ProjectionResult) {
    assert_eq!(positions.len(), candidate.vertices.len());
    assert_eq!(normals.len(), candidate.normals.len());
    let mut max_position = 0.0_f64;
    let mut sum_position_squared = 0.0_f64;
    let mut max_normal_angle = 0.0_f64;
    let mut zero_normal_mismatches = 0_usize;
    let mut position_outliers = 0_usize;
    let mut normal_outliers = 0_usize;
    let mut nonfinite = 0_usize;
    for index in 0..positions.len() / 3 {
        let base = index * 3;
        if positions[base..base + 3]
            .iter()
            .chain(&normals[base..base + 3])
            .chain(&candidate.vertices[base..base + 3])
            .chain(&candidate.normals[base..base + 3])
            .any(|value| !value.is_finite())
        {
            nonfinite += 1;
            continue;
        }
        let delta = [
            f64::from(positions[base] - candidate.vertices[base]),
            f64::from(positions[base + 1] - candidate.vertices[base + 1]),
            f64::from(positions[base + 2] - candidate.vertices[base + 2]),
        ];
        let distance = delta.iter().map(|value| value * value).sum::<f64>().sqrt();
        max_position = max_position.max(distance);
        sum_position_squared += distance * distance;
        position_outliers += usize::from(distance > 1.24e-3);

        let reference = &normals[base..base + 3];
        let actual = &candidate.normals[base..base + 3];
        let reference_length = reference
            .iter()
            .map(|value| f64::from(*value).powi(2))
            .sum::<f64>()
            .sqrt();
        let actual_length = actual
            .iter()
            .map(|value| f64::from(*value).powi(2))
            .sum::<f64>()
            .sqrt();
        if reference_length <= 1.0e-12 || actual_length <= 1.0e-12 {
            zero_normal_mismatches +=
                usize::from((reference_length <= 1.0e-12) != (actual_length <= 1.0e-12));
            continue;
        }
        let dot = reference
            .iter()
            .zip(actual)
            .map(|(left, right)| f64::from(*left) * f64::from(*right))
            .sum::<f64>()
            / (reference_length * actual_length);
        let angle = dot.clamp(-1.0, 1.0).acos().to_degrees();
        max_normal_angle = max_normal_angle.max(angle);
        normal_outliers += usize::from(angle > 0.006);
    }
    let rms = (sum_position_squared / (positions.len() / 3) as f64).sqrt();
    println!("{label}_max_position_error_vox={max_position:.9}");
    println!("{label}_rms_position_error_vox={rms:.9}");
    println!("{label}_max_normal_angle_deg={max_normal_angle:.9}");
    println!("{label}_zero_normal_mismatches={zero_normal_mismatches}");
    println!("{label}_position_outliers={position_outliers}");
    println!("{label}_normal_outliers={normal_outliers}");
    println!("{label}_nonfinite_count={nonfinite}");
    println!(
        "{label}_weld_pass={}",
        nonfinite == 0 && max_position < 0.25
    );
    println!(
        "{label}_strict_parity_pass={}",
        max_position <= 1.24e-3
            && rms <= 2.2e-4
            && max_normal_angle <= 0.006
            && zero_normal_mismatches == 0
            && nonfinite == 0
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 9 && args.len() != 13 {
        return Err(format!(
            "usage: {} BACKEND INPUT.pre.positions_zyx.f32 OUTPUT_PREFIX PASSES RADIUS ORIGIN_Z ORIGIN_Y ORIGIN_X [CPU_POS CPU_NORMALS CLEAN_HIP_POS CLEAN_HIP_NORMALS]",
            args[0]
        )
        .into());
    }
    let backend = parse_backend(&args[1])?;
    let input = read_f32(&args[2])?;
    if input.is_empty() || !input.len().is_multiple_of(3) {
        return Err("input must contain one or more zyx triples".into());
    }
    let passes: usize = args[4].parse()?;
    let radius: f32 = args[5].parse()?;
    let origin = [args[6].parse()?, args[7].parse()?, args[8].parse()?];
    if passes == 0 {
        return Err("PASSES must be greater than zero".into());
    }
    let repetitions = std::env::var("MLS_BENCH_REPS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(3)
        .max(1);
    let compute_normals = std::env::var("MLS_COMPUTE_NORMALS").map_or(true, |value| value != "0");

    let _warm = run_passes(&input, radius, origin, passes, backend, compute_normals)?;
    let mut timings = Vec::with_capacity(repetitions);
    let mut result = None;
    for repetition in 0..repetitions {
        let start = Instant::now();
        let projected = run_passes(&input, radius, origin, passes, backend, compute_normals)?;
        let elapsed = start.elapsed().as_secs_f64() * 1000.0;
        println!("cubecl_rep_{repetition}_ms={elapsed:.3}");
        timings.push(elapsed);
        result = Some(projected);
    }
    let result = result.expect("at least one repetition");
    println!("backend={:?}", backend);
    println!("count={}", input.len() / 3);
    println!("passes={passes}");
    println!("compute_normals={compute_normals}");
    println!("cubecl_median_ms={:.3}", median(timings));

    let positions_path = format!("{}.cubecl.positions_zyx.f32", args[3]);
    let normals_path = format!("{}.cubecl.normals_zyx.f32", args[3]);
    write_f32(&positions_path, &result.vertices)?;
    write_f32(&normals_path, &result.normals)?;
    println!("cubecl_positions={positions_path}");
    println!("cubecl_normals={normals_path}");

    if args.len() == 13 {
        compare("cpu", &read_f32(&args[9])?, &read_f32(&args[10])?, &result);
        compare(
            "clean_hip",
            &read_f32(&args[11])?,
            &read_f32(&args[12])?,
            &result,
        );
    }
    Ok(())
}
