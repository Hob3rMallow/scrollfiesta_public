use herculaneum_mls_cubecl::{
    MLS_cubecl_arena_create_from_points, MLS_cubecl_arena_destroy, MLS_cubecl_arena_set_backend,
    MLS_cubecl_last_error, MLS_cubecl_project_with_arena, MLS_project_verts,
};
use std::ffi::CStr;

#[test]
fn explicit_arena_api_projects_a_plane() {
    let mut points = Vec::new();
    for y in -6..=6 {
        for x in -6..=6 {
            points.extend_from_slice(&[x as f32, y as f32, 0.0]);
        }
    }
    let origin = [-10.0_f32, -10.0, -10.0];
    let arena = unsafe {
        MLS_cubecl_arena_create_from_points(points.as_ptr(), points.len() / 3, 4.0, origin.as_ptr())
    };
    assert!(!arena.is_null(), "{}", unsafe {
        CStr::from_ptr(MLS_cubecl_last_error()).to_string_lossy()
    });
    assert_eq!(unsafe { MLS_cubecl_arena_set_backend(arena, 1) }, 0);

    let input = [0.25_f32, -0.5, 1.75];
    let mut output = [f32::NAN; 3];
    let mut normals = [f32::NAN; 3];
    unsafe {
        MLS_cubecl_project_with_arena(
            arena,
            input.as_ptr(),
            1,
            4.0,
            output.as_mut_ptr(),
            normals.as_mut_ptr(),
        );
    }
    let error = unsafe { CStr::from_ptr(MLS_cubecl_last_error()).to_string_lossy() };
    assert_eq!(error, "no error");
    assert!(output[2].abs() < 1.0e-5, "output={output:?}");
    assert!(normals[2] > 0.9999, "normals={normals:?}");
    unsafe { MLS_cubecl_arena_destroy(arena) };
}

#[test]
fn legacy_symbol_accepts_foreign_scratch_and_rebuilds_current_support() {
    let mut input = Vec::new();
    for y in -3..=3 {
        for x in -3..=3 {
            let z = if x == 0 && y == 0 { 0.75 } else { 0.0 };
            input.extend_from_slice(&[z, y as f32, x as f32]);
        }
    }
    let center = (3 * 7 + 3) * 3;
    let mut output = vec![f32::NAN; input.len()];
    let mut normals = vec![f32::NAN; input.len()];
    let origin = [16128.0_f32, 2560.0, 7680.0];
    unsafe {
        MLS_project_verts(
            std::ptr::dangling_mut(),
            input.as_ptr(),
            input.len() / 3,
            4.0,
            origin.as_ptr(),
            output.as_mut_ptr(),
            normals.as_mut_ptr(),
        );
    }
    let error = unsafe { CStr::from_ptr(MLS_cubecl_last_error()).to_string_lossy() };
    assert_eq!(error, "no error");
    assert!(
        output[center].abs() < input[center],
        "output={:?}",
        &output[center..center + 3]
    );
    assert!(
        normals[center].abs() > 0.9,
        "normal={:?}",
        &normals[center..center + 3]
    );

    let mut in_place = input.clone();
    let mut in_place_normals = vec![f32::NAN; input.len()];
    unsafe {
        MLS_project_verts(
            std::ptr::dangling_mut(),
            in_place.as_ptr(),
            in_place.len() / 3,
            4.0,
            origin.as_ptr(),
            in_place.as_mut_ptr(),
            in_place_normals.as_mut_ptr(),
        );
    }
    assert_eq!(in_place, output);
    assert_eq!(in_place_normals, normals);
}
