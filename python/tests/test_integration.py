"""Integration tests that drive the native binaries (marked requires_binary)."""

import subprocess
from pathlib import Path

import numpy as np
import pytest
import tifffile
from _helpers import grid_disk, make_local_zarr

from scrollunwrap.cube_planner import plan_cubes
from scrollunwrap.flatten import run_flatboi
from scrollunwrap.mesher import run_mesher
from scrollunwrap.meshprep import write_obj_vf
from scrollunwrap.tifxyz import run_obj2tifxyz
from scrollunwrap.zarr_source import open_volume

pytestmark = pytest.mark.requires_binary


def _slab_cube(p, halo):
    """p^3 padded buffer with a thin recto sheet in the owned interior
    (owned-local z[62:65], y/x[16:112]) — geometry that survives the pipeline."""
    v = np.zeros((p, p, p), np.uint8)
    v[halo + 62:halo + 65, halo + 16:halo + 112, halo + 16:halo + 112] = 255
    return v


def test_stdin_raw_world_offset(cube_mesh_bin, tmp_path):
    """cube_mesh --stdin-raw ingests an in-memory padded cube and applies the
    world offset from (oz,oy,ox). Owned origin (128,128,128) + slab at owned-z
    ~63 -> vertex Z ~191."""
    halo = 13
    p = 128 + 2 * halo
    buf = _slab_cube(p, halo)
    dump = tmp_path / "dump"
    cmd = [cube_mesh_bin, "--stdin-raw", str(p), "128", "128", "128",
           "--halo", str(halo), "--dump-obj", str(dump)]
    proc = subprocess.run(cmd, input=buf.tobytes(), capture_output=True)
    assert proc.returncode == 0, proc.stderr.decode()[-800:]
    rs = dump / "z00128_y00128_x00128" / "z00128_y00128_x00128_raw_snap" / \
        "z00128_y00128_x00128_raw_snap_all.obj"
    assert rs.exists(), "raw_snap not produced"
    zs = [float(line.split()[1]) for line in open(rs) if line.startswith("v ")]
    assert len(zs) > 0
    # vertex col1 is Z in world voxels = 128 (origin) + owned-local slab z (~63)
    assert 185.0 <= min(zs) and max(zs) <= 198.0, (min(zs), max(zs))


def test_stream_mesh_and_weld(cube_mesh_bin, grid_weld_bin, tmp_path):
    """Stream a local zarr cube directly into the mesher (no TIFF) and weld."""
    uri, _ = make_local_zarr(tmp_path / "v.zarr", shape=(128, 128, 128),
                             chunks=(64, 64, 64),
                             blob=(slice(62, 65), slice(16, 112), slice(16, 112)))
    vol = open_volume(uri, 0)
    cubes = plan_cubes((0, 128, 0, 128, 0, 128), level=0)
    assert len(cubes) == 1 and cubes[0].cube_id == "z00000_y00000_x00000"
    res = run_mesher(vol, cubes, tmp_path / "mesh",
                     cube_mesh_bin=cube_mesh_bin, grid_weld_bin=grid_weld_bin,
                     halo=13, max_concurrent=1)
    assert res.n_obj >= 1
    assert res.welded_obj.exists()
    nverts = sum(1 for ln in open(res.welded_obj) if ln.startswith("v "))
    assert nverts > 0


def test_flatten_and_tifxyz_coordinate_fidelity(flatboi_bin, obj2tifxyz_bin, tmp_path):
    """A planar disk at z=100, x,y in [0,19] -> flatboi -> tifxyz. z.tif must be
    ~100 everywhere valid; x.tif must span ~[0,19] (proves coords survive the
    flatten + legacy-converter chain unchanged)."""
    m = grid_disk(20, translate=(0.0, 0.0, 100.0))   # x,y in [0,19], z=100
    prepped = tmp_path / "prepped.obj"
    write_obj_vf(prepped, np.asarray(m.vertices), np.asarray(m.faces))
    flat = run_flatboi(prepped, binary=flatboi_bin, iters=20)
    assert flat.exists()
    out, meta, _ = run_obj2tifxyz(flat, tmp_path / "tifxyz", binary=obj2tifxyz_bin,
                                  step_size=2)
    assert meta["format"] == "tifxyz"
    z = tifffile.imread(out / "z.tif").astype(np.float64)
    x = tifffile.imread(out / "x.tif").astype(np.float64)
    valid = z > 0
    assert valid.sum() > 0
    assert np.allclose(z[valid], 100.0, atol=1.0)        # constant-z plane preserved
    assert -1.0 <= x[valid].min() and x[valid].max() <= 21.0
