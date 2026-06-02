"""The most important correctness test: scrollfiesta (z,y,x) at level L must
become (x,y,z) in level-0 voxels, since vc_obj2tifxyz_legacy does no axis swap
and maps OBJ col1->x.tif, col2->y.tif, col3->z.tif."""

import numpy as np

from scrollunwrap.meshprep import to_level0_xyz, write_obj_vf


def test_to_level0_xyz_reverse_and_scale():
    V = np.array([[3.0, 5.0, 7.0], [1.0, 2.0, 4.0]])
    out = to_level0_xyz(V, level=2)            # reverse triplet, then x4
    assert np.allclose(out, [[28, 20, 12], [16, 8, 4]])


def test_to_level0_xyz_level0_is_just_reverse():
    V = np.array([[3.0, 5.0, 7.0]])
    assert np.allclose(to_level0_xyz(V, 0), [[7, 5, 3]])


def test_write_obj_vf_roundtrip(tmp_path):
    V = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    F = np.array([[0, 1, 2]])
    p = tmp_path / "m.obj"
    write_obj_vf(p, V, F)
    Vr, Fr = [], []
    for line in open(p):
        if line.startswith("v "):
            Vr.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("f "):
            Fr.append([int(t.split("/")[0]) for t in line.split()[1:]])
    assert np.allclose(Vr, V)
    assert Fr == [[1, 2, 3]]                    # OBJ is 1-indexed


def test_full_chain_zyx_voxel_to_xyz_obj(tmp_path):
    """A scrollfiesta vertex (z=100,y=50,x=10) at level 1 -> OBJ 'v 20 100 200'
    so x.tif=20, y.tif=100, z.tif=200 in level-0 scroll voxels."""
    V = np.array([[100.0, 50.0, 10.0]])
    out = to_level0_xyz(V, level=1)
    p = tmp_path / "v.obj"
    write_obj_vf(p, out, np.zeros((0, 3), dtype=int))
    line = next(ln for ln in open(p) if ln.startswith("v "))
    x, y, z = (float(v) for v in line.split()[1:4])
    assert (x, y, z) == (20.0, 100.0, 200.0)
