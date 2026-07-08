import numpy as np
import trimesh
from _helpers import grid_disk

from scrollunwrap import geometry_report as gr


def test_analyze_uv_flip_detection():
    # tri0 CCW (area>0), tri1 CW (area<0) -> exactly one orientation flip
    face_uv = np.array([
        [[0, 0], [1, 0], [0, 1]],
        [[0, 0], [0, 1], [1, 0]],
    ], dtype=float)
    r = gr.analyze_uv(face_uv)
    assert r["n_faces"] == 2
    assert r["flipped_uv_triangles"] == 1


def test_analyze_uv_all_consistent():
    face_uv = np.array([
        [[0, 0], [1, 0], [0, 1]],
        [[1, 1], [0, 1], [1, 0]],
    ], dtype=float)
    r = gr.analyze_uv(face_uv)
    assert r["flipped_uv_triangles"] == 0


def test_load_obj_uv(tmp_path):
    p = tmp_path / "f.obj"
    p.write_text("v 0 0 0\nv 1 0 0\nv 0 1 0\n"
                 "vt 0 0\nvt 1 0\nvt 0 1\n"
                 "f 1/1 2/2 3/3\n")
    V, F, fuv = gr.load_obj_uv(p)
    assert V.shape == (3, 3) and F.shape == (1, 3) and fuv.shape == (1, 3, 2)
    assert np.allclose(fuv[0], [[0, 0], [1, 0], [0, 1]])


def test_analyze_welded_counts_components():
    big = grid_disk(6)
    small = grid_disk(3, translate=(50, 0, 0))
    combo = trimesh.util.concatenate([big, small])
    a = gr.analyze_welded(combo)
    assert a["n_components"] == 2
    assert a["boundary_loops"] >= 2
