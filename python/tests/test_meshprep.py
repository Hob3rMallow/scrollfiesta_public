import numpy as np
import pytest
import trimesh
from _helpers import grid_disk

from scrollunwrap import meshprep as mp
from scrollunwrap.meshprep import PrepError


def test_disk_classification():
    d = grid_disk(6)
    assert mp.euler_genus(d) == (1, 0, 1)
    assert mp.is_disk(d)


def test_sphere_not_disk():
    s = trimesh.creation.icosphere(subdivisions=1)
    chi, g, b = mp.euler_genus(s)
    assert (chi, g, b) == (2, 0, 0)
    assert not mp.is_disk(s)


def test_holed_then_filled():
    holed = grid_disk(7, hole=True)
    assert mp.euler_genus(holed)[2] == 2          # two boundary loops
    assert not mp.is_disk(holed)
    filled = mp.fill_small_interior_holes(mp.make_manifold(holed), max_edges=30)
    assert mp.is_disk(filled)


def test_select_components_filters_and_ranks():
    big = grid_disk(8)
    small = grid_disk(4, translate=(100, 0, 0))
    combo = trimesh.util.concatenate([big, small])
    comps = mp.select_components(combo, min_faces=1)
    assert len(comps) == 2
    assert len(comps[0].faces) >= len(comps[1].faces)
    assert len(mp.select_components(combo, min_faces=len(big.faces))) == 1
    assert len(mp.select_components(combo, mode="largest")) == 1


def test_prep_component_disk_and_level_scaling():
    pp = mp.prep_component(grid_disk(6), level=2, fill_max_edges=30)
    assert pp.info["genus"] == 0 and pp.info["boundary_loops"] == 1
    assert pp.mesh.vertices.shape[1] == 3


def test_prep_component_rejects_sphere():
    with pytest.raises(PrepError):
        mp.prep_component(trimesh.creation.icosphere(subdivisions=1))
