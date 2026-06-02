from _helpers import make_local_zarr

from scrollunwrap.roi_finder import auto_pick_rois


def test_auto_pick_covers_blob(tmp_path):
    uri, _ = make_local_zarr(tmp_path / "v.zarr", shape=(40, 40, 40),
                             blob=(slice(8, 24), slice(8, 24), slice(8, 24)))
    rois = auto_pick_rois(uri, n=1, scan_level=0, roi_cubes=(1, 1, 1),
                          target_level=0, cube=16, min_separation_cells=1)
    assert len(rois) == 1
    z0, z1, y0, y1, x0, x1 = rois[0].bbox_l0
    # the chosen 16-voxel block must overlap the blob [8:24]
    assert z0 < 24 and z1 > 8 and y0 < 24 and y1 > 8 and x0 < 24 and x1 > 8
    assert rois[0].density > 0


def test_auto_pick_empty_volume_returns_nothing(tmp_path):
    uri, _ = make_local_zarr(tmp_path / "e.zarr", shape=(40, 40, 40), blob=None)
    rois = auto_pick_rois(uri, n=3, scan_level=0, roi_cubes=(1, 1, 1),
                          target_level=0, cube=16)
    assert rois == []


def test_auto_pick_distinct_and_capped(tmp_path):
    uri, _ = make_local_zarr(tmp_path / "v.zarr", shape=(48, 48, 48),
                             blob=(slice(0, 48), slice(0, 48), slice(0, 48)))
    rois = auto_pick_rois(uri, n=3, scan_level=0, roi_cubes=(1, 1, 1),
                          target_level=0, cube=16, min_separation_cells=1,
                          density_max=1.0)   # whole-volume blob: disable the cap
    assert 1 <= len(rois) <= 3
    assert len({r.bbox_l0 for r in rois}) == len(rois)   # distinct


def test_auto_pick_density_band_excludes_blob(tmp_path):
    # a fully-solid volume is a "blob"; the default band must reject it
    uri, _ = make_local_zarr(tmp_path / "blob.zarr", shape=(48, 48, 48),
                             blob=(slice(0, 48), slice(0, 48), slice(0, 48)))
    assert auto_pick_rois(uri, n=3, scan_level=0, roi_cubes=(1, 1, 1),
                          target_level=0, cube=16) == []
