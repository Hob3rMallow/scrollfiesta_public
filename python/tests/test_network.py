"""Live S3 tests against the real PHerc0139 surface-prediction zarr (marked
network; deselect with -m 'not network')."""

import os

import numpy as np
import pytest

from scrollunwrap.roi_finder import auto_pick_rois
from scrollunwrap.zarr_source import (aws_credentials_available, open_volume,
                                      read_block, read_padded_cube)

pytestmark = pytest.mark.network

ZARR = ("s3://vesuvius-challenge-open-data/PHerc0139/representations/predictions/"
        "surfaces/20250728140407-surface-20260413222639-surface-m7-L0-th0.2.zarr")


def test_open_levels():
    v0 = open_volume(ZARR, 0)
    assert v0.dtype == np.uint8 and len(v0.shape) == 3
    assert v0.shape == (20974, 6621, 6621)        # confirmed level-0 shape
    v5 = open_volume(ZARR, 5)
    assert all(s > 0 for s in v5.shape) and v5.scale == 32


def test_auto_pick_real():
    rois = auto_pick_rois(ZARR, n=2, scan_level=5, roi_cubes=(2, 2, 2), target_level=0)
    assert len(rois) >= 1
    assert rois[0].density > 0
    v0 = open_volume(ZARR, 0)
    z0, z1, y0, y1, x0, x1 = rois[0].bbox_l0
    assert 0 <= z0 < z1 <= v0.shape[0]
    assert 0 <= y0 < y1 <= v0.shape[1]
    assert 0 <= x0 < x1 <= v0.shape[2]


def test_read_padded_cube_real():
    rois = auto_pick_rois(ZARR, n=1, scan_level=5, roi_cubes=(1, 1, 1), target_level=0)
    assert rois
    z0, _, y0, _, x0, _ = rois[0].bbox_l0
    v0 = open_volume(ZARR, 0)
    buf = read_padded_cube(v0, z0, y0, x0, size=128, halo=13)
    assert buf.shape == (154, 154, 154) and buf.dtype == np.uint8
    assert int((buf != 0).sum()) > 0     # a dense ROI must contain surface


def test_private_s3_signed_access():
    """Signed access to a PRIVATE OME-Zarr. Skips unless AWS creds (e.g. STS)
    are in the environment AND SCROLLUNWRAP_PRIVATE_ZARR points at a private
    store — so this never hardcodes credentials or a private URI."""
    if not aws_credentials_available():
        pytest.skip("no AWS credentials in environment")
    uri = os.environ.get("SCROLLUNWRAP_PRIVATE_ZARR")
    if not uri:
        pytest.skip("set SCROLLUNWRAP_PRIVATE_ZARR to a private OME-Zarr URI")
    v = open_volume(uri, 0)              # anon=None -> auto -> signed
    assert v.dtype == np.uint8 and len(v.shape) == 3
    s = min(64, *v.shape)
    blk = read_block(v, 0, s, 0, s, 0, s)   # exercises a signed data fetch
    assert blk.size > 0
