import dataclasses

import numpy as np
from _helpers import make_local_zarr

from scrollunwrap.zarr_source import (apply_threshold, open_volume,
                                      read_block, read_padded_cube)


def test_open_volume_metadata(tmp_path):
    uri, _ = make_local_zarr(tmp_path / "v.zarr")
    vol = open_volume(uri, 0)
    assert vol.shape == (40, 40, 40)
    assert vol.chunks == (16, 16, 16)
    assert vol.dtype == np.uint8
    assert vol.level == 0 and vol.scale == 1


def test_read_padded_cube_owned_and_edge_pad(tmp_path):
    uri, vol = make_local_zarr(tmp_path / "v.zarr")  # blob at [8:24]^3, value 200
    zv = open_volume(uri, 0)
    buf = read_padded_cube(zv, 8, 8, 8, size=16, halo=4)
    assert buf.shape == (24, 24, 24) and buf.dtype == np.uint8
    # owned region [4:20] maps to world [8:24] = the whole blob -> 16^3 nonzero
    assert int((buf[4:20, 4:20, 4:20] != 0).sum()) == 16 ** 3
    # low halo extends below 0 -> stays zero
    buf0 = read_padded_cube(zv, 0, 0, 0, size=16, halo=4)
    assert buf0[0:4, :, :].sum() == 0


class _Spy:
    """Records slicing so we can prove we never read the whole array."""

    def __init__(self, arr):
        self.arr = arr
        self.keys = []

    def __getitem__(self, k):
        self.keys.append(k)
        return self.arr[k]


def test_streaming_never_reads_full_array(tmp_path):
    uri, _ = make_local_zarr(tmp_path / "v.zarr")
    vol = open_volume(uri, 0)
    spy = _Spy(vol.array)
    vol = dataclasses.replace(vol, array=spy)
    read_padded_cube(vol, 8, 8, 8, size=16, halo=4)
    assert len(spy.keys) == 1
    key = spy.keys[0]
    # each axis slice must be a strict sub-range of the 40-voxel extent
    for sl, dim in zip(key, vol.shape):
        span = sl.stop - sl.start
        assert 0 < span < dim, f"slice {sl} should be a strict sub-range of {dim}"


def test_read_block_clips_to_bounds(tmp_path):
    uri, _ = make_local_zarr(tmp_path / "v.zarr")
    vol = open_volume(uri, 0)
    blk = read_block(vol, -5, 5, 0, 10, 0, 10)   # clipped to [0:5]
    assert blk.shape == (5, 10, 10)
    empty = read_block(vol, 100, 200, 0, 10, 0, 10)
    assert empty.size == 0


def test_s3_auth_auto_detection(monkeypatch):
    """anon=None auto-selects signed vs anonymous from the AWS env; explicit
    True/False forces; storage_options merge through (e.g. private endpoints)."""
    from scrollunwrap.zarr_source import (_AWS_CRED_ENV, _s3_storage_options,
                                          aws_credentials_available)
    for k in _AWS_CRED_ENV:
        monkeypatch.delenv(k, raising=False)
    assert aws_credentials_available() is False
    assert _s3_storage_options(None, None) == {"anon": True}        # no creds -> anon
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "AKIA_test")            # incl. STS sets this
    assert aws_credentials_available() is True
    assert _s3_storage_options(None, None) == {"anon": False}       # creds -> signed
    assert _s3_storage_options(True, None) == {"anon": True}        # force anon
    assert _s3_storage_options(False, {"endpoint_url": "http://minio:9000"}) == \
        {"anon": False, "endpoint_url": "http://minio:9000"}


def test_apply_threshold():
    arr = np.array([[0, 1, 127, 128, 255]], dtype=np.uint8)
    assert (apply_threshold(arr, None) == np.array([[0, 255, 255, 255, 255]])).all()
    assert (apply_threshold(arr, ">=128") == np.array([[0, 0, 0, 255, 255]])).all()
    assert (apply_threshold(arr, "==255") == np.array([[0, 0, 0, 0, 255]])).all()
    assert apply_threshold(arr, ">=1").dtype == np.uint8
