# conftest anchors the tests dir on sys.path so `import _helpers` works, and
# provides fixtures for the native binaries (skipping when they're not built).
from pathlib import Path

import pytest

from scrollunwrap import _paths


def _need(path, name):
    if not path or not Path(path).exists():
        pytest.skip(f"{name} not built/available at {path!r}")
    return str(path)


@pytest.fixture
def cube_mesh_bin():
    return _need(_paths.default_cube_mesh(), "cube_mesh")


@pytest.fixture
def grid_weld_bin():
    return _need(_paths.default_grid_weld(), "grid_weld")


@pytest.fixture
def flatboi_bin():
    return _need(_paths.default_flatboi(), "flatboi")


@pytest.fixture
def obj2tifxyz_bin():
    return _need(_paths.default_obj2tifxyz(), "vc_obj2tifxyz_legacy")
