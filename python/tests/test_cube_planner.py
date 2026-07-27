from scrollunwrap.cube_planner import (cube_id_str, plan_cubes,
                                       snap_bbox_to_level_grid)


def test_cube_id_format():
    assert cube_id_str(128, 256, 384) == "z00128_y00256_x00384"


def test_snap_level0():
    assert snap_bbox_to_level_grid((5, 30, 5, 30, 5, 30), 0, 8) == (0, 32, 0, 32, 0, 32)


def test_plan_origins_multiple_of_cube():
    cubes = plan_cubes((8, 24, 8, 24, 8, 24), level=0, cube=8)
    assert len(cubes) == 8
    for c in cubes:
        assert c.oz % 8 == 0 and c.oy % 8 == 0 and c.ox % 8 == 0
    assert cubes[0].cube_id == "z00008_y00008_x00008"


def test_level_downsample_divides_bbox():
    # level-1: a 256-voxel level-0 extent -> 128 at level-1 -> exactly one cube
    cubes = plan_cubes((0, 256, 0, 256, 0, 256), level=1, cube=128)
    assert len(cubes) == 1 and cubes[0].oz == 0 and cubes[0].oy == 0 and cubes[0].ox == 0


def test_origins_never_negative():
    cubes = plan_cubes((-50, 300, 0, 300, 0, 300), level=0, cube=128)
    assert all(c.oz >= 0 and c.oy >= 0 and c.ox >= 0 for c in cubes)
