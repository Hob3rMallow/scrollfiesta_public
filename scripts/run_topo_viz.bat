@echo off
REM Run topology error heatmap visualization for all 7 training cubes.
REM Requires --metrics with TRAIN_LABELS directory and --dump-obj for OBJ output.

set EXE=build\Release\cube_mesh.exe
set CUBES=86701140 118632705 1006462223 1013184726 3290306825 3294954456 3394433588

for %%c in (%CUBES%) do (
    echo === %%c ===
    %EXE% nnunet-preds\%%c.tif output\topo_viz\%%c.tif ^
        --dump-obj output\topo_viz --metrics TRAIN_LABELS --no-timeout
    echo.
)

echo Done.
