@echo off
REM Build the seam-weld unit test. Run from any dir.
REM   build_seam_weld_test.bat        -- compile only
REM Output: output\test_seam_weld\seam_weld_test.exe (one dir per experiment)
setlocal
call "C:\Program Files\Microsoft Visual Studio\2022\Professional\Common7\Tools\VsDevCmd.bat" -arch=x64 -no_logo
cd /d D:\work\vesuvius-c
if not exist output\test_seam_weld mkdir output\test_seam_weld
cl /nologo /W4 /Zi /Fe:output\test_seam_weld\seam_weld_test.exe /Fo:output\test_seam_weld\ ^
   scripts\step0-mesh-extract\seam_weld_test.c ^
   src\remesh\seam_weld.c src\remesh\ball_pivot.c src\remesh\pinhole_fill.c ^
   src\remesh\orient_weld.c ^
   src\common\arena.c src\common\except.c
exit /b %errorlevel%
