@echo off
REM run_all.bat — Batch runner for vesuvius-c pipeline on Windows.
REM
REM Usage:
REM   run_all.bat [input_dir] [output_dir]
REM
REM Defaults:
REM   input_dir  = nnunet-preds (relative to repo root)
REM   output_dir = output       (relative to repo root)
REM
REM Runs cube_mesh.exe on each .tif file in input_dir.
REM Results logged to output_dir\results.csv.

setlocal enabledelayedexpansion

REM Resolve repo root (parent of scripts/)
set "SCRIPT_DIR=%~dp0"
set "REPO_ROOT=%SCRIPT_DIR%.."

REM Find executable
set "EXE=%REPO_ROOT%\build\Release\cube_mesh.exe"
if not exist "%EXE%" (
    echo ERROR: %EXE% not found.
    echo Build with: /vs build Release
    exit /b 1
)

REM Parse arguments
set "INPUT_DIR=%~1"
set "OUTPUT_DIR=%~2"
if "%INPUT_DIR%"=="" set "INPUT_DIR=%REPO_ROOT%\nnunet-preds"
if "%OUTPUT_DIR%"=="" set "OUTPUT_DIR=%REPO_ROOT%\output"

if not exist "%INPUT_DIR%" (
    echo ERROR: Input directory not found: %INPUT_DIR%
    exit /b 1
)
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

REM CSV header
set "CSV=%OUTPUT_DIR%\results.csv"
echo cube_id,total_time,s0,s1,s2,s3,s4,s5,s6,status> "%CSV%"

set TOTAL=0
set OK=0
set FAIL=0

for %%f in ("%INPUT_DIR%\*.tif") do (
    set /a TOTAL+=1
    set "CUBE_ID=%%~nf"
    set "OUT_TIF=%OUTPUT_DIR%\!CUBE_ID!.tif"
    set "LOG=%OUTPUT_DIR%\!CUBE_ID!.log"

    echo === [!TOTAL!] Processing !CUBE_ID! ===

    "%EXE%" "%%f" "!OUT_TIF!" --dump-obj "%OUTPUT_DIR%" 2>"!LOG!"
    set "RC=!errorlevel!"

    if !RC! equ 0 (
        set /a OK+=1
        set "STATUS=OK"
    ) else (
        set /a FAIL+=1
        set "STATUS=FAILED"
    )

    REM Parse step timings from log
    set "S0=" & set "S1=" & set "S2=" & set "S3=" & set "S4=" & set "S5=" & set "S6=" & set "TOTAL_TIME="
    for /f "tokens=2-8 delims= " %%a in ('findstr "Step timings:" "!LOG!" 2^>nul') do (
        for /f "tokens=2 delims==" %%x in ("%%a") do set "S0=%%x"
        for /f "tokens=2 delims==" %%x in ("%%b") do set "S1=%%x"
        for /f "tokens=2 delims==" %%x in ("%%c") do set "S2=%%x"
        for /f "tokens=2 delims==" %%x in ("%%d") do set "S3=%%x"
        for /f "tokens=2 delims==" %%x in ("%%e") do set "S4=%%x"
        for /f "tokens=2 delims==" %%x in ("%%f") do set "S5=%%x"
        for /f "tokens=2 delims==" %%x in ("%%g") do set "S6=%%x"
    )
    for /f "tokens=2 delims= " %%t in ('findstr "Total:" "!LOG!" 2^>nul') do (
        set "TOTAL_TIME=%%t"
        set "TOTAL_TIME=!TOTAL_TIME:s=!"
    )

    echo !CUBE_ID!,!TOTAL_TIME!,!S0!,!S1!,!S2!,!S3!,!S4!,!S5!,!S6!,!STATUS!>> "%CSV%"
    echo   -^> !STATUS!  !TOTAL_TIME!
)

echo.
echo ===============================
echo Summary
echo ===============================
echo   Total cubes:   !TOTAL!
echo   OK:            !OK!
echo   Failed:        !FAIL!
echo   Results:       %CSV%
echo ===============================
