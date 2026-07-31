[CmdletBinding()]
param(
    [string]$Toolchain = "1.92.0"
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$toolsRoot = Join-Path $repo ".toolchains"
$cargoHome = Join-Path $toolsRoot "cargo"
$rustupHome = Join-Path $toolsRoot "rustup"
$rustupInit = Join-Path $toolsRoot "rustup-init.exe"

New-Item -ItemType Directory -Force -Path $toolsRoot | Out-Null
if (-not (Test-Path -LiteralPath $rustupInit)) {
    $uri = "https://static.rust-lang.org/rustup/dist/x86_64-pc-windows-msvc/rustup-init.exe"
    Write-Host "Downloading rustup-init from $uri"
    Invoke-WebRequest -UseBasicParsing -Uri $uri -OutFile $rustupInit
}

$env:CARGO_HOME = $cargoHome
$env:RUSTUP_HOME = $rustupHome
& $rustupInit -y --profile minimal --default-host x86_64-pc-windows-msvc `
    --default-toolchain $Toolchain --no-modify-path
if ($LASTEXITCODE -ne 0) { throw "rustup-init failed with exit code $LASTEXITCODE" }

$rustup = Join-Path $cargoHome "bin\rustup.exe"
& $rustup component add rustfmt clippy --toolchain $Toolchain
if ($LASTEXITCODE -ne 0) { throw "rustup component install failed with exit code $LASTEXITCODE" }

Write-Host "Project-local Rust installed:"
& (Join-Path $cargoHome "bin\rustc.exe") --version
Write-Host "CARGO_HOME=$cargoHome"
Write-Host "RUSTUP_HOME=$rustupHome"
