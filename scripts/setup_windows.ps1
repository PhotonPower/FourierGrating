<#
Installs build dependencies of FourierGrating on Windows (MSYS2 / MinGW-w64, UCRT64).
  .\scripts\setup_windows.ps1              installs MSYS2 (if missing) and the UCRT64 toolchain
  .\scripts\setup_windows.ps1 -SkipMsys2Install   assumes MSYS2 is already installed at C:\msys64

Native Windows build, CPU only (no CUDA backend). The reference target system per AGENTS.md
remains AlmaLinux 10 / GCC 14; this is an unofficial CPU-only path for Windows contributors.
#>
[CmdletBinding()]
param(
    [switch]$SkipMsys2Install
)
$ErrorActionPreference = "Stop"

$msys2Root = "C:\msys64"
$msys2Bash = "$msys2Root\usr\bin\bash.exe"

if (-not (Test-Path $msys2Bash)) {
    if ($SkipMsys2Install) {
        throw "MSYS2 not found at $msys2Root and -SkipMsys2Install was given."
    }
    Write-Host "MSYS2 not found - installing via winget ..."
    winget install -e --id MSYS2.MSYS2 --accept-package-agreements --accept-source-agreements
}

if (-not (Test-Path $msys2Bash)) {
    throw "MSYS2 still not found at $msys2Root after install attempt. Install it manually from https://www.msys2.org/ and re-run this script."
}

Write-Host "Updating MSYS2 package database ..."
& $msys2Bash -lc "pacman -Sy --noconfirm"

$packages = @(
    "mingw-w64-ucrt-x86_64-gcc",
    "mingw-w64-ucrt-x86_64-cmake",
    "mingw-w64-ucrt-x86_64-ninja",
    "mingw-w64-ucrt-x86_64-pkgconf",
    "mingw-w64-ucrt-x86_64-eigen3",
    "mingw-w64-ucrt-x86_64-nlohmann-json",
    "mingw-w64-ucrt-x86_64-gtest",
    "mingw-w64-ucrt-x86_64-python",
    "mingw-w64-ucrt-x86_64-python-numpy",
    "mingw-w64-ucrt-x86_64-python-matplotlib"
)
Write-Host "Installing: $($packages -join ' ')"
& $msys2Bash -lc "pacman -S --noconfirm --needed $($packages -join ' ')"

Write-Host @"

Done. Next (PATH must include $msys2Root\ucrt64\bin for the mingw runtime DLLs):
  `$env:PATH = "$msys2Root\ucrt64\bin;" + `$env:PATH
  cmake --preset win-release
  cmake --build --preset win-release
  ctest --preset win-release

  .\build\win-release\apps\fgsolve.exe --info
"@
