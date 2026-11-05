@echo off
rem Runs a command inside the MSVC 2022 x64 developer environment (vcvars64), e.g.
rem   scripts\msvc.cmd cmake --preset win-msvc
rem   scripts\msvc.cmd cmake --build --preset win-msvc
rem   scripts\msvc.cmd ctest --preset win-msvc -L unit
setlocal
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSDIR=%%i"
if not defined VSDIR (
  echo MSVC 2022 with the x64 tools not found 1>&2
  exit /b 1
)
rem vcvars64.bat itself calls vswhere.exe without a path
set "PATH=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer;%PATH%"
call "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" >nul || exit /b 1
%*
