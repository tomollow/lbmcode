@echo off
setlocal

if "%~1"=="" (
  echo Usage: scripts\build_one.cmd source-file [output-exe]
  exit /b 1
)

set "SOURCE=%~1"
if not exist "%SOURCE%" (
  echo Source file not found: %SOURCE%
  exit /b 1
)

if "%~2"=="" (
  set "OUTPUT=build\bin\%~n1.exe"
) else (
  set "OUTPUT=%~2"
)

for %%I in ("%OUTPUT%") do set "OUTPUT_DIR=%%~dpI"
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

set "OBJECT=build\obj\%~n1.obj"
for %%I in ("%OBJECT%") do set "OBJECT_DIR=%%~dpI"
if not exist "%OBJECT_DIR%" mkdir "%OBJECT_DIR%"

where cl >nul 2>nul
if errorlevel 1 (
  if exist "%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" (
    for /f "usebackq delims=" %%I in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -find VC\Auxiliary\Build\vcvars64.bat`) do (
      call "%%I"
      goto :compiler_ready
    )
  )

  if exist "%ProgramFiles%\Microsoft Visual Studio\Installer\vswhere.exe" (
    for /f "usebackq delims=" %%I in (`"%ProgramFiles%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -find VC\Auxiliary\Build\vcvars64.bat`) do (
      call "%%I"
      goto :compiler_ready
    )
  )

  echo cl.exe was not found. Start a Developer Command Prompt or install Visual Studio C++ tools.
  exit /b 1
)

:compiler_ready
cl /nologo /D_USE_MATH_DEFINES /Fo"%OBJECT%" /Fe"%OUTPUT%" "%SOURCE%" /link /STACK:8388608
