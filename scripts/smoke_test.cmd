@echo off
setlocal

call "%~dp0build_lbmtv.cmd"
if errorlevel 1 exit /b %errorlevel%

if not exist outputs\lbmtv mkdir outputs\lbmtv

pushd outputs\lbmtv
echo Running ..\..\build\bin\lbmtv.exe
..\..\build\bin\lbmtv.exe
if errorlevel 1 exit /b %errorlevel%
popd

if not exist outputs\lbmtv\error (
  echo Expected output file was not generated: error
  exit /b 1
)

findstr /b /c:"u " outputs\lbmtv\error >nul 2>nul
if errorlevel 1 (
  echo The error file does not contain the expected u entry.
  type outputs\lbmtv\error
  exit /b 1
)
