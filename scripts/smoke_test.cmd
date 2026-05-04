@echo off
setlocal

call "%~dp0build_lbmtv.cmd"
if errorlevel 1 exit /b %errorlevel%

echo Running build\lbmtv.exe
build\lbmtv.exe
if errorlevel 1 exit /b %errorlevel%

if not exist error (
  echo Expected output file was not generated: error
  exit /b 1
)

findstr /b /c:"u " error >nul 2>nul
if errorlevel 1 (
  echo The error file does not contain the expected u entry.
  type error
  exit /b 1
)
