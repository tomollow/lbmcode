@echo off
setlocal

if "%~1"=="" (
  echo Usage: scripts\run_one.cmd source-file [output-dir]
  exit /b 1
)

set "SOURCE=%~1"
if not exist "%SOURCE%" (
  echo Source file not found: %SOURCE%
  exit /b 1
)

call "%~dp0build_one.cmd" "%SOURCE%"
if errorlevel 1 exit /b %errorlevel%

if "%~2"=="" (
  set "RUN_DIR=outputs\%~n1"
) else (
  set "RUN_DIR=%~2"
)

if not exist "%RUN_DIR%" mkdir "%RUN_DIR%"

for %%I in ("build\bin\%~n1.exe") do set "EXE_FULL=%%~fI"

pushd "%RUN_DIR%"
echo Running %EXE_FULL%
"%EXE_FULL%"
set "RUN_EXIT=%ERRORLEVEL%"
popd

exit /b %RUN_EXIT%