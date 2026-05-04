@echo off
setlocal

call :run_lbmtv
if errorlevel 1 exit /b %errorlevel%

call :run_fdmadv
if errorlevel 1 exit /b %errorlevel%

call :run_fdlbm
if errorlevel 1 exit /b %errorlevel%

exit /b 0

:run_lbmtv
call "%~dp0build_lbmtv.cmd"
if errorlevel 1 exit /b %errorlevel%

if not exist outputs\sec1\lbmtv mkdir outputs\sec1\lbmtv

pushd outputs\sec1\lbmtv
echo Running ..\..\..\build\bin\lbmtv.exe
..\..\..\build\bin\lbmtv.exe > run.log
if errorlevel 1 (
  set "TEST_EXIT=%ERRORLEVEL%"
  popd
  exit /b %TEST_EXIT%
)
popd

if not exist outputs\sec1\lbmtv\error (
  echo Expected output file was not generated: error
  exit /b 1
)

findstr /b /c:"u " outputs\sec1\lbmtv\error >nul 2>nul
if errorlevel 1 (
  echo The error file does not contain the expected u entry.
  type outputs\sec1\lbmtv\error
  exit /b 1
)

exit /b 0

:run_fdmadv
call "%~dp0build_one.cmd" "src\sec2\fdmadv.c"
if errorlevel 1 exit /b %errorlevel%

if not exist outputs\sec2\fdmadv mkdir outputs\sec2\fdmadv

pushd outputs\sec2\fdmadv
echo Running ..\..\..\build\bin\fdmadv.exe with Upwind scheme
echo 1| ..\..\..\build\bin\fdmadv.exe > run.log
if errorlevel 1 (
  set "TEST_EXIT=%ERRORLEVEL%"
  popd
  exit /b %TEST_EXIT%
)
popd

if not exist outputs\sec2\fdmadv\fdmadv (
  echo Expected output file was not generated: fdmadv
  exit /b 1
)

exit /b 0

:run_fdlbm
call "%~dp0build_one.cmd" "src\sec2\fdlbm.c"
if errorlevel 1 exit /b %errorlevel%

if not exist outputs\sec2\fdlbm mkdir outputs\sec2\fdlbm

pushd outputs\sec2\fdlbm
echo Running ..\..\..\build\bin\fdlbm.exe
..\..\..\build\bin\fdlbm.exe > run.log
if errorlevel 1 (
  set "TEST_EXIT=%ERRORLEVEL%"
  popd
  exit /b %TEST_EXIT%
)
popd

if not exist outputs\sec2\fdlbm\data (
  echo Expected output file was not generated: data
  exit /b 1
)

exit /b 0
