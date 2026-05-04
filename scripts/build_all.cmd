@echo off
setlocal

for /r src %%F in (*.c) do (
  echo Building %%F
  call "%~dp0build_one.cmd" "%%F" "build\%%~nF.exe"
  if errorlevel 1 exit /b %errorlevel%
)
