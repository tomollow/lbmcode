@echo off
setlocal

if not exist build mkdir build

call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
if errorlevel 1 exit /b %errorlevel%

cl /nologo /D_USE_MATH_DEFINES /Febuild\lbmtv.exe src\sec1\lbmtv.c /link /STACK:8388608