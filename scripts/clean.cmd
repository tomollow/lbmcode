@echo off
setlocal

if exist build rmdir /s /q build
if exist outputs rmdir /s /q outputs

del /q *.obj 2>nul
del /q *.exe 2>nul
del /q error 2>nul
del /q data* 2>nul