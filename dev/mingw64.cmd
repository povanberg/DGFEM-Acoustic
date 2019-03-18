@echo off

:: This file opens a terminal which allows you to compile the code with a 64-bit mingw compiler
::
::     https://mingw-w64.org/
::     https://www.msys2.org/
::
:: HOW TO USE THIS FILE ?
::
::   [check the PATHs below]
::   [run this file]
::   mkdir build
::   cd build
::   cmake -G "MinGW Makefiles" ..
::   mingw32-make
::   ctest
::
:: How to clean the "build" folder using cmd line ?
::
::   cd build
::   rd /q /s .

echo setting MinGW64 environment...

:: set the location of gmsh SDK ( **MODIFY THIS LINE FOR YOUR SYSTEM** )
set GMSHSDK=C:\local\gmsh-4.1.4-Windows64-sdk

:: where is gmsh.exe and gmsh-**.dll ? (HINT: copy gmsh-**.dll to the bin folder)
set PATH=%GMSHSDK%\bin;%PATH%
set PATH=C:\local\bin;%PATH%

:: where is gmsh.h ? (rename gmsh.h_cwrap => gmsh.h)
set INCLUDE=%GMSHSDK%\include;%INCLUDE%
set INCLUDE=C:\Eigen;%INCLUDE%
set INCLUDE=C:\local\include;%INCLUDE%

:: where is gmsh.lib ?
set LIB=%GMSHSDK%\lib;%LIB%
set LIB=C:\local\lib;%LIB%

:: where is gmsh.py ? (required only if you want to use the python API)
set PYTHONPATH=%GMSHSDK%\lib;%PYTHONPATH%

::set PATH=C:\mingw-w64\mingw64\bin;%PATH%
set PATH=C:\mingw-w64\mingw64\bin;%PATH%

:: open terminal
CD /d "%~dp0"
CD ..
%comspec% /K
