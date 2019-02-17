@echo off

:: This file opens a terminal which allows you to compile the code with Microsoft Visual Studio 2015 (x64) on Windows
::
:: HOW TO USE THIS FILE ?
::
::   [run this file]
::   mkdir build
::   cd build
::   cmake -A x64 ..
::   cmake --build . --config Release
::   [executables are built in the bin/Release folder]
::   ctest -C Release

:: set the location of gmsh SDK ( **MODIFY THIS LINE FOR YOUR SYSTEM** )
set GMSHSDK=F:\local\gmsh-sdk

:: where is gmsh.exe and gmsh-**.dll ? (HINT: copy gmsh-**.dll to the bin folder)
set PATH=%GMSHSDK%\bin;%PATH%

:: where is gmsh.h ? (rename gmsh.h_cwrap => gmsh.h)
set INCLUDE=%GMSHSDK%\include;%INCLUDE%

:: where is gmsh.lib ?
set LIB=%GMSHSDK%\lib;%LIB%

:: where is gmsh.py ? (required only if you want to use the python API)
set PYTHONPATH=%GMSHSDK%\lib;%PYTHONPATH%

:: set the environment of the msvc compiler
CD /d "%~dp0"
CD ..
%comspec% /K ""C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" %processor_architecture%"