@echo off

:: This file opens a terminal which allows you to compile the code with the 32-bits mingw compiler provided with Code::Blocks on Windows
::
:: HOW TO USE THIS FILE ?
::
::   [run this file]
::   mkdir build
::   cd build
::   cmake -G "CodeBlocks - MinGW Makefiles" ..
::   [open the generated project in Code::Blocks]

:: set the location of gmsh SDK ( **MODIFY THIS LINE FOR YOUR SYSTEM** )
set GMSHSDK=C:\local\gmsh-sdk32

:: where is gmsh.exe and gmsh-**.dll ? (HINT: copy gmsh-**.dll to the bin folder)
set PATH=%GMSHSDK%\bin;%GMSHSDK%\lib;%PATH%

:: where is gmsh.h ? (rename gmsh.h_cwrap => gmsh.h)
set INCLUDE=%GMSHSDK%\include;%INCLUDE%

:: where is gmsh.lib ?
set LIB=%GMSHSDK%\lib;%LIB%

:: where is gmsh.py ? (required only if you want to use the python API)
set PYTHONPATH=%GMSHSDK%\lib;%PYTHONPATH%

:: set the environment of the Code::Blocks compiler (mingw)
CD /d "%~dp0"
CD ..
%comspec% /K ""C:\Program Files (x86)\CodeBlocks\MinGW\mingwvars.bat""
