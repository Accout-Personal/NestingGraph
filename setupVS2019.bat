@echo off
setlocal EnableExtensions

set "PRESET=vs-static"
set "ROOT=%~dp0"
set "ROOT=%ROOT:~0,-1%"
pushd "%ROOT%" >nul

rem --- Find vswhere ---
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" (
  echo ERROR: vswhere.exe not found. Visual Studio Installer is missing.
  goto :fail
)

rem --- Find VS2019-bundled CMake directly (no guessing, no VSROOT concat) ---
set "CMAKE_EXE="
for /f "delims=" %%i in ('
  "%VSWHERE%" -latest -version 16 -products * ^
    -find Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe 2^>nul
') do set "CMAKE_EXE=%%i"

rem --- Fallback: use standalone CMake if VS CMake isn't installed ---
if not defined CMAKE_EXE (
  for /f "delims=" %%i in ('where cmake 2^>nul') do (
    set "CMAKE_EXE=%%i"
    goto :have_cmake
  )
)

:have_cmake
if not defined CMAKE_EXE (
  echo ERROR: Could not find CMake.
  echo Install "C++ CMake tools for Windows" in VS2019 OR install standalone CMake.
  goto :fail
)

echo Using CMake:
echo   %CMAKE_EXE%
echo.

rem --- Configure = generate .sln/.vcxproj using your preset ---
"%CMAKE_EXE%" --preset "%PRESET%"
if errorlevel 1 goto :fail

rem Your preset binaryDir is build-static
set "BUILD_DIR=%ROOT%\build-static"

rem --- Open the generated solution (prefer .sln if present) ---
set "SLN="
for /f "delims=" %%s in ('dir /b "%BUILD_DIR%\*.sln" 2^>nul') do (
  set "SLN=%%s"
  goto :open_it
)

:open_it
if defined SLN (
  echo Opening solution:
  echo   %BUILD_DIR%\%SLN%
  start "" "%BUILD_DIR%\%SLN%"
) else (
  echo No .sln found, opening build folder:
  explorer "%BUILD_DIR%"
)

popd >nul
endlocal
exit /b 0

:fail
echo.
echo FAILED.
popd >nul
endlocal
pause
exit /b 1
