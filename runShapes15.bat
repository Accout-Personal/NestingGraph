@echo off
setlocal EnableExtensions

REM Root folder = directory containing this .bat
set "ROOT=%~dp0"
REM Remove trailing backslash for cleaner joins (optional)
if "%ROOT:~-1%"=="\" set "ROOT=%ROOT:~0,-1%"

set "EXE=%ROOT%\build-static\Release\nesting_graph.exe"
set "DATASETS=%ROOT%\Inputs\dataset"
set "INSTANCES=%ROOT%\Inputs\Instances\shapes15.txt"

REM Existence checks
if not exist "%EXE%" (
  echo EXE not found: "%EXE%"
  exit /b 1
)
if not exist "%DATASETS%" (
  echo Datasets path not found: "%DATASETS%"
  exit /b 1
)
if not exist "%INSTANCES%" (
  echo Instances file not found: "%INSTANCES%"
  exit /b 1
)

echo Running: "%EXE%"
echo   --datasets  "%DATASETS%"
echo   --instances "%INSTANCES%"
echo   --typeOriented 1
echo   --CliqueCovering 1
echo.

REM Run (merge stderr into stdout)
"%EXE%" --datasets "%DATASETS%" --instances "%INSTANCES%" --typeOriented 0 --cliqueCovering 1 --apply_cuts 1 2>&1

set "CODE=%ERRORLEVEL%"
echo.
echo Exit code: %CODE%
pause
exit /b %CODE%
