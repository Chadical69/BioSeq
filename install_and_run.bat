@echo off
setlocal
cd /d "%~dp0"

echo ========================================
echo   BioSeq - DNA Sequence Analysis Tool
echo ========================================
echo.

set "PYTHON="
where py >nul 2>nul
if not errorlevel 1 set "PYTHON=py"
if not defined PYTHON (
    where python >nul 2>nul
    if not errorlevel 1 set "PYTHON=python"
)
if not defined PYTHON goto :no_python

if not exist ".venv\Scripts\python.exe" (
    echo [1/3] Creating virtual environment...
    %PYTHON% -m venv .venv
    if errorlevel 1 goto :error
)

echo [2/3] Installing or updating dependencies...
".venv\Scripts\python.exe" -m pip install --upgrade pip
if errorlevel 1 goto :error
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto :error

echo.
echo [3/3] Starting BioSeq...
echo The browser should open automatically. Keep this window open while using the app.
echo.
".venv\Scripts\python.exe" -m streamlit run app.py
if errorlevel 1 goto :error
exit /b 0

:no_python
echo Python was not found. Install Python 3.10+ from https://www.python.org/downloads/
echo Be sure to select "Add Python to PATH" during installation.
pause
exit /b 1

:error
echo.
echo Something went wrong. Check the error messages above and try again.
pause
exit /b 1
