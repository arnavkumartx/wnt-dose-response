@echo off
setlocal
REM ---------------------------------------------------------------------------
REM Launch JupyterLab rooted at THIS repository, whatever directory you run it
REM from - including double-clicking this file.
REM
REM Why this exists: `python -m jupyterlab` with no arguments roots the file
REM browser at the shell's current directory. Run it from another project and
REM Jupyter shows that project instead of this one, which looks like Jupyter
REM opening "the wrong app". %~dp0 is this script's own folder, so the root is
REM pinned to the model regardless of where you are standing.
REM ---------------------------------------------------------------------------

set "REPO=%~dp0"
if "%REPO:~-1%"=="\" set "REPO=%REPO:~0,-1%"

REM Prefer the known interpreter over whatever `python` resolves to - on Windows
REM that is often the Microsoft Store stub, which is not a working Python.
set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if not exist "%PY%" set "PY=python"

echo Serving from : %REPO%
echo Opening in   : %REPO%\notebooks
echo.

"%PY%" -m jupyterlab ^
  --ServerApp.root_dir="%REPO%" ^
  --ServerApp.preferred_dir="%REPO%\notebooks" ^
  %*

endlocal
