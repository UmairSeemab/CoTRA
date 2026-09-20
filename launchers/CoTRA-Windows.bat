@echo off
setlocal enabledelayedexpansion

where docker >nul 2>&1 || (
  echo Docker was not found. Install Docker Desktop and start Docker first.
  pause
  exit /b 1
)

docker info >nul 2>&1 || (
  echo Docker is installed but the Docker engine is not running.
  echo Start Docker Desktop and run this file again.
  pause
  exit /b 1
)

for %%I in ("%~dp0..") do set "ROOT_DIR=%%~fI"
set "DATA_DIR=%ROOT_DIR%\CoTRA_data"
set "RESULTS_DIR=%ROOT_DIR%\CoTRA_results"

if not exist "%DATA_DIR%" mkdir "%DATA_DIR%"
if not exist "%RESULTS_DIR%" mkdir "%RESULTS_DIR%"

docker pull ghcr.io/umairseemab/cotra:latest || goto :error

docker rm -f cotra >nul 2>&1

docker run -d --rm --name cotra -p 3838:3838 ^
  -v "%DATA_DIR%:/data" ^
  -v "%RESULTS_DIR%:/results" ^
  ghcr.io/umairseemab/cotra:latest >nul || goto :error

echo Starting CoTRA...
for /L %%G in (1,1,90) do (
  for /F "delims=" %%S in ('docker inspect --format "{{.State.Health.Status}}" cotra 2^>nul') do set "STATUS=%%S"
  if "!STATUS!"=="healthy" goto :ready
  if "!STATUS!"=="unhealthy" goto :unhealthy
  timeout /t 2 /nobreak >nul
)

goto :unhealthy

:ready
start "" http://localhost:3838
echo CoTRA is running at http://localhost:3838
echo Results: %RESULTS_DIR%
echo Stop with: docker stop cotra
pause
exit /b 0

:unhealthy
docker logs cotra
echo CoTRA did not become healthy. See the log above.
pause
exit /b 1

:error
echo Docker could not start CoTRA.
pause
exit /b 1
