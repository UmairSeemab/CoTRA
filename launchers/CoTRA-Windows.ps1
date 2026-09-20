$ErrorActionPreference = "Stop"

if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
    throw "Docker was not found. Install Docker Desktop and start Docker first."
}

docker info | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "Docker is installed but the Docker engine is not running. Start Docker Desktop first."
}

$RootDir = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$DataDir = Join-Path $RootDir "CoTRA_data"
$ResultsDir = Join-Path $RootDir "CoTRA_results"

New-Item -ItemType Directory -Force -Path $DataDir | Out-Null
New-Item -ItemType Directory -Force -Path $ResultsDir | Out-Null

docker pull ghcr.io/umairseemab/cotra:latest

docker rm -f cotra 2>$null | Out-Null

docker run -d --rm --name cotra `
    -p 3838:3838 `
    -v "${DataDir}:/data" `
    -v "${ResultsDir}:/results" `
    ghcr.io/umairseemab/cotra:latest | Out-Null

Write-Host "Starting CoTRA..."
for ($i = 0; $i -lt 90; $i++) {
    $status = docker inspect --format '{{.State.Health.Status}}' cotra 2>$null
    if ($status -eq "healthy") {
        Start-Process "http://localhost:3838"
        Write-Host "CoTRA is running at http://localhost:3838"
        Write-Host "Results: $ResultsDir"
        Write-Host "Stop with: docker stop cotra"
        exit 0
    }
    if ($status -eq "unhealthy") {
        docker logs cotra
        throw "The CoTRA container became unhealthy. See the log above."
    }
    Start-Sleep -Seconds 2
}

docker logs cotra
throw "CoTRA did not become healthy. See the log above."
