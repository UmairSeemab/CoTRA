#!/usr/bin/env bash
set -euo pipefail

command -v docker >/dev/null 2>&1 || {
  echo "Docker was not found. Install Docker Engine first."
  exit 1
}

docker info >/dev/null 2>&1 || {
  echo "Docker is installed but the Docker engine is not running, or your user cannot access it."
  exit 1
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$ROOT_DIR/CoTRA_data"
RESULTS_DIR="$ROOT_DIR/CoTRA_results"

mkdir -p "$DATA_DIR" "$RESULTS_DIR"
docker pull ghcr.io/umairseemab/cotra:latest

docker rm -f cotra >/dev/null 2>&1 || true

docker run -d --rm --name cotra \
  -p 3838:3838 \
  -v "$DATA_DIR:/data" \
  -v "$RESULTS_DIR:/results" \
  ghcr.io/umairseemab/cotra:latest >/dev/null

echo "Starting CoTRA..."
for _ in $(seq 1 90); do
  status="$(docker inspect --format '{{.State.Health.Status}}' cotra 2>/dev/null || true)"
  if [ "$status" = "healthy" ]; then
    if command -v xdg-open >/dev/null 2>&1; then
      xdg-open http://localhost:3838 >/dev/null 2>&1 || true
    fi
    echo "CoTRA is running at http://localhost:3838"
    echo "Results: $RESULTS_DIR"
    echo "Stop with: docker stop cotra"
    exit 0
  fi
  if [ "$status" = "unhealthy" ]; then
    docker logs cotra
    exit 1
  fi
  sleep 2
done

docker logs cotra
exit 1
