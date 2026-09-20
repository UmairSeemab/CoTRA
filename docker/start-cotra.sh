#!/usr/bin/env bash
set -euo pipefail

HOST="${COTRA_HOST:-0.0.0.0}"
PORT="${COTRA_PORT:-3838}"
DATA_DIR="${COTRA_DATA_DIR:-/data}"
RESULTS_DIR="${COTRA_RESULTS_DIR:-/results}"
APP_DIR="/tmp/CoTRA_app"

mkdir -p "$DATA_DIR" "$RESULTS_DIR"
rm -rf "$APP_DIR"

# CoTRA currently defaults output to ~/CoTRA_Results. Setting HOME to the
# mounted results directory makes that default persistent without changing
# normal non-container CoTRA behaviour.
export HOME="$RESULTS_DIR"

exec R --vanilla -q -e "options(shiny.host='${HOST}', shiny.port=as.integer('${PORT}')); library(CoTRA); CoTRA::runCoTRA(app_dir='${APP_DIR}', install_missing=FALSE, ask=FALSE, launch.browser=FALSE)"
