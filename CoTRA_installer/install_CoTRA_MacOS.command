#!/bin/bash

echo "====================================="
echo "Starting CoTRA package installation"
echo "====================================="

cd "$(dirname "$0")"

Rscript install_cotra_packages.R

echo ""
echo "====================================="
echo "Installation finished"
echo "====================================="
echo ""

read -p "Press Enter to close..."
