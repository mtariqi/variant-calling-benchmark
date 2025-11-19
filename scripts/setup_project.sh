#!/bin/bash
# Complete project setup script

set -euo pipefail

echo "Setting up variant calling benchmark project..."

# Create conda environment
echo "Creating conda environment..."
conda env create -f environment.yml

# Activate environment
echo "Activating environment..."
conda activate variant-benchmark

# Make scripts executable
chmod +x scripts/*.sh

echo "Setup complete! Next steps:"
echo "1. Download data: bash scripts/download_data.sh"
echo "2. Test workflow: snakemake --cores 2 --use-conda --dry-run"
echo "3. Run full analysis: snakemake --cores all --use-conda"
