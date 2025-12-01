#!/bin/bash

# Script to generate aggregated disease and drug scores with averages

# Activate the correct virtual environment (.venv as per repo guidelines)
source ../.venv/bin/activate

# Generate aggregated summary scores
echo "[INFO] Generating aggregated scores from simplified data..."
python3 aggregate_scores.py

echo "[INFO] Done! Aggregated scores available at ../data/results/aggregated_scores.csv"
