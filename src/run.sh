#!/bin/bash

# Script to run analyse_data.py 10 times, logging output and error with timestamp

for i in {1..10}
do
    timestamp=$(date +"%Y%m%d_%H%M%S")
    logfile="analyse_data_run_${timestamp}_$i.log"
    echo "[INFO] Run $i - Logging to $logfile"
    python3 analyse_data.py > "$logfile" 2>&1
    sleep 1  # Ensure unique timestamps

done
