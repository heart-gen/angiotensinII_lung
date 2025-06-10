#!/bin/bash

# Set the log file paths
LOG_FILE="summary.log"
ERROR_LOG="summary.error"

# Print some information about the job
echo "Start time: $(date)" > $LOG_FILE

# Activate the virtual environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate r44

# Check if the environment was activated successfully
if [ $? -ne 0 ]; then
    echo "Failed to activate virtual environment. Exiting." >> $ERROR_LOG
    exit 1
fi

echo "Virtual environment 'R version 4.4' activated successfully." >> $LOG_FILE

# Run script
Rscript ../_h/01.detect_power.R >> $LOG_FILE 2>&1

# Check if the script ran successfully
if [ $? -ne 0 ]; then
    echo "R script failed. Check the error logs." >> $ERROR_LOG
    exit 1
fi

# Deactivate the virtual environment
conda deactivate

echo "Job finished at: $(date)" >> $LOG_FILE
