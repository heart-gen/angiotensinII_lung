#!/bin/bash

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Local info ****"
echo "User: ${USER}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility

log_message "**** Run analysis ****"

Rscript ../_h/05.colocalization_plot.R

if [ $? -ne 0 ]; then
    log_message "Error: Rscript execution failed"
    exit 1
fi

log_message "**** Job ends ****"
