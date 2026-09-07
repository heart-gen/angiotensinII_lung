#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-small
#SBATCH --job-name=peri_stats
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH --output=logs/analyze_pericytes.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Bridges-2 info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## The batch shell does not source the login profile, so `module` and `conda` must
## be bootstrapped here rather than inherited from the submitting shell.
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** Run analysis ****"
Rscript ../_h/04.pericytes_disease_analysis.R

if [ $? -ne 0 ]; then
    ## A `_FAILED` MARKER, not just a message. This step writes four TSVs and a
    ## figure BEFORE its final return statement, so when it threw on 2026-01-05
    ## the output directory looked complete and nothing on disk said otherwise --
    ## a reader could, and did, treat those tables as finished results (P1-17).
    ## An exit message only helps someone reading the log; the marker sits next to
    ## the outputs, where the person quoting them is looking.
    log_message "Error: Rscript execution failed (04.pericytes_disease_analysis.R)"
    printf '%s\n' \
        "This directory's job EXITED NON-ZERO. Any files here were written BEFORE" \
        "the failure and are INCOMPLETE. Do not cite them." \
        "job: ${SLURM_JOBID:-<interactive>}   when: $(date '+%Y-%m-%d %H:%M:%S')" \
        "script: ../_h/04.pericytes_disease_analysis.R" \
        > ./pericyte_subclusters/_FAILED 2>/dev/null \
        || printf 'job %s failed at %s\n' "${SLURM_JOBID:-<interactive>}" "$(date)" > ./_FAILED
    exit 1
fi

## Clear a stale marker from a previous failed run once the step succeeds, so the
## marker means "the LAST run failed" rather than "some run once failed".
rm -f ./pericyte_subclusters/_FAILED ./_FAILED

conda deactivate
log_message "Job finished at: $(date)"
