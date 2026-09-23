#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ras_circuit_2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=16:00:00
#SBATCH --output=logs/ras_circuit_step2.log

## Step 2 of ras_circuit (Figures 5D, 5E). Submit from ras_circuit/_m:
##     sbatch -D ras_circuit/_m ras_circuit/_h/step_2.sh
##   00   fetch + cache decoupler networks (compute nodes have internet)
##   00b  derive the AngII signature from E-MTAB-8810 (rule: _h/signatures/README.md)
##   00c  mouse -> human 1:1 orthologs; FREEZE into _h/signatures/ if not frozen
##   03   score pericytes: signature + PROGENy/CollecTRI + legacy comparison
##   04   1,000 detection-matched null panels
##   05   Figure 5D statistics (count-model arbiter + null, scVI lens, continuum, subclusters)
##   06   Figure 5E statistics (BM - fibrillar contrast + null)

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"; echo "Node: ${SLURM_NODENAME}"
source /etc/profile.d/modules.sh
module purge; module load anaconda3/2024.10-1; module list
eval "$(conda shell.bash hook)"
mkdir -p logs stats_data networks

SIG=../_h/signatures/angii_response_signature.tsv
PERI=../../pericyte_states/_m/pericyte_states.h5ad

conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
log_message "**** 00 networks ****"
python ../_h/00.fetch_networks.py --outdir ./networks || { log_message "Error: 00"; exit 1; }
if [ ! -f "$SIG" ]; then
    log_message "**** 00b derive AngII signature (no frozen copy yet) ****"
    python ../_h/00b.derive_angii_signature.py --outdir ./ || { log_message "Error: 00b"; exit 1; }
    conda deactivate
    conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
    Rscript ../_h/00c.map_orthologs.R ./angii_signature/mouse_signature.tsv "$SIG" \
        || { log_message "Error: 00c"; exit 1; }
    conda deactivate
    conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
else
    log_message "frozen signature present, not re-derived: $SIG"
fi

log_message "**** 03 score ****"
python ../_h/03.at1r_response_score.py --adata "$PERI" --signature "$SIG" \
       --networks ./networks \
       --continuum ../../pericyte_states/_m/continuum_metadata.tsv.gz \
       --bm ../../basement_membrane/_m/bm_metadata.tsv.gz \
       --outdir ./ || { log_message "Error: 03"; exit 1; }
log_message "**** 04 null panels ****"
python ../_h/04.at1r_null_panels.py --adata "$PERI" \
       --audit ./stats_data/at1r_signature_audit.tsv --outdir ./ --n-null 1000 \
       || { log_message "Error: 04"; exit 1; }
conda deactivate

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
log_message "**** 05 Figure 5D statistics ****"
Rscript ../_h/05.at1r_response_stats.R --outdir ./stats_data || { log_message "Error: 05"; exit 1; }
log_message "**** 06 Figure 5E statistics ****"
Rscript ../_h/06.matrix_consequence_stats.R --outdir ./stats_data || { log_message "Error: 06"; exit 1; }
conda deactivate
log_message "**** Job ends ****"
