#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ras_circuit_1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --output=logs/ras_circuit_step1.log

## Step 1 of ras_circuit (Figure 5). Submit from ras_circuit/_m:
##     sbatch -D ras_circuit/_m ras_circuit/_h/step_1.sh
##   (a) receiver programmes + the module gene panel
##   (b) donor x cell-type pseudobulk of that panel over the CCC niche (Figure 5C nodes)
##   (c) Figure 5B: four niche-affinity axes + detection-matched gene pseudobulk
##   (d) Figure 5B statistics

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"; echo "Node: ${SLURM_NODENAME}"
source /etc/profile.d/modules.sh
module purge; module load anaconda3/2024.10-1; module list
eval "$(conda shell.bash hook)"
mkdir -p logs stats_data

conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** (a) receiver programmes + module panel ****"
python ../_h/receiver_programs.py --out ./receiver_programs.tsv \
       --targets-out ./outgoing_targets.tsv || { log_message "Error: programmes"; exit 1; }
python - <<'EOF' || { echo "Error: panel build"; exit 1; }
import pandas as pd
ras = pd.read_csv("../../agt_axis/_h/ras_panel.tsv", sep="\t")
prog = pd.read_csv("./receiver_programs.tsv", sep="\t")
rows = [(g, r) for g, r in zip(ras["gene"], ras.iloc[:, 1] if ras.shape[1] > 1 else ["ras"] * len(ras))]
rows += [(g, "program:" + p) for g, p in zip(prog["gene"], prog["program"])]
panel = pd.DataFrame(rows, columns=["gene", "role"]).drop_duplicates("gene")
panel.to_csv("./ras_circuit_panel.tsv", sep="\t", index=False)
print(f"ras_circuit_panel.tsv: {len(panel)} genes")
EOF

log_message "**** (b) donor x cell-type pseudobulk ****"
python ../../basement_membrane/_h/02.niche_pseudobulk.py \
       --adata ../../cell_communication/_m/ccc_niche.h5ad \
       --genes ./ras_circuit_panel.tsv \
       --outfile ./ras_circuit_pseudobulk.tsv.gz || { log_message "Error: pseudobulk"; exit 1; }

log_message "**** (c) niche affinity ****"
python ../_h/01.niche_affinity.py \
       --airspace ../../localization/airspace_analysis/_m/pericytes_with_airspace_score.h5ad \
       --pericytes ../../pericyte_states/_m/pericyte_states.h5ad \
       --outdir ./ || { log_message "Error: niche affinity"; exit 1; }
conda deactivate

log_message "**** (d) niche-affinity statistics ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/02.niche_affinity_stats.R --outdir ./stats_data --min-cells 10 \
    || { log_message "Error: niche-affinity stats"; exit 1; }
conda deactivate
log_message "**** Job ends ****"
