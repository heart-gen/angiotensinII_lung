#!/bin/bash
# Recovery script: submits everything still outstanding from the BM/AGT expansion.
#
# Only needed if the two background submitters died before finishing
# (check with: pgrep -f 'cascade2.sh|cross_species.sh').
# Safe to run more than once ONLY if you first confirm the job is not already
# queued -- it does not deduplicate. Check `squeue -u $USER` first.
#
# Already submitted as of the pause (do NOT resubmit these):
#   42496703 cogaps_validate      42496704 ccc_prepare
#   42496705 ccc_liana            42496714 ccc_nichenet
#   42496715 ccc_alluvial         42496716 ccc_donor_val
#
# Outstanding, in dependency order:
#   basement_membrane/step_2  (refresh BM pseudobulk against rebuilt ccc_niche)
#   agt_axis/step_0           (refresh RAS pseudobulk likewise)
#   basement_membrane/step_4  (BM NicheNet targets, if 42496332 did not finish)
#   agt_axis/step_1           (AGT vs ligands, if 42496123 did not finish)
#   figures/step_figures
#   cross_species/step_1 -> step_2 -> step_3   (mouse scVI re-run with BM panel)

set -u
cd /ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung
ACCT=bio250020p
CC0=42496704   # ccc_prepare  -- rebuilds ccc_niche.h5ad
CC1=42496705   # ccc_liana
CC2=42496714   # ccc_nichenet

wait_slot() { while [ "$(squeue -u "$USER" -h -o %i | wc -l)" -ge 9 ]; do sleep 60; done; }
sub() {  # sub <module> <step> [dep_jobid]
  wait_slot
  local id
  id=$(sbatch --parsable --account=$ACCT -D "$1/_m" \
       ${3:+--dependency=afterok:$3} "$1/_h/$2" 2>&1)
  echo "$1/$2 -> $id"
  echo "$id"
}

BM2=$(sub basement_membrane step_2.sh "$CC0" | tail -1)
AG0=$(sub agt_axis          step_0.sh "$CC0" | tail -1)
sub basement_membrane step_4.sh "$CC1" >/dev/null
sub agt_axis          step_1.sh "$CC2" >/dev/null
sub figures           step_figures.sh "$CC2" >/dev/null

X1=$(sub cross_species step_1.sh          | tail -1)
X2=$(sub cross_species step_2.sh "$X1"    | tail -1)
sub cross_species step_3.sh "$X2" >/dev/null

echo "ALL REMAINING SUBMITTED"
