#!/usr/bin/env bash
#
set -euo pipefail

if [ -f "04a-feature_counts.log" ]
then
    RESUME_PARAM="--batch-resume 04a-feature_counts.log"
else
    RESUME_PARAM=""
fi

tail -n+2 input/sra_rnaseq.tsv \
| awk -F '\t' '$14 == "TRUE" || $14 == "NETWORK_ONLY"' \
| ../code/slurm_scripts/bin/pt --file - "code/feature_count.sh {0} {6} {7} input/Sscl1980-nuclear.fasta work/genes.gtf"  \
| ../code/slurm_scripts/bin/sbatch_jobarray.sh \
  --cpus-per-task 4 \
  --mem 8G \
  --batch-dry-run \
  --batch-module system/Miniconda3 \
  --batch-setup 'source /home/djones/.bashrc && eval "$(conda shell.bash hook)"' \
  --batch-condaenv ${PWD}/condaenv \
  --partition workq \
  --account djones \
  --time 4:00:00 \
  --job-name feature_counts \
  ${RESUME_PARAM}
