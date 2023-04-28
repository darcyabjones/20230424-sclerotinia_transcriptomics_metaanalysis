#!/bin/bash --login

set -euo pipefail

if [ -f "02a-align_sra.log" ]
then
    RESUME_PARAM="--batch-resume ./02a-align_sra.log"
else
    RESUME_PARAM=""
fi
IDS=( $(tail -n+2 input/resequencing.tsv | awk -F'\t' '!($1 ~ /^\s*$/) {print $1}') )

../code/slurm_scripts/bin/pt 'code/align_sra.sh {} work/genome.fasta output' "${IDS[@]}" \
| ../code/slurm_scripts/bin/sbatch_jobarray.sh \
  --cpus-per-task 4 \
  --mem 8G \
  --batch-dry-run \
  --account djones \
  --batch-module system/Miniconda3 \
  --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' \
  --batch-condaenv ${PWD}/condaenv \
  --partition workq \
  --time 4:00:00 \
  --job-name align_sra \
  ${RESUME_PARAM}
