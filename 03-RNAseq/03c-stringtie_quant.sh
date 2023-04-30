#!/usr/bin/env bash
#
set -euo pipefail

if [ -f "03c-stringtie_quant.log" ]
then
    RESUME_PARAM="--batch-resume 03c-stringtie_quant.log"
else
    RESUME_PARAM=""
fi

tail -n+2 input/sra_rnaseq.tsv \
| awk -F '\t' '$14 == "TRUE"' \
| ../code/slurm_scripts/bin/pt --file - "code/stringtie_quant.sh {0} {7} input/Sscl1980-nuclear.fasta output/stringtie_merged.gtf"  \
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
  --job-name stringtie_quant \
  ${RESUME_PARAM}

