#!/usr/bin/env bash
#
set -euo pipefail

if [ -f "03a-stringtie_indiv.log" ]
then
    RESUME_PARAM="--batch-resume 03a-stringtie_indiv.log"
else
    RESUME_PARAM=""
fi

code/slurm_scripts/bin/pt "code/stringtie.sh {be} input/Sscl1980-nuclear.fasta input/Sscl1980-mRNA.gtf" output/crams/*.cram  \
| code/slurm_scripts/bin/sbatch_jobarray.sh \
  --cpus-per-task 4 \
  --mem 8G \
  --batch-dry-run \
  --batch-module system/Miniconda3 \
  --batch-setup 'source /home/djones/.bashrc && eval "$(conda shell.bash hook)"' \
  --batch-condaenv ${PWD}/condaenv \
  --partition work \
  --time 4:00:00 \
  --job-name stringtie \
  ${RESUME_PARAM}

