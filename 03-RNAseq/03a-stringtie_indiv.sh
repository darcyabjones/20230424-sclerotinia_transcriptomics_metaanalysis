#!/usr/bin/env bash
#
set -euo pipefail

if [ -f "03a-stringtie_indiv.log" ]
then
    RESUME_PARAM="--batch-resume 03a-stringtie_indiv.log"
else
    RESUME_PARAM=""
fi

tail -n+2 input/sra_rnaseq.tsv \
| awk -F '\t' '$8 != "unstranded" && $14 == "TRUE"' \
| ../code/slurm_scripts/bin/pt --file - --group '{14}' \
  'code/stringtie.sh {0@j/_/} {7@uj/_/} input/Sscl1980-nuclear.fasta work/genes.gtf {0@s~^~output/crams/~s/$/.cram/}' \
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
  --job-name stringtie \
  ${RESUME_PARAM}

