#!/usr/bin/env bash
#
set -euo pipefail

if [ -f "02-align_sra.log" ]
then
    RESUME_PARAM="--batch-resume 02-align_sra.log"
else
    RESUME_PARAM=""
fi

code/slurm_scripts/bin/pt "code/align_sra.sh {} input/Sscl1980-nuclear.fasta input/Sscl1980-mRNA.gtf input/Sscl1980-cds.fna work/hisat2_index work/kallisto_index" $(tail -n+2 input/sra_rnaseq.tsv | cut -f 1 | sort -u)  \
| code/slurm_scripts/bin/sbatch_jobarray.sh \
  --cpus-per-task 4 \
  --mem 8G \
  --batch-dry-run \
  --batch-max-simultaneous 8 \
  --account djones \
  --batch-module system/Miniconda3 \
  --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' \
  --batch-condaenv ${PWD}/condaenv \
  --partition workq \
  --time 4:00:00 \
  --job-name align_sra \
  ${RESUME_PARAM}

