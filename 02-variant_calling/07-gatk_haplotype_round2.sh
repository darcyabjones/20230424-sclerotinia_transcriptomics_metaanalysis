#!/bin/bash --login

set -euo pipefail

source ./file_locations.sh

if [ -f "07-gatk_haplotype_round2.log" ]
then
    RESUME_PARAM="--batch-resume ./07-gatk_haplotype_round2.log"
else
    RESUME_PARAM=""
fi

../code/slurm_scripts/bin/pt 'code/gatk_haplotype_caller.sh input/Sscl1980-nuclear.fasta work/gatk_haplotype_round2/{be}.g.vcf.gz {}' output/crams_bqsr/*.cram \
| ../code/slurm_scripts/bin/sbatch_jobarray.sh \
  --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' \
  --batch-condaenv ${PWD}/condaenv \
  --batch-module system/Miniconda3 \
  --batch-dry-run \
  --batch-max-simultaneous 8 \
  --cpus-per-task 4 \
  --account djones \
  --partition workq \
  --time 6:00:00 \
  --job-name haplotype_caller \
  ${RESUME_PARAM}
