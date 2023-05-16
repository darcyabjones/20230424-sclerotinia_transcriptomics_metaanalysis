#!/bin/bash --login

set -euo pipefail

if [ -f "06-gatk_bqsr_round2.log" ]
then
    RESUME_PARAM="--batch-resume ./06-gatk_bqsr_round2.log"
else
    RESUME_PARAM=""
fi

../code/slurm_scripts/bin/pt 'code/gatk_bqsr.sh work/genome.fasta output/crams_bqsr/{be} output/round1/filtered.g.vcf.gz {}' output/crams/*.cram \
| ../code/slurm_scripts/bin/sbatch_jobarray.sh \
  --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' \
  --batch-condaenv ${PWD}/condaenv \
  --batch-module system/Miniconda3 \
  --batch-dry-run \
  --cpus-per-task 1 \
  --batch-max-simultaneous 8 \
  --mem 4G \
  --partition workq \
  --account djones \
  --time 3:00:00 \
  --job-name bqsr \
  ${RESUME_PARAM}
