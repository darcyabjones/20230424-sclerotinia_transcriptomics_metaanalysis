#!/bin/bash --login

set -euo pipefail


if [ -f "07-gatk_haplotype_round2.log" ]
then
    RESUME_PARAM="--batch-resume ./07-gatk_haplotype_round2.log"
else
    RESUME_PARAM=""
fi

../code/slurm_scripts/bin/pt 'code/gatk_haplotype_caller.sh work/genome.fasta work/gatk_haplotype_round2/{be}.g.vcf.gz {}' output/crams_bqsr/*.cram \
| ../code/slurm_scripts/bin/sbatch_jobarray.sh \
  --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' \
  --batch-condaenv ${PWD}/condaenv \
  --batch-module system/Miniconda3 \
  --batch-dry-run \
  --cpus-per-task 8 \
  --account djones \
  --partition workq \
  --time 24:00:00 \
  --job-name haplotype_caller \
  ${RESUME_PARAM}
