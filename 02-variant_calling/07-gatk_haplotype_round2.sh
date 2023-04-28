#!/bin/bash --login

set -euo pipefail

source ./file_locations.sh

if [ -f "07-gatk_r2_haplotype_caller.log" ]
then
    RESUME_PARAM="--batch-resume ./07-gatk_r2_haplotype_caller.log"
else
    RESUME_PARAM=""
fi

code/slurm_scripts/bin/pt 'code/haplotype_caller.sh input/ArME14.fasta work/gatk_r2_haplotype_caller/{be}.g.vcf.gz {}' output/recalibrated_illumina_alignments/*.cram \
  | code/slurm_scripts/bin/sbatch_jobarray.sh --batch-module miniconda3/latest --cpus-per-task 2 --partition work --time 3:00:00 --job-name haplotype_caller ${RESUME_PARAM}
