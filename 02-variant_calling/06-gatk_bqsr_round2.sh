#!/bin/bash --login

set -euo pipefail

source ./file_locations.sh

if [ -f "06-gatk_r2_bqsr.log" ]
then
    RESUME_PARAM="--batch-resume ./06-gatk_r2_bqsr.log"
else
    RESUME_PARAM=""
fi

code/slurm_scripts/bin/pt 'code/bqsr.sh input/ArME14.fasta output/recalibrated_illumina_alignments/{b} output/gatk_r1_joint_genotypes-filtered.g.vcf.gz {}' output/illumina_alignments/*.cram \
  | code/slurm_scripts/bin/sbatch_jobarray.sh --batch-module miniconda3/latest --cpus-per-task 1 --mem 4G --partition work --time 3:00:00 --job-name bqsr ${RESUME_PARAM}
