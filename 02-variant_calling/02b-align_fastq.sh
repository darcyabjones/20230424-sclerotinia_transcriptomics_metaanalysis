#!/bin/bash --login

set -euo pipefail

if [ -f "02b-align_fastq.log" ]
then
    RESUME_PARAM="--batch-resume ./02b-align_fastq.log"
else
    RESUME_PARAM=""
fi

read -r -d '' TMP <<'EOF' || :
code/align_fastq.sh SSCLLC41 work/genome.fasta input/T20210808_00001_RUN305_1.fq.gz input/T20210808_00001_RUN305_2.fq.gz output
code/align_fastq.sh SSCLLANG work/genome.fasta input/L_angustifoli_sclerotinia_HGC7MBCXX_GCCAATAT_L001_R1.fastq.gz input/L_angustifoli_sclerotinia_HGC7MBCXX_GCCAATAT_L001_R2.fastq.gz output
code/align_fastq.sh SSCLLMUT work/genome.fasta input/L_mutab_Sclerotinia_HGC7MBCXX_CGATGTAT_L001_R1.fastq.gz input/L_mutab_Sclerotinia_HGC7MBCXX_CGATGTAT_L001_R2.fastq.gz output
EOF

echo "${TMP}" | ../code/slurm_scripts/bin/sbatch_jobarray.sh \
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
  --job-name align_fastq \
  ${RESUME_PARAM}
