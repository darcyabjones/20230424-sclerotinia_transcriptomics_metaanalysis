#!/bin/bash --login

module load system/Miniconda3
source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"

conda activate ./condaenv

set -euo pipefail

if [ -f "03-gatk_haplotype_round1.log" ]
then
    RESUME_PARAM="--batch-resume ./03-gatk_haplotype_round1.log"
else
    RESUME_PARAM=""
fi

if [ ! -f "work/genome.fasta.dict" ]
then
  picard CreateSequenceDictionary \
        -R work/genome.fasta \
        -O work/genome.fasta.dict
fi
cp work/genome.fasta.dict work/genome.dict

if [ ! -f "work/genome.fasta.fai" ]
then
  samtools faidx work/genome.fasta
fi

../code/slurm_scripts/bin/pt 'code/gatk_haplotype_caller.sh work/genome.fasta work/gatk_haplotype_round1/{be}.g.vcf.gz {}' output/crams/*.cram \
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

