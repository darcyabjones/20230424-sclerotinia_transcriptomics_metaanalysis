#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --export=NONE
#SBATCH --time=4:00:00
#SBATCH --job-name=gfeats
#SBATCH --account=djones

mkdir -p work

module load system/Miniconda3
#source /home/djones/.bashrc
unset PERL5LIB
eval "$(conda shell.bash hook)"
conda activate ./condaenv

set -euo pipefail

#code/tabix_gff3.sh output/pdb.gff3 output/pdb.gff3.gz
#rm -f output/pdb.gff3

#code/tabix_gff3.sh output/swissprot.gff3 output/swissprot.gff3.gz
#rm -f output/swissprot.gff3
#
#gunzip --stdout output/swissprot.gff3.gz > output/swissprot.gff3
#predutils prot_to_genome -o output/swissprot-genome.gff3 work/Sscl1980-mRNA.gff3 output/swissprot.gff3
#rm -f output/swissprot.gff3
#
#gunzip --stdout output/pdb.gff3.gz > output/pdb.gff3
#predutils prot_to_genome -o output/pdb-genome.gff3 work/Sscl1980-mRNA.gff3 output/pdb.gff3
#rm -f output/pdb.gff3

code/tabix_gff3.sh output/pdb-genome.gff3 output/pdb-genome.gff3.gz
code/tabix_gff3.sh output/swissprot-genome.gff3 output/swissprot-genome.gff3.gz
