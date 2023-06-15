#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --export=NONE
#SBATCH --time=4:00:00
#SBATCH --job-name=gfeats
#SBATCH --account=djones

mkdir -p work/fungismash

#cp output/fungismash.zip work/fungismash
#cd work/fungismash
#unzip fungismash.zip
#cd -

module load system/Python-3.11.1

source gffpal_env/bin/activate

gffpal antismash2gff -g work/Sscl1980-mRNA.gff3 -o output/fungismash.gff3 work/fungismash/chr*.gbk

deactivate


module load system/Miniconda3
unset PERL5LIB
eval "$(conda shell.bash hook)"
conda activate ./condaenv

code/tabix_gff3.sh output/fungismash.gff3 output/fungismash.gff3.gz
