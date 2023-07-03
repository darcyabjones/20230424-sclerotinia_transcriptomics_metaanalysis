#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --export=NONE
#SBATCH --time=24:00:00
#SBATCH --job-name=index

mkdir -p work

awk -F'\t' -v OFS="\t" '$12 !~ /^[[:space:]]*$/ {print $1, $12}' input/effectordb.tsv \
| tail -n+2 \
> work/effectordb_to_split.tsv

for i in {0..29}
do
  awk -F'\t' -v I=${i} 'NR % 30 == I {printf ">%s\n%s\n", $1, $2}' work/effectordb_to_split.tsv > "work/effectordb_${i}.fasta"
done


#wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
#mv uniref50.fasta.gz work
#gunzip work/uniref50.fasta.gz
