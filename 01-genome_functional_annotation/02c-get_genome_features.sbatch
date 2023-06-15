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

echo "here"
module load system/Miniconda3
echo "here"
#source /home/djones/.bashrc
unset PERL5LIB
eval "$(conda shell.bash hook)"
echo "here"
conda activate ./condaenv

echo "here"
set -euo pipefail

#cp input/Sscl1980-mRNA.gff3.gz work
#gunzip --force work/Sscl1980-mRNA.gff3.gz

#predutils ipr_to_gff -o output/interproscan-protein.gff3 work/ips.xml
#cp work/ips.tsv output/interproscan-protein.tsv
#cp work/ips.xml output/interproscan-protein.xml

#mkdir output/interproscan_genome
#predutils prot_to_genome -o output/interproscan_genome/ --split work/Sscl1980-mRNA.gff3 output/interproscan-protein.gff3

#for f in output/interproscan_genome/*.gff3
#do
#  code/tabix_gff3.sh "output/interproscan_genome/$(basename $f)" "output/interproscan_genome/$(basename $f).gz"
#done

#mkdir -p output/predector_{protein,genome}
#predutils r2js tmbed --software-version 1.0.0 output/tmbed.3line work/Sscl1980.faa | predutils gff - | awk -F '\t' -v OFS="\t" '{$1=gensub(/^(\S+).*/, "\\1", "g", $1); print}' > output/predector_protein/tmbed:1.0.0.gff3

#predutils prot_to_genome -o "output/predector_genome/tmbed:1.0.0.gff3" work/Sscl1980-mRNA.gff3 "output/predector_protein/tmbed:1.0.0.gff3"
#code/tabix_gff3.sh "output/predector_genome/tmbed:1.0.0.gff3" "output/predector_genome/tmbed:1.0.0.gff3.gz"


#cp -r work/predector_results/ output/predector_protein

predutils map_to_genome --split "source" -o "output/predector_genome/" work/Sscl1980-mRNA.gff3 output/predector_protein/Scscl_1980-proteins.gff3
predutils scores_to_genome --outfile output/predector_genome/Sscl1980-proteins.bedgraph  work/Sscl1980-mRNA.gff3 output/predector_protein/Scscl_1980-proteins-ranked.tsv

for f in output/predector_genome/*.gff3
do
  code/tabix_gff3.sh "$f" "${f}.gz"
done
code/multibedgraph_to_bigwig.sh input/Sscl1980-nuclear.fasta.gz.fai output/predector_genome/Sscl1980-proteins.bedgraph "output/predector_genome/"
