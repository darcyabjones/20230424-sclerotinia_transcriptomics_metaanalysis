#!/bin/bash --login

unset PERL5LIB
module load system/Miniconda3
source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
conda activate ./caretta_condaenv

set -euo pipefail

mkdir -p work
mkdir -p output/08c-hmmer_search_clusters

# Download uniref100 FUNGI
#sbatch -N 1 -n 1 --time 12:00:00 --account=djones --export=ALL <<EOF
##!/bin/bash --login
##wget -c -O work/uniref.fasta.gz "https://rest.uniprot.org/uniref/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A4751%29%29+AND+%28identity%3A0.9%29"
#
#wget -c -O "work/uniref.tsv.gz" "https://rest.uniprot.org/uniref/stream?compressed=true&fields=id%2Cname%2Ctypes%2Ccount%2Clength%2Cidentity%2Corganism_id&format=tsv&query=%28%28taxonomy_id%3A4751%29%29+AND+%28identity%3A0.9%29"
#
#gunzip work/uniref.fasta.gz
#EOF


#if [ -s 08c-hmmer_search_clusters_jackhmmer.log ]
#then
#  RESUME_COMMAND="--batch-resume 08c-hmmer_search_clusters_jackhmmer.log"
#else
#  RESUME_COMMAND=""
#fi
#
#pt 'code/hmmsearch_iterative.sh {} work/uniref.fasta output/08c-hmmer_search_clusters' output/08b-aligned_clusters/*/*.fasta \
#| ../code/slurm_scripts/bin/sbatch_jobarray.sh \
#  --account djones \
#  --export=ALL \
#  --ntasks 1 \
#  --cpus-per-task 2 \
#  --mem 32G \
#  --partition workq \
#  --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' \
#  --batch-module system/Miniconda3 \
#  ${RESUME_COMMAND} \
#  --time 6:00:00 \
#  --batch-condaenv ${PWD}/caretta_condaenv 
#
#exit 0



sbatch -N 1 -n 1 --mem 16G --time 12:00:00 --account=djones --export=ALL <<'EOF'
#!/bin/bash --login

set -euo pipefail

#zcat work/uniref.tsv.gz \
#| gawk -F'\t' -v OFS="\t" '!/^Cluster/ {split($7, t, "[[:space:]]*;[[:space:]]*"); for (ti in t) {print $1, t[ti]}}' \
#| sort -u \
#> work/uniref_taxids.tsv

gawk -v OFS='\t' -v FULL=0.00001 -v DOMAIN=0.001 '(!/^#/) && ($5 <= FULL) && ($8 <= DOMAIN) {print $1, $3}' output/08c-hmmer_search_clusters/*.tblout \
| sort -u \
> work/hmmer_matches.tblout

gawk -F '\t' -v OFS="\t" 'NR==FNR {if ($1 in TIDS) {TIDS[$1]=TIDS[$1]";"$2} else {TIDS[$1]=$2}; next} $1 in TIDS {split(TIDS[$1], tids, ";"); for (ti in tids) {print $1, $2, tids[ti]}}' work/uniref_taxids.tsv work/hmmer_matches.tblout \
| sort -u \
> output/08c-hmmer_search_clusters_taxids.tsv
EOF
