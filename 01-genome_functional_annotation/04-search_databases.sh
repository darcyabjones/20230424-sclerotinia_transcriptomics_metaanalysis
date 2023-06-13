#!/usr/bin/env bash

set -euo pipefail


#sbatch -n 3 -N 1 --time 04:00:00 --account djones --partition workq <<EOF
##!/usr/bin/bash --login
#
#wget -c -O work/pdb.faa.gz https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz &
#wget -c -O work/swissprot.faa.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz &
##wget -c -O work/alphafold.faa https://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta &
#wait
#gunzip work/pdb.faa.gz
#gunzip work/swissprot.faa.gz
#EOF

#if [ -f "04b-search_databases_index.log" ]
#then
#    RESUME_PARAM="--batch-resume 04b-search_databases_index.log"
#else
#    RESUME_PARAM=""
#fi
#
#../code/slurm_scripts/bin/pt 'mkdir -p work/{be}_db; mmseqs createdb {} work/{be}_db/db; cp {} work/{be}_db/db.fasta' work/{Sscl1980,pdb,swissprot}.faa work/foldseek_results_all_seqs.fasta \
#| ../code/slurm_scripts/bin/sbatch_jobarray.sh \
#  --cpus-per-task 1 \
#  --account djones \
#  ${RESUME_PARAM} \
#  --batch-module system/Miniconda3 \
#  --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' \
#  --batch-condaenv ${PWD}/condaenv \
#  --partition workq \
#  --time 1:00:00 \
#  --job-name index_mmseqs2
#
#exit 0

if [ -f "04c-search_databases_search.log" ]
then
    RESUME_PARAM="--batch-resume 04c-search_databases_search.log"
else
    RESUME_PARAM=""
fi

../code/slurm_scripts/bin/pt 'code/run_mmseqs_search.sh work/Sscl1980_db/db {} work/{dbr/_db/}-matches.tsv' work/{pdb_db,swissprot_db}/db work/foldseek_results_all_seqs_db/db \
| ../code/slurm_scripts/bin/sbatch_jobarray.sh \
  --cpus-per-task 16 \
  --account djones \
  --batch-dry-run \
  ${RESUME_PARAM} \
  --batch-module system/Miniconda3 \
  --batch-module bioinfo/mmseqs2-v14-7e284 \
  --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' \
  --batch-condaenv ${PWD}/condaenv \
  --partition workq \
  --time 2:00:00 \
  --job-name search_mmseqs2
