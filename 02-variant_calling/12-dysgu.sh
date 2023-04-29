#!/bin/bash --login
module load miniconda3/latest
conda activate ./condaenv/

set -euo pipefail

REFERENCE_GENOME="$(realpath "input/Sscl1980-nuclear.fasta")"

WORK="${PWD}/work/dysgu_hap"
mkdir -p "${WORK}"

COMMANDS=$(../code/slurm_scripts/bin/pt "dysgu run --clean --procs 4 --min-support 4 -v 2 -o '${WORK}/{be}.vcf' '${REFERENCE_GENOME}' '${WORK}/tmp_{be}' '{}'" output/crams/*.cram)
#JOBID=$(echo "${COMMANDS}" | code/slurm_scripts/bin/sbatch_jobarray.sh --batch-dry-run --batch-module miniconda3/latest --batch-condaenv "${PWD}/condaenv" --cpus-per-task 4 --account djones --partition workq --time 3:00:00 --job-name dysgu_hap)
SITES="${PWD}/work/dysgu_sites.vcf"

read -d '' COMMAND <<EOF || :
#!/bin/bash --login

module load miniconda3/latest
conda activate "${PWD}/condaenv"

cd work/dysgu_hap
dysgu merge -o '${SITES}' *.vcf
EOF

JOBID=$(echo "${COMMAND}" | sbatch --mem 8G --partition workq --parsable --job-name dysgu_sites --account djones) #  --dependency=afterok:${JOBID}

WORK="${PWD}/work/dysgu_geno"
mkdir -p "${WORK}"

COMMANDS=$(../code/slurm_scripts/bin/pt "dysgu run --clean --procs 4 --min-support 4 -v 2 -o '${WORK}/{be}.vcf' --sites '${SITES}' '${REFERENCE_GENOME}' '${WORK}/tmp_{be}' '{}'" output/crams/*.cram)
JOBID=$(echo "${COMMANDS}" | ../code/slurm_scripts/bin/sbatch_jobarray.sh --dependency=afterok:${JOBID} --batch-module miniconda3/latest --batch-condaenv "${PWD}/condaenv" --cpus-per-task 4 --account djones --partition workq --time 3:00:00 --job-name dysgu_geno)

read -d '' COMMAND <<EOF || :
#!/bin/bash --login

module load miniconda3/latest
conda activate "${PWD}/condaenv"

dysgu merge -o work/dysgu_merged_genos.vcf  work/dysgu_geno/*.vcf
bcftools filter -Ou --include 'COUNT(GT="AA") > 0' work/dysgu_merged_genos.vcf \
| bcftools view -f 'PASS,.' -O z -o output/dysgu.vcf.gz
bcftools index --tbi output/dysgu.vcf.gz
EOF

JOBID=$(echo "${COMMAND}" | sbatch --mem 8G --partition workq --dependency=afterok:${JOBID} --parsable --job-name dysgu_filter --account djones)

