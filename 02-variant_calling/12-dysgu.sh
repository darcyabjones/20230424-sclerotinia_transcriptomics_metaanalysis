#!/bin/bash --login
module load system/Miniconda3
source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"
conda activate ./condaenv/

set -euo pipefail

REFERENCE_GENOME="$(realpath "work/genome.fasta")"

WORK="${PWD}/work/dysgu_hap"
mkdir -p "${WORK}"

shopt -s extglob
export GLOBIGNORE="SRR14280050.cram"

COMMANDS=$(../code/slurm_scripts/bin/pt "dysgu run --clean --procs 4 --min-support 4 -v 2 -o '${WORK}/{be}.vcf' '${REFERENCE_GENOME}' '${WORK}/tmp_{be}' '{}'" output/crams/*.cram)
#JOBID=$(echo "${COMMANDS}" | ../code/slurm_scripts/bin/sbatch_jobarray.sh --batch-resume ./12-dysgu.log --batch-module system/Miniconda3 --batch-condaenv "${PWD}/condaenv" --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' --cpus-per-task 4 --account djones --partition workq --time 1-00:00:00 --job-name dysgu_hap)
SITES="${PWD}/work/dysgu_sites.vcf"

read -d '' COMMAND <<EOF || :
#!/bin/bash --login

module load system/Miniconda3
source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"
conda activate "${PWD}/condaenv"

cd work/dysgu_hap/
dysgu merge -o '${SITES}' *.vcf
EOF

if [ ! -z "${JOBID:-}" ]
then
  DEPSTRING="--dependency=afterok:${JOBID}"
else
  DEPSTRING=""
fi

unset JOBID

#JOBID=$(echo "${COMMAND}" | sbatch --mem 32G --partition workq --parsable --job-name dysgu_sites --account djones ${DEPSTRING})


WORK="${PWD}/work/dysgu_geno"
mkdir -p "${WORK}"

if [ ! -z "${JOBID:-}" ]
then
  DEPSTRING="--dependency=afterok:${JOBID}"
else
  DEPSTRING=""
fi
COMMANDS=$(../code/slurm_scripts/bin/pt "dysgu run --clean --procs 4 --min-support 4 -v 2 -o '${WORK}/{be}.vcf' --sites '${SITES}' '${REFERENCE_GENOME}' '${WORK}/tmp_{be}' '{}'" output/crams/*.cram)
JOBID=$(echo "${COMMANDS}" | ../code/slurm_scripts/bin/sbatch_jobarray.sh ${DEPSTRING} --batch-module system/Miniconda3 --batch-condaenv "${PWD}/condaenv" --batch-setup 'source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"' --cpus-per-task 4 --account djones --partition workq --time 3:00:00 --job-name dysgu_geno)
exit 0
read -d '' COMMAND <<EOF || :
#!/bin/bash --login

module load system/Miniconda3
source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"
conda activate "${PWD}/condaenv"

dysgu merge -o work/dysgu_merged_genos.vcf  work/dysgu_geno/*.vcf
bcftools filter -Ou --include 'COUNT(GT="AA") > 0' work/dysgu_merged_genos.vcf \
| bcftools view -f 'PASS,.' -O z -o output/dysgu.vcf.gz
bcftools index --tbi output/dysgu.vcf.gz
EOF

JOBID=$(echo "${COMMAND}" | sbatch --mem 8G --partition workq --dependency=afterok:${JOBID} --parsable --job-name dysgu_filter --account djones)

