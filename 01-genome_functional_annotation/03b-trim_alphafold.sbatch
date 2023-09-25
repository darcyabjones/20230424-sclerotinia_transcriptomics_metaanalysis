#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --export=NONE
#SBATCH --time=12:00:00
#SBATCH --job-name=trim

unset PERL5LIB
module load system/Miniconda3
source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
conda activate ./condaenv

set -euo pipefail

export OMP_NUM_CPUS=4

if [ -s 03b-trim_alphafold.log ]
then
  RESUME_ARG="--resume-failed"
else
  RESUME_ARG=""
fi

CMD="srun -N 1 -n 1 --cpus-per-task '${SLURM_CPUS_PER_TASK:-1}' --exact --export=ALL code/trim_alphafold_cifs.py --compress --outdir work/alphafold_structures_trimmed -w 5 --targetp ./envs/targetp-2.0/bin/targetp -"
find work/alphafold_structures/ -name "*.cif.gz" | parallel --delay 0.5 -j "${SLURM_NTASKS:-1}" --pipe -N1000 --joblog 03b-trim_alphafold.log ${RESUME_ARG} "${CMD}" 
