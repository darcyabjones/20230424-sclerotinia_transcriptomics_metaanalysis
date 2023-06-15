#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --time=24:00:00
#SBATCH --job-name=deeploc2

module load system/Python-3.11.1

source deeploc2_env/bin/activate

set -euo pipefail

mkdir -p work/deeploc2_results
export TORCH_HOME="${PWD}/work/torch"
srun -n 1 -N 1 -c 16 --export=ALL deeploc2 --fasta work/Sscl1980.faa --output work/deeploc2_results --model Fast
