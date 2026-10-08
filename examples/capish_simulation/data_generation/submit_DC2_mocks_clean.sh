#!/usr/bin/bash
#SBATCH --time=6:00:00
#SBATCH --partition=hpc,lsst
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --array=0-999%50
#SBATCH --job-name=DC2_mocks_clean
#SBATCH --output=logs/DC2_mocks_clean_%A_%a.out

module load conda
conda activate /sps/lsst/users/ebarroso/conda_envs/firecrown_clp

mkdir -p logs
mkdir -p mocks_clean


s1=$((SLURM_ARRAY_TASK_ID + 1))
s2=$((SLURM_ARRAY_TASK_ID * 100000 + 1))

python build_DC2_mocks_clean.py "$s1" "$s2"
