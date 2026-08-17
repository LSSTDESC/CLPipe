#!/usr/bin/bash
#SBATCH --time=20:00:00
#SBATCH --partition=hpc,lsst
#SBATCH --cpus-per-task=1
#SBATCH --mem=80gb
#SBATCH --ntasks=51

module load conda
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/firecrown_developer_clp
export PYTHONPATH=/sps/lsst/users/ebarroso/CLPipe:$PYTHONPATH

cd ./outputs_both
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mpirun -n ${SLURM_NTASKS} cosmosis --mpi sampler_file.ini
