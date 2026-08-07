#!/usr/bin/bash
#SBATCH --time=25:00:00
#SBATCH --partition=hpc,lsst
#SBATCH --cpus-per-task=1
#SBATCH --mem=80gb
#SBATCH --ntasks=30

module load conda
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/firecrown_developer_clp
export PYTHONPATH=/sps/lsst/users/ebarroso/CLPipe:$PYTHONPATH

ceci TJPCov.yml
ceci Firecrown.yml

cd ./outputs_both
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mpirun -n ${SLURM_NTASKS} cosmosis --mpi cluster_counts_mean_mass_redshift_richness.ini

