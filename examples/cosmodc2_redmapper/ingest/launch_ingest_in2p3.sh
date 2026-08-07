#!/usr/bin/bash
#SBATCH --time=7:00:00
#SBATCH --partition=hpc,lsst
#SBATCH --cpus-per-task=1
#SBATCH --mem=300gb
#SBATCH --ntasks=1

module load conda
export HDF5_DO_MPI_FILE_SYNC=0
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/txpipe_clp
export PYTHONPATH=/sps/lsst/users/ebarroso/CLPipe:$PYTHONPATH
export PYTHONPATH=/sps/lsst/groups/clusters/cl_pipeline_project/TXPipe/:$PYTHONPATH

ceci TXPipe_ingest.yml
