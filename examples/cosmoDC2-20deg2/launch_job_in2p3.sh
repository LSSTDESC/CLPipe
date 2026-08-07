#!/usr/bin/bash
#SBATCH --time=1:00:00
#SBATCH --partition=hpc
#SBATCH --ntasks=30
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb

module load conda
#source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/txpipe_clp
export HDF5_DO_MPI_FILE_SYNC=0
export PYTHONPATH=/sps/lsst/groups/clusters/cl_pipeline_project/TXPipe:$PYTHONPATH
export PYTHONPATH=../../:$PYTHONPATH
ceci TXPipe.yml
conda deactivate
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/firecrown_developer_clp
export PYTHONPATH=../../:$PYTHONPATH
ceci TJPCov.yml

ceci Firecrown.yml
cd outputs
cosmosis sampler_file.ini
