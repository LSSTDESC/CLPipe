#!/usr/bin/bash
#SBATCH --time=15:00:00
#SBATCH --partition=hpc,lsst
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --ntasks=30 
module load conda
export HDF5_DO_MPI_FILE_SYNC=0
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/txpipe_clp
export PYTHONPATH=../../../:$PYTHONPATH
ceci CL_cosmoDC2-full_concat_mor.yml --yamlId TJPCov
conda deactivate
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/firecrown_developer_clp
export PYTHONPATH=../../../:$PYTHONPATH
ceci CL_cosmoDC2-full_concat_mor.yml --yamlId Firecrown
cd outputs_mor
cosmosis sampler_file.ini
