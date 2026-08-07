#!/bin/bash
#SBATCH -A m1727
#SBATCH -C cpu
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

module load conda
conda activate /global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/conda_envs/txpipe_clp
export HDF5_DO_MPI_FILE_SYNC=0
export PYTHONPATH=/global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/TXPipe:$PYTHONPATH
export PYTHONPATH=../../:$PYTHONPATH
ceci CL_cosmoDC2-20deg2_concat_nersc.yml --yamlId TXPipe 
ceci CL_cosmoDC2-20deg2_concat_nersc.yml --yamlId TJPCov
conda deactivate
conda activate /global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/conda_envs/firecrown_clp
export PYTHONPATH=../../:$PYTHONPATH
ceci CL_cosmoDC2-20deg2_concat_nersc.yml --yamlId Firecrown
cd outputs
cosmosis sampler_file.ini
