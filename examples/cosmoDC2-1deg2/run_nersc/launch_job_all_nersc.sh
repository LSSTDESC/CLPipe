#!/usr/bin/bash
#SBATCH -A m1727
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

module load conda
export HDF5_DO_MPI_FILE_SYNC=0
conda activate /global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/conda_envs/txpipe_clp
# Use the CLPipe checkout this run directory lives in
export PYTHONPATH=$(cd ../../.. && pwd):/global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/TXPipe:$PYTHONPATH

ceci TXPipe.yml

conda deactivate
conda activate /global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/conda_envs/firecrown_developer_clp

ceci TJPCov.yml
ceci Firecrown.yml

cd ./outputs
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

srun -n ${SLURM_NTASKS} cosmosis --mpi sampler_file.ini
