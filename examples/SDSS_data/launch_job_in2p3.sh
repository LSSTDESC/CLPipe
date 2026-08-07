#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --partition=hpc,lsst
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --ntasks=1



module load conda
#source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/firecrown_clp
export PYTHONPATH=../../:$PYTHONPATH
ceci CL_SDSS_concat.yml --yamlId Firecrown
cd outputs
cosmosis sampler_file.ini
