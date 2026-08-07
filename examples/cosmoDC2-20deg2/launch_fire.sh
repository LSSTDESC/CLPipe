#!/usr/bin/bash
#SBATCH --time=00:10:00
#SBATCH --partition=hpc,lsst
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb

module load conda
export PYTHONPATH=../../:$PYTHONPATH
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/firecrown_developer_clp
export PYTHONPATH=../../:$PYTHONPATH
ceci CL_cosmoDC2-20deg2_concat_in2p3.yml --yamlId Firecrown
#cd outputs
#cosmosis sampler_file.ini
