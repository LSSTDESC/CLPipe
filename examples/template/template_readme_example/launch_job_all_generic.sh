#!/usr/bin/bash
$BATCH_ALL

module load conda
conda activate $CONDA_DIRECTORY/conda_envs/txpipe_clp
export HDF5_DO_MPI_FILE_SYNC=0
export PYTHONPATH=$CONDA_DIRECTORY/TXPipe:$PYTHONPATH
export PYTHONPATH=../../:$PYTHONPATH

cd $COMPUTINGDIR
ceci TXPipe.yml
ceci TJPCov.yml
conda deactivate
conda activate $CONDA_DIRECTORY/conda_envs/firecrown_developer_clp
export PYTHONPATH=../../:$PYTHONPATH
ceci Firecrown.yml
cd $OUTPUTDIR
cosmosis sampler_file.ini
