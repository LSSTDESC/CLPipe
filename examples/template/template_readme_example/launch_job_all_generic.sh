#!/usr/bin/bash
$BATCH_ALL

$CONDA_SETUP
export HDF5_DO_MPI_FILE_SYNC=0
conda activate $TXPIPE_ENV
# Use the CLPipe checkout this run directory lives in
export PYTHONPATH=$(cd ../../.. && pwd):$TXPIPE_DIRECTORY:$PYTHONPATH

ceci TXPipe.yml

conda deactivate
conda activate $FIRECROWN_ENV

ceci TJPCov.yml
ceci Firecrown.yml

cd $OUTPUTDIR
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

$COSMOSIS_COMMAND
