#!/usr/bin/bash
$BATCH_ALL

module load conda
export HDF5_DO_MPI_FILE_SYNC=0
#export PYTHONPATH=/sps/lsst/groups/clusters/cl_pipeline_project/TXPipe:$PYTHONPATH
conda activate $CONDA_DIRECTORY/conda_envs/txpipe_clp
export PYTHONPATH=$CLPIPELINE_DIR:$PYTHONPATH
export PYTHONPATH=$TXPIPE_DIRECTORY:$PYTHONPATH

cd $COMPUTATIONDIR

ceci TXPipe.yml

conda deactivate
conda activate $CONDA_DIRECTORY/conda_envs/firecrown_developer_clp
export PYTHONPATH=$CLPIPELINE_DIR:$PYTHONPATH

ceci TJPCov.yml
ceci Firecrown.yml

cd $OUTPUTDIR
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mpirun -n ${SLURM_NTASKS} cosmosis --mpi cluster_counts_mean_mass_redshift_richness.ini

