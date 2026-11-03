**************
Installation
**************

Requirements
============


Installation
============
To install CLPipe you currently need to build it from source.
Note: Do not use pip directly to install the required packages,
as not all dependencies can be installed this way. Instead, use conda, then pip to install clpipe.::


  git clone https://github.com/LSSTDESC/CLPipe.git 
  cd CLPipe
  conda env update -f txpipe_environment.yml
  conda activate txpipe_clp
  pip install .
  conda deactivate


  conda env update -f firecrown_environment.yml
  conda activate firecrown_clp
  pip install .
  conda env config vars set CSL_DIR=${CONDA_PREFIX}/cosmosis-standard-library
  conda deactivate
  conda activate firecrown_clp

  cd ${CONDA_PREFIX}
  source ${CONDA_PREFIX}/bin/cosmosis-configure
  cosmosis-build-standard-library main
  

To run the tests you can do::

  
  pytest
  
  
