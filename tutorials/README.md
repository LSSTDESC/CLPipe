# CLPipe tutorials

CLPipe runs DESC tools one after the other to go from cluster and galaxy
catalogs to cosmological constraints. For now, the analysis uses binned cluster counts and stacked weak
lensing profiles. These tutorials show what each step does and how to run it.
Installation and the example analyses are described in the main
[README](https://github.com/LSSTDESC/CLPipe#readme).

## The pipeline

```
cluster and galaxy catalogs
        |   TXPipe cluster stages
        v
SACC file with the counts and the stacked DeltaSigma profiles
        |   CLPCovariance (TJPCov)
        v
SACC file with the covariance
        |   CLPFirecrown
        v
Firecrown likelihood and CosmoSIS files
        |   CosmoSIS, run as a separate job
        v
chains
```

| Step | Code | Tutorial |
|---|---|---|
| Measurements | cluster stages of TXPipe (in the TXPipe repository) | [txpipe.md](txpipe.md) |
| Covariance | `clpipe/clp_covariance.py`, built on TJPCov | [covariance_tjpcov.ipynb](covariance_tjpcov.ipynb) |
| Likelihood and sampling | `clpipe/clp_firecrown.py`, with Firecrown and CosmoSIS | [likelihood_firecrown.ipynb](likelihood_firecrown.ipynb) |
| Predictions for given parameters | `clpipe/cluster_theory_pred.py` | [theory_prediction.ipynb](theory_prediction.ipynb) |
| Plots of the chains | `clpipe/utils/cosmosis_mcmc_plots.py` | used in `examples/cosmodc2_redmapper/plot_samples.ipynb` |

## Where to start

If you are new to CLPipe, start with
[running_an_analysis.md](running_an_analysis.md). It explains how a run is
configured and launched, and where its outputs end up. Then go through the
steps in the order of the pipeline.

To add your own step to the pipeline, see
[writing_a_stage.ipynb](writing_a_stage.ipynb).

## Environments and kernels

The notebooks in this folder run in the Firecrown environment of the pipeline.
The main README describes the two environments, where the shared ones are at
CC-IN2P3 and NERSC, and how to create them on your own machine.

To run the notebooks on a Jupyter platform, copy the shared kernel to your home
directory once, then select the `firecrown_clp` kernel in Jupyter. Despite its
name, this kernel starts the shared `firecrown_developer_clp` environment.

At CC-IN2P3 (https://notebook.cc.in2p3.fr):

```bash
mkdir -p ${HOME}/.local/share/jupyter/kernels
cp -r /sps/lsst/groups/clusters/cl_pipeline_project/jupyter_kernel/firecrown_clp ${HOME}/.local/share/jupyter/kernels/
```

At NERSC (https://jupyter.nersc.gov):

```bash
mkdir -p ${HOME}/.local/share/jupyter/kernels
cp -r /global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/jupyter_kernels/firecrown_clp ${HOME}/.local/share/jupyter/kernels/
```

The `txpipe_clp` kernel is in the same folders, if you want to work with the
TXPipe environment.

On your own machine, register the Firecrown environment as a kernel with
`python -m ipykernel install --user --name firecrown_clp`.

## Files in this folder

The notebooks are saved with their outputs, so you can read them on GitHub or
in the documentation without running them. When you run them, they read their
inputs from `data/` and write their files to `outputs/`, which git ignores.

The files in `tutorials/data/` are small and tracked by git. The large files of
the example analyses are in the top-level `data/` folder instead, downloaded
from the NERSC portal (see the main README).
