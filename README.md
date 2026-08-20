# CLPipe

Repository dedicated to the Cluster Working Group of the DESC-LSST collaboration.

## Installation

To run the examples, there is no need to create new Conda environments — shared
environments are already available on both IN2P3 and NERSC.

Two environments are provided:

- **Firecrown** — for MCMC inference
- **TXPipe and TJPCov** — for data computation and covariance

### Firecrown

```bash
# IN2P3
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/firecrown_clp

# NERSC
conda activate /global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/conda_envs/firecrown_clp
```

### TXPipe and TJPCov

```bash
# IN2P3
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/txpipe_clp

# NERSC
conda activate /global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/conda_envs/txpipe_clp
```

### Local Installation

To reproduce these environments locally:

```bash
conda env update -f txpipe_environment.yml
conda env update -f firecrown_environment.yml

conda activate firecrown_clp
conda env config vars set CSL_DIR=${CONDA_PREFIX}/cosmosis-standard-library
conda deactivate
conda activate firecrown_clp

cd ${CONDA_PREFIX}
source ${CONDA_PREFIX}/bin/cosmosis-configure
cosmosis-build-standard-library main
```

Switch between environments with:

```bash
conda activate firecrown_clp
# or
conda activate txpipe_clp
```

## Pipeline Structure

This repository connects three DESC tools through [ceci](https://github.com/LSSTDESC/ceci):

- **TXPipe** — data computation, produces a SACC file
- **TJPCov** — covariance matrix computation, updates the SACC file
- **Firecrown** — prepares the inputs required for MCMC inference

Running the pipeline requires two configuration files:

1. **`CL_concat.yml`** — the *pipeline* config. Defines which stages run,
   their inputs/outputs, parallelization, and other global settings (e.g.
   whether MPI is used).
2. **`config.yml`** — the *stage* config. Contains parameters specific to
   each stage; any modeling change should be made here.

TXPipe, TJPCov, and Firecrown each need their own stage configuration, since
they produce different outputs and run in different Conda environments.

**Pipeline flow:**

```
TXPipe    → SACC file (data vectors)
TJPCov    → SACC file + covariance matrix
Firecrown → inference inputs
```

A separate job must be run to execute the final Firecrown inference — it is
not launched automatically as part of the ceci pipeline.

> **Note:** due to missing components in TJPCov, the Firecrown stage
> currently performs additional computations to rescale the covariance
> matrix.

### Running a single stage

Every stage can also be run standalone, without a full pipeline file, since
each `PipelineStage` subclass exposes its own CLI:

```bash
python -m clpipe <StageName> --config=<config.yml> --<tag>=<path> ...
```

This is mainly useful for debugging a single stage in isolation.

## Example Runs

The `examples` directory contains sample pipeline executions:

- `examples/cosmodc2_remapper` — multiple example analyses
- `examples/cosmodc_halos` — a manually constructed pipeline, useful as a
  reference if you'd rather not use the template
- `examples/template` — scripts to generate the pipeline and stage
  configuration files (currently built for CosmoDC2; support for other
  analyses is under construction)

After generating the necessary files, run the pipeline with ceci, e.g.:

```bash
ceci examples/cosmodc2_remapper/baseline/run_in2p3_mor
```

### Downloading example data

Large chains, sacc files, and (for capish) mock realizations are not
committed to git -- the notebooks under `examples/cosmodc2_redmapper/` and
`examples/capish_simulation/` load them from `CLPipe/data/`, hosted on the
NERSC Science Gateway portal:

https://portal.nersc.gov/cfs/lsst/clpipe/data/

Download the whole tree into the right place with:

```bash
wget -r -np -nH --cut-dirs=3 -R "index.html*" \
  https://portal.nersc.gov/cfs/lsst/clpipe/data/ \
  -P CLPipe/
```

This lands files at `CLPipe/data/cosmodc2_redmapper/{chains,sacc}/` and
`CLPipe/data/capish_simulation/{chains,*.sacc,data_generation/mocks_seeds}/`,
matching the relative paths the notebooks expect. See
[`data/README.md`](data/README.md) for what is hosted there and how to
refresh it.

## Testing

Tests live under `tests/` and run via `pytest`, split into two tiers:

- **Fast** — only needs `ceci` installed. Checks file generation and CLI
  wiring, no real computation.

  ```bash
  pytest tests -m "not slow"
  ```

- **Full** — needs the full DESC stack (pyccl, tjpcov, crow, firecrown,
  cosmosis). Runs real covariance computation and cosmosis sampler runs.

  ```bash
  pytest tests
  ```

CI (`.github/workflows/ci.yml`) runs the fast tier first and only starts the
full-stack tier once it passes.

## Contributing

- Open a pull request against `main`; CI runs automatically (skipped for
  draft PRs).
- The fast tier (`pytest tests -m "not slow"`) gates the full-stack tier —
  if your PR fails fast, the full-stack job never starts.
- Add or update tests under `tests/` for any change to stage behavior,
  generated config files, or the Conda environments.
- For questions or to get involved with the Cluster Working Group, reach
  out to the maintainers listed in `pyproject.toml`.