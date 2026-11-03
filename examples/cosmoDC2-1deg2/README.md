# cosmoDC2 1 deg² example

Full CLPipe run on the 1 deg² cosmoDC2 sample used by the TXPipe tests:
TXPipe (source selection, BPZ photo-z, cluster binning, stacked ΔΣ profiles) →
TJPCov (covariance) → Firecrown (likelihood) → CosmoSIS.
It is small enough to run in a few minutes before the sampling.

The stage configuration follows the cosmoDC2 redMaPPer analysis (`examples/cosmodc2_redmapper`),
with a binning adapted to the 75 clusters of this sample.

## Directories

`run_in2p3` and `run_nersc` are the runs at CC-IN2P3 and NERSC, submitted with
`sbatch launch_job_all_<site>.sh` from the directory.

Both are generated from `examples/template/template_readme_example`. To change the configuration,
edit the template and regenerate the runs from `examples/template` :

	python configure_pipeline.py template_readme_example/configuration.yml '{"site":"in2p3"}'
	python configure_pipeline.py template_readme_example/configuration.yml '{"site":"nersc"}'

## Data

The inputs come from the TXPipe example data :

	wget https://portal.nersc.gov/cfs/lsst/txpipe/data/example-v10.tar.gz

The files used are `data/example/inputs/metadetect_shear_catalog.hdf5`, `mock_spectroscopic_catalog.hdf5`,
`cluster_catalog.hdf5` and the BPZ templates in `data/example/rail-bpz-inputs`.
The `data/example` directory of the archive goes in the `TXPipe_data` directory of the
`cl_pipeline_project` area, i.e. `TXPipe_data/example/inputs/...`.
The metadetect catalog must have the `shear/ns` group read by TXPipe (the `v10` file has it, older copies only have `shear/00`).
