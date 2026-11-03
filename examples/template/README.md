### This script was developed to simplify the configuration of the analysis pipeline's files and ensure their consistency.

It consists in defining a set of configuration parameters (in a `configuration.yml` file) that can be used to parametrize any file used in the analysis pipeline (config files, job submission script, TXPipe, TJPCov or Firecrown files).

The example in `template_readme_example` generates the runs of `examples/cosmoDC2-1deg2`, for CC-IN2P3 (`run_in2p3`) and NERSC (`run_nersc`). The templates used for the cosmoDC2 redMaPPer analysis are in `examples/cosmodc2_redmapper/templates`.


### >> `configuration.yml` file

Each key defined in the `configuration.yml` file will be initialized as a parameter that can be used in the definition of any file.

For example :

	id: configure

	site: in2p3
	runidentifier: 1deg2

	name: ../cosmoDC2-$RUNIDENTIFIER/run_$SITE
	templatedir: template_readme_example
	computationdir: $PWD/$NAME

	jobconcatyamlfile: $TEMPLATEDIR/job_concat.yml

	datadir: example/inputs
	fiducial_cosmology: ../../../data/fiducial_cosmology_cosmodc2.yml

	txpipe_datadir:
	    nersc: /global/cfs/projectdirs/lsst/groups/CL/cl_pipeline_project/TXPipe_data
	    in2p3: /sps/lsst/groups/clusters/cl_pipeline_project/TXPipe_data
	...
	batch_id: cosmodc2-$RUNIDENTIFIER
	batch_param_file: ./$TEMPLATEDIR/batch_param.yml
	batch_files: ./$TEMPLATEDIR/launch_job_all_generic.sh

will setup a set of parameters that can be used in any file used to prepare the run.<br>
For each parameter defined in `configuration.yml`, a variable $PARAMETER_NAME is defined :

	$SITE = in2p3
	$RUNIDENTIFIER = 1deg2
	$NAME = ../cosmoDC2-1deg2/run_in2p3
	$TEMPLATEDIR = template_readme_example
	...
	$TXPIPE_DATADIR = /sps/lsst/groups/clusters/cl_pipeline_project/TXPipe_data   (dictionaries are read following the $SITE value)
	...

Mandatory parameters :
- `site` : selects the value of every parameter defined per site (here `in2p3` or `nersc`)
- `name` : directory where the configured files are written, relative to where `configure_pipeline.py` is launched. Jobs are submitted from this directory and the outputs and logs are stored in it. Use a different name for each run.
- `templatedir` : directory that contains the template files
- `computationdir` : absolute path of the `name` directory

**It is possible to add/remove parameters in `configuration.yml` depending on the type of pipeline to configure.**
The values are substituted with `sed`, so they must not contain the characters `&` or `'`.


### >> ceci yaml files

The yaml files used to define the different ceci pipeline processes are concatenated in an unique file using the yaml multiple document feature (files are separated using a `---`).<br>
The name of this concatenated file is defined as *jobconcatyamlfile* in `configuration.yml`.<br>
This yaml file is parametrized based on the values defined in `configuration.yml` and then split into smaller yaml files based on the "id" defined in each section.<br>

In our example, the file `template_readme_example/job_concat.yml` is split into 3 files : `TXPipe.yml`, `TJPCov.yml` and `Firecrown.yml`.

The stage configuration file `template_readme_example/cosmodc2_config_1deg2.yml` (`config_file_generic`) is configured the same way and written as `config_<site>.yml` (`config_file`).


### >> `batch_param.yml` file

Instead of redefining the batch parameters (CPU, RAM, tasks, ...) for each job, one can store them in a `batch_param.yml` file and retrieve and configure the batch scripts automatically.
For example :

	id : cosmodc2-1deg2

	# Batch parameters
	batch_all:
	    nersc: -A m1727;-C cpu;--qos=regular;--time=04:00:00;--nodes=1;--ntasks-per-node=10
	    in2p3: --time=4:00:00;--partition=hpc,lsst;--ntasks=10;--cpus-per-task=1;--mem=16gb

will define a `$BATCH_ALL` value, with one `#SBATCH` line per option, used to configure the job submission script `template_readme_example/launch_job_all_generic.sh`.
The configured script is written as `launch_job_all_<site>.sh`.


### >> creation and initialization of the run files

#### --- configuration command

Run from this directory :

	python configure_pipeline.py template_readme_example/configuration.yml


#### --- where are the output files

The configured files are all stored in the `$NAME` directory defined in the `configuration.yml` file, here `../cosmoDC2-1deg2/run_in2p3`.<br>
The job is then submitted from that directory :

	cd ../cosmoDC2-1deg2/run_in2p3
	sbatch launch_job_all_in2p3.sh


#### --- how to overwrite a value defined in the `configuration.yml` file

It is possible to overwrite parameters defined in `configuration.yml` by adding a dictionary of values at the end of the command line.<br>
For example :

	python configure_pipeline.py template_readme_example/configuration.yml '{"site":"nersc"}'

will overwrite the "site" parameter defined in `configuration.yml` and write the NERSC run in `../cosmoDC2-1deg2/run_nersc`.
